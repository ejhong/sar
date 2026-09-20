"""Portable field-atlas publishing with provenance and model-integrity checks."""
import hashlib
import html
import json
from pathlib import Path
import shutil
from string import Template

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def confined(root, relative):
    path = (root/relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'Asset must remain within {root.name}: {relative}')
    return path


def validate_model(model):
    if model['native_target_mask'] is not None or model['field_detection'] is not None:
        raise ValueError('Historical reference model cannot contain a radar result')
    if not model['assumptions'] or not model['datum'] or not model['source_url'].startswith('https://'):
        raise ValueError('A survey model needs a source, datum and disclosed simplifications')
    depth = model['depth_m']
    if not np.isfinite(depth) or depth <= 0:
        raise ValueError('Model depth must be finite and positive')
    for part in model['parts']:
        polygon = np.asarray(part['polygon'], float)
        top, bottom = part['top_m'], part['bottom_m']
        if (polygon.ndim != 2 or polygon.shape[1] != 2 or len(polygon) < 3
                or not np.isfinite(polygon).all() or not np.isfinite([top,bottom]).all()
                or not -depth-1e-9 <= bottom < top <= 0):
            raise ValueError('Invalid metric survey volume')
        area = abs(np.sum(polygon[:,0]*np.roll(polygon[:,1],-1)-polygon[:,1]*np.roll(polygon[:,0],-1)))/2
        if area <= 0:
            raise ValueError('Empty survey footprint')


def obj_model(model):
    """Independent prism envelopes; intentional junction overlap, no boolean union."""
    validate_model(model)
    lines = ['# Historical survey reference, NOT a SAR reconstruction',
             '# Units: metres; axes x=east y=north z=up; separate simplified volumes',
             '# Not georeferenced, not a watertight engineering model',
             '# '+model['source_url'], '# '+model['datum']]
    lines += ['# Simplification: '+a for a in model['assumptions']]
    offset = 1
    for i, part in enumerate(model['parts']):
        n = len(part['polygon'])
        lines.append(f'o {model["id"]}_{i}_{part["kind"]}')
        for z in (part['top_m'], part['bottom_m']):
            lines += [f'v {x:.6f} {y:.6f} {z:.6f}' for x, y in part['polygon']]
        faces = [list(range(n)), list(reversed(range(n,2*n)))]
        faces += [[j,(j+1)%n,(j+1)%n+n,j+n] for j in range(n)]
        lines += ['f '+' '.join(str(offset+j) for j in face) for face in faces]
        offset += 2*n
    return '\n'.join(lines)+'\n'


def build_field_site(docs):
    assets = ROOT/'site/assets'
    for source, destination in [('main.css','style.css'),('main.js','site.js'),
                                ('field.css','field.css'),('model-viewer.js','model-viewer.js')]:
        shutil.copyfile(assets/source, docs/destination)
    js = (assets/'field.js').read_text().replace("'./model-viewer.js'",f"'./model-viewer.js?v={sha(assets/'model-viewer.js')[:12]}'")
    (docs/'field.js').write_text(js)
    catalog = read(ROOT/'catalog/sites.json')
    site_ids = [site['id'] for site in catalog['sites']]
    acquisitions = [a for site in catalog['sites'] for a in site['acquisitions']]
    for values in (site_ids, [a['id'] for a in acquisitions], [a['report'] for a in acquisitions]):
        if len(values) != len(set(values)):
            raise ValueError('Site IDs, acquisition IDs and report paths must be unique')
    out = docs/'data/field'
    out.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(ROOT/'catalog/sites.json',out/'catalog.json')
    cards = []
    for site in catalog['sites']:
        survey_file = confined(ROOT/'catalog',site['surveys'])
        surveys = read(survey_file)
        destination = confined(out,site['surveys'])
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(survey_file,destination)
        for model in surveys['models']:
            destination.with_name(model['id']+'.obj').write_text(obj_model(model))
        for acquisition in site['acquisitions']:
            audit_file=confined(ROOT/'research',acquisition['geometry_audit'])
            audit=read(audit_file)
            for path,digest in audit['source_sha256'].items():
                if sha(confined(ROOT,path)) != digest:
                    raise ValueError(f'Acquisition geometry audit is stale: {path}')
            audit_destination=confined(docs/'data',acquisition['geometry_audit'])
            audit_destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(audit_file,audit_destination)
            report_file = confined(ROOT/'results/field',acquisition['report'])
            report = read(report_file)
            if report['acquisition'] != acquisition or report['site'] != site['id']:
                raise ValueError('Catalog and measured acquisition report disagree; rerun field checks')
            for path,digest in report['source_sha256'].items():
                if sha(confined(ROOT,path)) != digest:
                    raise ValueError(f'Field-check source changed: {path}; rerun scripts/run_field_checks.py')
            design_file = confined(ROOT/'catalog',acquisition['design'])
            if sha(design_file) != report['design_sha256']:
                raise ValueError('Field design changed after measurement')
            for path,digest in report['output_sha256'].items():
                source = confined(report_file.parent,path)
                if sha(source) != digest:
                    raise ValueError(f'Derived field artifact changed: {path}')
                dest = confined(out,acquisition['report']).parent/path
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(source,dest)
            dest_report = confined(out,acquisition['report'])
            dest_report.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(report_file,dest_report)
            dest_design = confined(out,acquisition['design'])
            dest_design.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(design_file,dest_design)
        # The first acquisition is the stable cover; all are selectable in the atlas.
        acquisition=site['acquisitions'][0]
        report=read(ROOT/'results/field'/acquisition['report'])
        image=report['views'][0]
        image_path=(Path('data/field')/acquisition['report']).parent/image['image']
        e=html.escape
        cards.append(f'<a class="field-card" href="field.html?site={e(site["id"])}"><img src="{image_path}" width="{image["width"]}" height="{image["height"]}" alt="Actual SAR intensity of {e(site["name"])}; surface returns and radar shadows" loading="lazy"><div><span class="eyebrow">{e(site["country"])} · {e(acquisition["date"])}</span><h3>{e(site["name"])} <span aria-hidden="true">↗</span></h3><p>{len(surveys["models"])} dimensioned survey models · real-image checks</p></div></a>')
    template=Template((ROOT/'site/field.html').read_text())
    (docs/'field.html').write_text(template.substitute(
        stylesheet_version=sha(docs/'style.css')[:12],
        field_stylesheet_version=sha(docs/'field.css')[:12],
        field_script_version=sha(docs/'field.js')[:12]))
    return '<div class="field-atlas-cards">'+''.join(cards)+'</div>'
