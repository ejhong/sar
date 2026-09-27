"""M1-01 · How much does the ground near Giza move on its own? Measured.

    uv run python experiments/m1_01_ambient_levels.py

The nearest open broadband record to Giza is MedNet station KEG at Kottamya
(29.9275 N, 31.8292 E, 460 m, a Streckeisen STS-1 at 7 m depth; about 67 km
east of Khafre), archived openly at INGV. One day in the middle of every
month of 1996, three components at 20 samples per second, is corrected for
the instrument response and turned into power spectral densities the way
network operators do (McNamara & Buland 2004, as implemented in ObsPy's
PPSD: one-hour windows overlapping by half). The percentiles of those
spectra, compared with Peterson's (1993) low- and high-noise models, give
the ambient ground motion in each band, measured rather than assumed.

Kottamya is an observatory on a quiet limestone hill: its microseisms
(below 1 Hz) are those of the whole region, Giza included; above 1 Hz it is
a lower bound for Giza, which borders a city of twenty million.
"""
import io
import time
import urllib.request

import numpy as np

from katabasis.ambient import peterson
from katabasis.runs import RESULTS, Run

FDSN = 'https://webservices.ingv.it/fdsnws'
NET, STA = 'MN', 'KEG'
DAYS = [f'1996-{m:02d}-15' for m in range(1, 13)]
DATA = RESULTS.parent / 'data' / 'kottamya'
BANDS = [(0.05, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 8.0)]
PCTS = [10, 50, 90]
UTC_OFFSET_H = 2           # Egypt standard time


def fetch(url: str, path) -> bytes:
    if path.exists():
        return path.read_bytes()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                data = r.read()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return data
        except Exception as e:                      # noqa: BLE001
            print(f'  retry {attempt + 1}: {e}', flush=True)
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f'could not fetch {url}')


def band_rms(freq: np.ndarray, acc_db: np.ndarray, lo: float, hi: float) -> float:
    """RMS velocity (m/s) in a band from an acceleration PSD in dB."""
    m = (freq >= lo) & (freq <= hi) & np.isfinite(acc_db)
    f = freq[m][np.argsort(freq[m])]
    pv = 10 ** (acc_db[m][np.argsort(freq[m])] / 10) / (2 * np.pi * f) ** 2
    return float(np.sqrt(np.trapezoid(pv, f)))


def main():
    from obspy import UTCDateTime, read, read_inventory
    from obspy.signal import PPSD

    params = {'station': f'{NET}.{STA}', 'archive': FDSN, 'days': DAYS, 'channels': 'BHZ, BHN, BHE (20 sps)',
              'method': 'McNamara & Buland (2004) via obspy PPSD, 1 h windows, 50% overlap',
              'local_time': f'UTC+{UTC_OFFSET_H}', 'bands_hz': BANDS, 'percentiles': PCTS}
    with Run('m1_01_ambient_levels', 'Ambient ground motion measured near Giza', params) as run:
        inv = read_inventory(io.BytesIO(fetch(f'{FDSN}/station/1/query?net={NET}&sta={STA}&level=response',
                                              DATA / 'response.xml')), format='STATIONXML')
        out = {'components': {}}
        # local day 08-18 and night 22-04, in UTC hours for PPSD's time-of-day selection
        day = [(-1, 8 - UTC_OFFSET_H, 18 - UTC_OFFSET_H)]
        night = [(-1, 22 - UTC_OFFSET_H, 24), (-1, 0, 4 - UTC_OFFSET_H)]
        for cha in ('BHZ', 'BHN', 'BHE'):
            pp = None
            for d in DAYS:
                t0 = UTCDateTime(d)
                url = (f'{FDSN}/dataselect/1/query?net={NET}&sta={STA}&cha={cha}'
                       f'&start={t0.isoformat()}&end={(t0 + 86400).isoformat()}')
                raw = fetch(url, DATA / f'{cha}_{d}.mseed')
                if not raw:
                    print(f'  {cha} {d}: no data', flush=True)
                    continue
                st = read(io.BytesIO(raw))
                if pp is None:
                    pp = PPSD(st[0].stats, metadata=inv, ppsd_length=3600, overlap=0.5)
                pp.add(st)
                print(f'  {cha} {d}: {len(st)} traces', flush=True)
            comp = {'windows': len(pp.times_processed), 'hours_of_record': round(len(pp.times_processed) / 2, 1)}
            per, curves = None, {}
            for key, sel in (('all', None), ('day', day), ('night', night)):
                pp.calculate_histogram(time_of_weekday=sel)
                for q in PCTS:
                    per, db = pp.get_percentile(percentile=q)
                    curves[f'{key}_p{q}'] = np.asarray(db)
                comp[f'windows_{key}'] = int(pp.current_histogram_count)
            per = np.asarray(per)
            freq = 1 / per
            comp['period_s'] = per.round(4)
            comp['acc_db'] = {k: v.round(2) for k, v in curves.items()}
            comp['band_rms_um_s'] = {
                f'{lo}-{hi} Hz': {k: band_rms(freq, v, lo, hi) * 1e6 for k, v in curves.items()} for lo, hi in BANDS}
            out['components'][cha] = comp
            print(f"  {cha}: {comp['windows']} windows; median rms "
                  + ', '.join(f"{b} {v['all_p50']:.3g}" for b, v in comp['band_rms_um_s'].items()) + ' um/s', flush=True)
        out['peterson'] = {
            'source': 'Peterson (1993), USGS OFR 93-322, Tables 3-4',
            'band_rms_um_s': {f'{lo}-{hi} Hz': {'low': peterson.band_rms_velocity(lo, hi, 'low') * 1e6,
                                                 'high': peterson.band_rms_velocity(lo, hi, 'high') * 1e6}
                              for lo, hi in BANDS},
            'period_s': np.geomspace(0.1, 30, 120).round(4),
        }
        out['peterson']['low_db'] = peterson.acceleration_db(out['peterson']['period_s'], 'low').round(2)
        out['peterson']['high_db'] = peterson.acceleration_db(out['peterson']['period_s'], 'high').round(2)
        z = out['components']['BHZ']['band_rms_um_s']
        micro = z['0.1-0.3 Hz']['all_p50']
        cult_d, cult_n = z['1.0-3.0 Hz']['day_p50'], z['1.0-3.0 Hz']['night_p50']
        out['finding'] = (f"At Kottamya, 67 km east of Giza, the median vertical ground velocity is {micro:.2g} um/s in the "
                          f"microseism band (0.1-0.3 Hz) and {cult_d:.2g} um/s by day, {cult_n:.2g} um/s by night at 1-3 Hz "
                          f"(1996, twelve days, instrument response removed). Peterson's global models bound the "
                          f"microseism band between {out['peterson']['band_rms_um_s']['0.1-0.3 Hz']['low']:.2g} and "
                          f"{out['peterson']['band_rms_um_s']['0.1-0.3 Hz']['high']:.2g} um/s.")
        run.save(out)


if __name__ == '__main__':
    main()
