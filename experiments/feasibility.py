"""Supporting experiments adapted from the user-supplied SAR-Voids-Quick-Experiments.zip.

Reproducible synthetic feasibility probes; NO real SAR or field void data.

Run: python run_experiments.py
Outputs: JSON/CSV metrics, NPZ waveforms, and PNG figures in results/.
These are information-channel demonstrations, not estimates of satellite
underground detection performance. See README.md for assumptions and limits.
"""
from pathlib import Path
import csv
import json
import platform
import sys
import time
import numpy as np
import scipy
from scipy.interpolate import RegularGridInterpolator
import sklearn
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sarsim import viz

OUT = ROOT / 'results' / 'feasibility'
OUT.mkdir(parents=True, exist_ok=True)
viz.style()
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})


def save_csv(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def motion_channel():
    """Ideal dominant scatterer after perfect platform-motion removal.

    z[n] = exp(i * 4*pi*u[n]/lambda) + circular complex Gaussian noise.
    This is a pulse-time receiver model, NOT a focused SLC image simulator.
    Intensity comparison is pulse intensity, not formed-image texture.
    """
    rng = np.random.default_rng(20260919)
    n, fs, wavelength, snr_db = 4096, 512., .031, 20.
    t = np.arange(n) / fs
    f = np.fft.rfftfreq(n, 1 / fs)
    band = (f >= 1) & (f <= 30)
    sigma = 10 ** (-snr_db / 20)

    def measure(amplitude_um, count):
        phi = rng.uniform(-np.pi, np.pi, (count, 1))
        u = amplitude_um * 1e-6 * np.sin(2*np.pi*7.5*t + phi)
        z = np.exp(1j*4*np.pi*u/wavelength)
        z += sigma / np.sqrt(2) * (rng.normal(size=z.shape) + 1j*rng.normal(size=z.shape))
        phase = np.angle(z)
        power = np.abs(z)**2
        sp = np.abs(np.fft.rfft(phase-phase.mean(axis=1,keepdims=True),axis=1))[:,band]*2/n
        si = np.abs(np.fft.rfft(power-power.mean(axis=1,keepdims=True),axis=1))[:,band]*2/n
        return sp.max(axis=1), si.max(axis=1), f[band][sp.argmax(axis=1)]

    null_phase, null_intensity, _ = measure(0., 2000)
    thp, thi = np.quantile(null_phase, .99), np.quantile(null_intensity, .99)
    rows=[]
    for a in [0., 1., 3., 10., 30., 100., 300.]:
        p, intensity, freq = measure(a, 800)
        rows.append({'amplitude_um':a, 'phase_detection_rate':float(np.mean(p>thp)),
                     'intensity_detection_rate':float(np.mean(intensity>thi)),
                     'phase_detection_and_correct_frequency_rate':float(np.mean((p>thp)&(np.abs(freq-7.5)<=fs/n)))})
    save_csv('motion_detection.csv',rows)
    fig,ax=plt.subplots(figsize=(7.5,4.2),layout='constrained')
    rr=rows[1:]
    ax.semilogx([r['amplitude_um'] for r in rr],[r['phase_detection_rate'] for r in rr],'o-',label='Complex phase')
    ax.semilogx([r['amplitude_um'] for r in rr],[r['intensity_detection_rate'] for r in rr],'s--',label='Pulse intensity')
    ax.axhline(.01,color='gray',lw=1,ls=':',label='Nominal 1% false alarm')
    ax.set(xlabel='Injected surface displacement amplitude (micrometres)',ylabel='Detection fraction',ylim=(-.04,1.05),
           title='Synthetic receiver: phase retains tiny surface motion')
    ax.legend(loc='upper left')
    fig.savefig(OUT/'motion_detection.png');plt.close(fig)
    return {'model':'single ideal coherent scatterer; not a focused SAR image', 'wavelength_m':wavelength,
            'duration_s':n/fs,'pulse_rate_hz':fs,'pulse_snr_db':snr_db,'frequency_hz':7.5,
            'noise':'independent circular complex Gaussian; perfect platform compensation',
            'threshold_calibration_trials':2000,'test_trials_per_amplitude':800,'rows':rows}


def sh_wave(depth=None, f0=25., dx=5., c=2000., duration=.70):
    """2-D homogeneous anti-plane SH elastic-wave surrogate.

    u_tt = c^2 div(grad u) in rock; zero normal derivative at surface/cavity.
    A circular hole in x-z represents an INFINITE cylindrical void along y.
    Boundary faces adjacent to void have zero shear traction (zero flux).
    Side and bottom sponge, no intrinsic attenuation, no heterogeneity.
    Source is a prescribed Ricker body force at (-150 m, 10 m).
    """
    x=np.arange(-800.,800.+dx/2,dx)
    z=np.arange(0.,800.+dx/2,dx)
    dt=.4*dx/c  # c*dt/dx=.4 < 1/sqrt(2)
    t=np.arange(int(duration/dt))*dt
    X,Z=np.meshgrid(x,z)
    rock=np.ones_like(X,dtype=bool)
    if depth is not None: rock=(X**2+(Z-depth)**2)>20.**2
    mx=rock[:,1:] & rock[:,:-1]
    mz=rock[1:,:] & rock[:-1,:]
    # Sponge thickness fixed in metres, rate fixed per second for grid comparison.
    edge=np.maximum(np.maximum((np.abs(X)-650)/150,0),(Z-650)/150)
    gamma=250.*np.maximum(edge,0)**2
    src=(int(round(10/dx)),int(round((-150-x[0])/dx)))
    rix=np.arange(int(round((-300-x[0])/dx)),int(round((300-x[0])/dx))+1)
    old=np.zeros_like(X);u=np.zeros_like(X)
    traces=np.zeros((len(t),len(rix)))
    maxabs=0.
    for k,tk in enumerate(t):
        lap=np.zeros_like(u)
        fx=(u[:,1:]-u[:,:-1])*mx
        fz=(u[1:,:]-u[:-1,:])*mz
        lap[:,:-1]+=fx;lap[:,1:]-=fx
        lap[:-1,:]+=fz;lap[1:,:]-=fz
        new=(2*u-(1-gamma*dt/2)*old+(c*dt/dx)**2*lap)/(1+gamma*dt/2)
        a=np.pi*f0*(tk-1.5/f0)
        new[src]+=dt**2/dx**2*(1-2*a*a)*np.exp(-a*a)
        new*=rock
        old,u=u,new
        traces[k]=u[0,rix]
        maxabs=max(maxabs,float(np.max(np.abs(u))))
    assert np.isfinite(traces).all() and maxabs<1, 'unstable field'
    return t,x[rix],traces


def wave_channel():
    rows=[]; arrays={}; plot_data=[]
    for f0 in [8.,25.]:
        t,x,base=sh_wave(f0=f0)
        norm=float(np.linalg.norm(base)); peak=float(np.abs(base).max())
        arrays[f't_{int(f0)}']=t;arrays[f'x_{int(f0)}']=x
        arrays[f'baseline_{int(f0)}']=base
        for d in [60.,120.,240.]:
            _,_,v=sh_wave(depth=d,f0=f0)
            diff=v-base
            arrays[f'cavity_f{int(f0)}_d{int(d)}']=v
            rows.append({'frequency_hz':f0,'depth_m':d,'radius_m':20.,
                         'relative_waveform_L2_change':float(np.linalg.norm(diff)/norm),
                         'peak_difference_over_peak_baseline':float(np.abs(diff).max()/peak)})
            if f0==25.:plot_data.append((d,diff/peak))
    # One concrete numerical check: compare baseline and cavity perturbation on finer grid.
    tf,xf,bf=sh_wave(f0=25.,dx=2.5)
    _,_,vf=sh_wave(depth=120.,f0=25.,dx=2.5)
    tc,xc=arrays['t_25'],arrays['x_25']
    points=np.stack(np.meshgrid(tc,xc,indexing='ij'),axis=-1)
    bfi=RegularGridInterpolator((tf,xf),bf)(points)
    dfi=RegularGridInterpolator((tf,xf),vf-bf)(points)
    coarse=arrays['baseline_25']
    diffc=arrays['cavity_f25_d120']-coarse
    qa={'baseline_waveform_correlation_coarse_fine':float(np.corrcoef(coarse.ravel(),bfi.ravel())[0,1]),
        'cavity_difference_correlation_coarse_fine':float(np.corrcoef(diffc.ravel(),dfi.ravel())[0,1]),
        'fine_grid_relative_waveform_L2_change_depth120':float(np.linalg.norm(vf-bf)/np.linalg.norm(bf)),
        'coarse_dx_m':5.,'fine_dx_m':2.5}
    save_csv('wave_sensitivity.csv',rows)
    cache = ROOT / 'results' / 'cache'
    cache.mkdir(exist_ok=True)
    np.savez_compressed(cache/'feasibility_waveforms.npz',**arrays)
    fig,axes=plt.subplots(1,3,figsize=(11,4.4),sharex=True,sharey=True,layout='constrained')
    lim=max(np.abs(v).max() for _,v in plot_data)
    for ax,(d,v) in zip(axes,plot_data):
        im=ax.imshow(v,origin='lower',aspect='auto',extent=[xc[0],xc[-1],tc[0],tc[-1]],cmap='RdBu_r',vmin=-lim,vmax=lim)
        ax.set(title=f'Cavity centre: {d:.0f} m deep',xlabel='Surface receiver position (m)')
    axes[0].set_ylabel('Time (s)')
    fig.colorbar(im,ax=axes,label='Difference / baseline peak displacement',shrink=.85)
    fig.suptitle('Synthetic SH waves: cavity changes the surface response\n20 m radius; prescribed 25 Hz source; no intrinsic attenuation')
    fig.savefig(OUT/'cavity_wave_response.png');plt.close(fig)
    return {'model':'2D SH half-space; infinite cylindrical cavity; known Ricker source',
            'shear_speed_m_s':2000.,'radius_m':20.,'intrinsic_attenuation':False,
            'source_position_m':[-150.,10.],'window_s':.70,
            'rows':rows,'numerical_check':qa,
            'warning':'Relative noiseless changes are not satellite detection probabilities or validated depth limits.'}


def ai_controls():
    """Two deliberately synthetic controls, not trained on real radar imagery.

    Site-only control has random site labels and no transferable void signal.
    Matched-site control has both classes at every site, plus a declared planted
    signal; held-out sites show the benefit of real transferable information.
    """
    rows=[]
    for seed in range(10):
        rng=np.random.default_rng(9000+seed)
        nsites,npatch,dim=80,40,32
        groups=np.repeat(np.arange(nsites),npatch)
        site_y=np.array([0,1]*(nsites//2));rng.shuffle(site_y)
        y=site_y[groups]
        X=3*rng.normal(size=(nsites,dim))[groups]+.5*rng.normal(size=(len(groups),dim))
        tr,te=train_test_split(np.arange(len(y)),test_size=.3,random_state=seed,stratify=y)
        clf=KNeighborsClassifier(1).fit(X[tr],y[tr])
        rows.append({'seed':seed,'experiment':'site_identity_only','split':'random_patches',
                     'signal_shift_sd':0.,'auc':float(roc_auc_score(y[te],clf.predict_proba(X[te])[:,1]))})
        splitter=StratifiedKFold(5,shuffle=True,random_state=seed)
        pred=np.zeros(len(y))
        for sitetr,sitete in splitter.split(np.arange(nsites),site_y):
            tr=np.isin(groups,sitetr);te=np.isin(groups,sitete)
            clf=KNeighborsClassifier(1).fit(X[tr],y[tr])
            pred[te]=clf.predict_proba(X[te])[:,1]
        rows.append({'seed':seed,'experiment':'site_identity_only','split':'held_out_sites',
                     'signal_shift_sd':0.,'auc':float(roc_auc_score(y,pred))})
        # Both classes in each site. Stable signal independent of site nuisance.
        ym=np.tile(np.repeat([0,1],npatch//2),nsites)
        nuisance=3*rng.normal(size=(nsites,8))[groups]+rng.normal(size=(len(groups),8))
        common=rng.normal(size=len(groups))
        for strength in [0.,.3,1.,2.]:
            Xm=np.column_stack([common+strength*ym,nuisance])
            pred=np.zeros(len(ym))
            for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=seed).split(Xm,ym,groups):
                clf=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=1000))
                clf.fit(Xm[tr],ym[tr]);pred[te]=clf.predict_proba(Xm[te])[:,1]
            rows.append({'seed':seed,'experiment':'matched_sites_planted_signal','split':'held_out_sites',
                         'signal_shift_sd':strength,'auc':float(roc_auc_score(ym,pred))})
    save_csv('ai_control_runs.csv',rows)
    summary=[]
    keys=sorted(set((r['experiment'],r['split'],r['signal_shift_sd']) for r in rows))
    for k in keys:
        vals=[r['auc'] for r in rows if (r['experiment'],r['split'],r['signal_shift_sd'])==k]
        summary.append({'experiment':k[0],'split':k[1],'signal_shift_sd':k[2],
                        'mean_auc':float(np.mean(vals)),'sd_across_10_seeds':float(np.std(vals,ddof=1)),
                        'min_auc':float(min(vals)),'max_auc':float(max(vals))})
    save_csv('ai_control_summary.csv',summary)
    fig,axes=plt.subplots(1,2,figsize=(9,4.3),layout='constrained')
    leak=[next(r for r in summary if r['experiment']=='site_identity_only' and r['split']==s) for s in ['random_patches','held_out_sites']]
    axes[0].bar(['Random image patches','Entire sites held out'],[r['mean_auc'] for r in leak],
                yerr=[r['sd_across_10_seeds'] for r in leak],color=['#cc6b48','#347f87'],capsize=4)
    axes[0].set(title='No transferable void signal was supplied',ylabel='ROC AUC',ylim=(0,1.08))
    signal=[r for r in summary if r['experiment']=='matched_sites_planted_signal']
    axes[1].errorbar([r['signal_shift_sd'] for r in signal],[r['mean_auc'] for r in signal],
                     yerr=[r['sd_across_10_seeds'] for r in signal],fmt='o-',capsize=4,color='#347f87')
    axes[1].set(title='With a planted transferable signal',xlabel='Class signal difference / noise SD',ylim=(0,1.08))
    for ax in axes:ax.axhline(.5,color='gray',ls=':',label='Chance');ax.legend(loc='lower right')
    fig.suptitle('Synthetic AI controls: good scores need the right holdout')
    fig.savefig(OUT/'ai_controls.png');plt.close(fig)
    return {'note':'Artificial datasets designed to isolate confounding; NOT measured SAR performance.',
            'sites':80,'patches_per_site':40,'seeds':10,'summary':summary}



def phase_readout():
    """Isolate what phase and translation estimators measure on exactly the same patch.

    A uniform phase rotation is a controlled algebraic input, not a focused-SAR
    vibration simulator. No spatial translation is introduced.
    """
    from sarsim.track import patch_shifts
    rng = np.random.default_rng(73)
    ref = rng.normal(size=(96, 96)) + 1j * rng.normal(size=(96, 96))
    wavelength = 0.031
    amplitude_um = 30.0
    t = np.linspace(0, 1, 41)
    displacement = amplitude_um * np.sin(2 * np.pi * t)
    phase_estimate, az, rg = [], [], []
    for d in displacement:
        off = ref * np.exp(4j * np.pi * d * 1e-6 / wavelength)
        phase_estimate.append(np.angle(np.vdot(ref, off)) * wavelength / (4 * np.pi) * 1e6)
        dr, dc = patch_shifts(ref, off, np.array([32, 48, 64]), np.array([32, 48, 64]))
        az.append(float(np.median(dr)))
        rg.append(float(np.median(dc)))
    rows = [dict(cycle=float(cycle), imposed_phase_equivalent_um=float(d),
                 measured_phase_equivalent_um=float(p), az_shift_px=a, rg_shift_px=r)
            for cycle, d, p, a, r in zip(t, displacement, phase_estimate, az, rg)]
    save_csv('phase_readout.csv', rows)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), layout='constrained')
    axs[0].plot(t, displacement, color='#9a5b33', lw=2.5, label='Imposed phase equivalent')
    axs[0].plot(t, phase_estimate, '.', color='#26756b', ms=6, label='Measured from complex phase')
    axs[0].set(xlabel='Fraction of imposed cycle', ylabel='Phase-equivalent displacement (µm)', title='A phase measurement follows the input')
    axs[0].legend(fontsize=8)
    axs[1].plot(t, az, color='#26756b', label='Azimuth shift')
    axs[1].plot(t, rg, '--', color='#9a5b33', label='Range shift')
    axs[1].set(xlabel='Fraction of imposed cycle', ylabel='Patch translation (pixels)', ylim=(-.001,.001), title='The translation estimate stays at zero')
    axs[1].legend(fontsize=8)
    fig.savefig(OUT/'phase_readout.png'); plt.close(fig)
    error = float(np.max(np.abs(np.array(phase_estimate) - displacement)))
    shift = float(max(np.max(np.abs(az)), np.max(np.abs(rg))))
    assert error < 1e-9 and shift < 1e-9
    return {'amplitude_um': amplitude_um, 'max_phase_error_um': error,
            'max_translation_px': shift, 'steps': len(t),
            'scope': 'Uniform phase rotation only; no spatial motion or SAR focusing.'}

if __name__=='__main__':
    start=time.time()
    all_results={'scope':'SYNTHETIC ONLY. No real satellite scene or surveyed cavity labels were processed.',
                 'versions':{'python':platform.python_version(),'numpy':np.__version__,
                             'scipy':scipy.__version__,'sklearn':sklearn.__version__,'matplotlib':matplotlib.__version__}}
    for name,fn in [('motion',motion_channel),('waves',wave_channel),('ai',ai_controls)]:
        print('Running',name,flush=True)
        all_results[name]=fn()
        print(json.dumps({name:all_results[name]},indent=2),flush=True)
    all_results['elapsed_seconds']=time.time()-start
    all_results['phase_readout'] = phase_readout()
    all_results['date'] = time.strftime('%Y-%m-%d')
    all_results['source'] = 'Adapted from SAR-Voids-Quick-Experiments.zip supplied with this project.'
    (OUT/'metrics.json').write_text(json.dumps(all_results,indent=2,allow_nan=False))
    print('Complete in',round(time.time()-start,1),'seconds',flush=True)
