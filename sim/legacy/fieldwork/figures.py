"""Exportable scientific figures from saved native-image spectrum measurements."""
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np


def spectrum_figure(spectra, timing, output, compact=False):
    style={'font.family':'DejaVu Sans','font.size':9,'axes.labelcolor':'#46554b',
           'text.color':'#46554b','xtick.color':'#46554b','ytick.color':'#46554b',
           'axes.edgecolor':'#bac6b4','svg.fonttype':'none','svg.hashsalt':'sar-field-atlas',
           'figure.facecolor':'#f8f9f4','axes.facecolor':'#f8f9f4'}
    with plt.rc_context(style):
        fig, ax=plt.subplots(figsize=(3.5,3.4) if compact else (9,2.8))
        half=timing['processed_bandwidth_hz']/2000
        ax.axvspan(-half,half,color='#e7ece1',label='Zero-centred metadata band')
        for strip,color in zip(spectra,['#256c5d','#b67c37','#727c97']):
            frequencies=np.array([p['image_frequency_hz'] for p in strip['curve']])/1000
            power=np.array([p['mean_power'] for p in strip['curve']])
            ax.plot(frequencies,10*np.log10(power/power.max()),color=color,lw=1,
                    label=f'Columns {strip["cols"][0]}–{strip["cols"][1]-1}')
        ax.set(xlim=(-timing['processing_prf_hz']/2000,timing['processing_prf_hz']/2000),
               ylim=(-50,2),xlabel='Image Fourier frequency (kHz)',ylabel='Relative power (dB)')
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',color='#ccd6c5',lw=.5)
        if compact:
            ax.set_xticks([-50,0,50])
        handles,labels=ax.get_legend_handles_labels()
        # Keep the diagnostic band's status explicit in the portable artifact.
        ax.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,-.30 if compact else -.27),
                  ncol=1 if compact else 4,frameon=False,fontsize=7)
        fig.subplots_adjust(left=.18 if compact else .09,right=.97,top=.97,bottom=.43 if compact else .31)
        fig.savefig(output,metadata={'Date':None,'Description':'Actual native-row Fourier power, not ground-vibration frequency. Each strip normalized separately. No depth inferred.'})
        plt.close(fig)
        output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
