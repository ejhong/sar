/**
 * Terms from two communities that rarely read each other: seismology and
 * synthetic-aperture radar. Short, plain definitions; the Method page lists
 * them all, and <Term> shows one where it is used.
 */
export interface Entry {
  term: string;
  field: 'seismic' | 'radar' | 'both';
  text: string;
}

export const GLOSSARY: Record<string, Entry> = {
  'p-wave': { term: 'P wave', field: 'seismic', text: 'The compressional wave: the fastest, so the first to arrive. About 3,300 m/s in Giza limestone.' },
  's-wave': { term: 'S wave', field: 'seismic', text: 'The shear wave: slower than P (about 1,800 m/s in the same rock) and unable to cross air or water.' },
  rayleigh: { term: 'Rayleigh wave', field: 'seismic', text: 'A wave that rolls along the ground surface, slightly slower than S, carrying most of a hammer blow’s energy. It samples to about a third of its wavelength in depth.' },
  'first-break': { term: 'first break', field: 'seismic', text: 'The first arrival of energy on a seismic record. Its timing, picked on every record, is the data of travel-time tomography.' },
  traveltime: { term: 'travel-time tomography', field: 'seismic', text: 'Finding the speed of the ground from how long first arrivals take between many sources and receivers, by tracing their rays and adjusting the speeds until the times fit.' },
  eikonal: { term: 'eikonal equation', field: 'seismic', text: 'The high-frequency limit of the wave equation: it gives the earliest arrival time everywhere, and the rays are its steepest descent.' },
  crosshole: { term: 'crosshole', field: 'seismic', text: 'A survey with sources down one borehole and receivers down another, so rays pass through the ground between them rather than skimming the surface.' },
  geophone: { term: 'geophone', field: 'seismic', text: 'A ground-motion sensor recording velocity; three-component geophones record east, north and up.' },
  fwi: { term: 'full-waveform inversion', field: 'seismic', text: 'Fitting whole seismic records, not just arrival times, by adjusting a model of the ground; each step is a gradient computed by propagating the misfit backwards in time.' },
  rtm: { term: 'reverse-time migration', field: 'seismic', text: 'Imaging by sending recorded waves backwards through a model of the ground and marking where they meet the waves sent forward from the source. It is the first step of full-waveform inversion.' },
  'ambient-noise': { term: 'ambient-noise interferometry', field: 'seismic', text: 'Cross-correlating the background hum recorded at two sensors turns it into the wave that would travel from one to the other, so noise can stand in for a source.' },
  microseism: { term: 'microseism', field: 'seismic', text: 'The ground’s constant hum from ocean swell, strongest between 0.1 and 0.3 Hz, felt everywhere on Earth as surface waves kilometres long.' },
  'virtual-source': { term: 'virtual source', field: 'seismic', text: 'A sensor made to act as a source: the correlation of its noise record with another sensor’s is the wave a real source there would have sent. Noise sources spread over the ground give it as minus the lag derivative of the correlation.' },
  pml: { term: 'PML', field: 'seismic', text: 'Perfectly matched layer: an absorbing boundary that lets simulated waves leave the grid without reflecting.' },
  ppw: { term: 'points per wavelength', field: 'seismic', text: 'How many grid cells span the shortest wave simulated; too few and the waves travel at the wrong speed.' },
  tomogram: { term: 'tomogram', field: 'both', text: 'A three-dimensional image of a property (here, wave speed) reconstructed from measurements taken around it.' },
  slc: { term: 'SLC', field: 'radar', text: 'Single-look complex image: a focused radar image that keeps amplitude and phase for every pixel.' },
  dwell: { term: 'Spotlight Dwell', field: 'radar', text: 'An ICEYE mode that keeps the beam on one scene for tens of seconds, here 24.5 s, far longer than a normal image needs.' },
  aperture: { term: 'synthetic aperture', field: 'radar', text: 'The stretch of orbit over which echoes are combined into one image; a long aperture gives fine resolution along the track.' },
  'sub-aperture': { term: 'sub-aperture', field: 'radar', text: 'A part of the synthetic aperture processed on its own: a look from a slightly different direction and time.' },
  doppler: { term: 'Doppler frequency', field: 'radar', text: 'The frequency shift of an echo from the satellite’s motion. In a focused image it maps to the time and direction from which a target was seen.' },
  azimuth: { term: 'azimuth', field: 'radar', text: 'The direction along the satellite’s track, as opposed to range, the direction of the beam.' },
  coherence: { term: 'coherence', field: 'radar', text: 'How similar two radar looks of the same ground are, from 0 to 1. Phase can only be compared where coherence is high.' },
  los: { term: 'line of sight', field: 'radar', text: 'The direction from the satellite to the ground; radar measures motion only along it.' },
  xband: { term: 'X-band', field: 'radar', text: 'Radar near 9.6 GHz (3.1 cm wavelength). It enters dry limestone by centimetres, not metres.' },
  'doppler-tomography': { term: 'Doppler tomography', field: 'radar', text: 'The claimed method: reading ground vibration from the differences between sub-apertures of one image and focusing it into depth. Published 2022, retracted 2026.' },
  look: { term: 'look', field: 'radar', text: 'An image formed from part of the Doppler band. In a dwell the band is time, so a look sees the ground as it was during a few seconds, and averages any motion over them.' },
  speckle: { term: 'speckle', field: 'radar', text: 'The grainy texture of radar images of natural ground: many small scatterers in each cell interfering. Each Doppler bin holds an independent speckle pattern.' },
  fisher: { term: 'Fisher information', field: 'both', text: 'How sharply the probability of the data changes with a parameter. Its inverse bounds the variance of any unbiased estimate (the Cramér–Rao bound), and for a small change it fixes how well any test can tell the two cases apart.' },
  imprint: { term: 'imprint', field: 'both', text: 'The difference a buried chamber makes to how the surface moves as a wave passes: what any passive method would have to detect.' },
};
