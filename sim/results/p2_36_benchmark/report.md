# P2-36 benchmark: generated report

Frozen parameters hash `1c58de9ffe7c`; frozen frequency 69 Hz (the presence (room) pair's oracle on the measured disc, largest). SNR 30 dB assumed. States: near chance (TV < 0.05), 95% at 5% excluded (TV < 0.9), unresolved.

## quiet

| pair | L1 single image, TV | state | growth to 0.9 | L1 within 70 m | oracle, TV | state | SNR for 0.9 | oracle within 70 m | oracle within 39 m | achieved (score test) | tail share |
|---|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 2.2e-07 | near chance | 3.6e+06 | 2.1e-07 | 3.3e-06 | near chance | 142 dB | 3.1e-06 | 3.1e-06 | 2.8e-09 | 9% |
| presence (L tunnel) | 2.5e-07 | near chance | 3.1e+06 | 2.4e-07 | 3.8e-06 | near chance | 141 dB | 3.6e-06 | 3.5e-06 | 3.1e-09 | 9% |
| location (6 m) | 1.9e-07 | near chance | 4e+06 | 1.9e-07 | 2.8e-06 | near chance | 143 dB | 2.8e-06 | 2.8e-06 | 3e-09 | 1% |
| shape (room or L tunnel) | 1.8e-07 | near chance | 4.2e+06 | 1.8e-07 | 2.7e-06 | near chance | 144 dB | 2.7e-06 | 2.7e-06 | 2.9e-09 | 1% |

SNR sweep (L1 and oracle, TV):

| pair | 20 dB | 30 dB | 40 dB | 50 dB | 60 dB |
|---|---|---|---|---|---|
| presence (room) | 2.2e-07 / 1e-06 | 2.2e-07 / 3.3e-06 | 2.2e-07 / 1e-05 | 2.2e-07 / 3.3e-05 | 2.2e-07 / 0.0001 |
| presence (L tunnel) | 2.5e-07 / 1.2e-06 | 2.5e-07 / 3.8e-06 | 2.5e-07 / 1.2e-05 | 2.5e-07 / 3.8e-05 | 2.5e-07 / 0.00012 |
| location (6 m) | 1.9e-07 / 8.9e-07 | 1.9e-07 / 2.8e-06 | 1.9e-07 / 8.9e-06 | 1.9e-07 / 2.8e-05 | 1.9e-07 / 8.9e-05 |
| shape (room or L tunnel) | 1.8e-07 / 8.6e-07 | 1.8e-07 / 2.7e-06 | 1.8e-07 / 8.6e-06 | 1.8e-07 / 2.7e-05 | 1.8e-07 / 8.6e-05 |

## strong

| pair | L1 single image, TV | state | growth to 0.9 | L1 within 70 m | oracle, TV | state | SNR for 0.9 | oracle within 70 m | oracle within 39 m | achieved (score test) | tail share |
|---|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 0.033 | near chance | 26.9 | 0.006 | 0.466 | 95% at 5% excluded | 38.4 dB | 0.087 | 0.075 | 0.002 | 97% |
| presence (L tunnel) | 0.062 | 95% at 5% excluded | 14.5 | 0.009 | 0.750 | 95% at 5% excluded | 33.1 dB | 0.137 | 0.114 | 0.004 | 98% |
| location (6 m) | 0.044 | near chance | 20.2 | 0.008 | 0.590 | 95% at 5% excluded | 36 dB | 0.120 | 0.104 | 0.003 | 97% |
| shape (room or L tunnel) | 0.056 | 95% at 5% excluded | 16.1 | 0.009 | 0.702 | 95% at 5% excluded | 34 dB | 0.137 | 0.116 | 0.004 | 97% |

SNR sweep (L1 and oracle, TV):

| pair | 20 dB | 30 dB | 40 dB | 50 dB | 60 dB |
|---|---|---|---|---|---|
| presence (room) | 0.033 / 0.156 | 0.033 / 0.466 | 0.033 / 0.951 | 0.033 / 1.000 | 0.033 / 1 |
| presence (L tunnel) | 0.062 / 0.284 | 0.062 / 0.750 | 0.062 / 1.000 | 0.062 / 1 | 0.062 / 1 |
| location (6 m) | 0.044 / 0.206 | 0.044 / 0.590 | 0.044 / 0.991 | 0.044 / 1.000 | 0.044 / 1 |
| shape (room or L tunnel) | 0.056 / 0.258 | 0.056 / 0.702 | 0.056 / 0.999 | 0.056 / 1 | 0.056 / 1 |

## strong: allowances

| pair | allowance | L1 | oracle |
|---|---|---|---|
| presence (room) | grid allowance (x1.11 amplitude, empirical) | 0.037 (near chance) | 0.510 (95% at 5% excluded) |
| presence (room) | tail envelope x2 | 0.047 (near chance) | 0.617 (95% at 5% excluded) |
| presence (room) | tail envelope x5 | 0.074 (95% at 5% excluded) | 0.830 (95% at 5% excluded) |
| presence (room) | attenuation Q = 50 (illustrative) | 0.011 (near chance) | 0.156 (95% at 5% excluded) |
| presence (room) | attenuation Q = 20 (illustrative) | 0.008 (near chance) | 0.119 (95% at 5% excluded) |
| presence (L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.069 (95% at 5% excluded) | 0.799 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x2 | 0.087 (95% at 5% excluded) | 0.895 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x5 | 0.138 (95% at 5% excluded) | 0.989 (unresolved) |
| presence (L tunnel) | attenuation Q = 50 (illustrative) | 0.019 (near chance) | 0.273 (95% at 5% excluded) |
| presence (L tunnel) | attenuation Q = 20 (illustrative) | 0.014 (near chance) | 0.204 (95% at 5% excluded) |
| location (6 m) | grid allowance (x1.11 amplitude, empirical) | 0.049 (near chance) | 0.640 (95% at 5% excluded) |
| location (6 m) | tail envelope x2 | 0.062 (95% at 5% excluded) | 0.752 (95% at 5% excluded) |
| location (6 m) | tail envelope x5 | 0.098 (95% at 5% excluded) | 0.931 (unresolved) |
| location (6 m) | attenuation Q = 50 (illustrative) | 0.014 (near chance) | 0.207 (95% at 5% excluded) |
| location (6 m) | attenuation Q = 20 (illustrative) | 0.011 (near chance) | 0.161 (95% at 5% excluded) |
| shape (room or L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.062 (95% at 5% excluded) | 0.752 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x2 | 0.079 (95% at 5% excluded) | 0.856 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x5 | 0.124 (95% at 5% excluded) | 0.979 (unresolved) |
| shape (room or L tunnel) | attenuation Q = 50 (illustrative) | 0.017 (near chance) | 0.253 (95% at 5% excluded) |
| shape (room or L tunnel) | attenuation Q = 20 (illustrative) | 0.013 (near chance) | 0.193 (95% at 5% excluded) |

## strong: checks

- presence (room): fixed scenes exact/formula 0.994 ± 0.006 (common wave ×10: 0.994 ± 0.006); one line, exact KL / certificate ×1: 0.157, ×10: 0.153, ×100: 0.115, ×1000: none; line floor 0.9952; σ(r) slope over 40–70 m -0.43, σ(67.5)/σ(42.5) 0.83
- presence (L tunnel): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.169, ×10: 0.164, ×100: 0.120, ×1000: none; line floor 0.9911; σ(r) slope over 40–70 m -0.34, σ(67.5)/σ(42.5) 0.86
- location (6 m): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.125, ×10: 0.121, ×100: 0.083, ×1000: none; line floor 0.9941; σ(r) slope over 40–70 m -0.36, σ(67.5)/σ(42.5) 0.85
- shape (room or L tunnel): fixed scenes exact/formula 0.996 ± 0.006 (common wave ×10: 0.996 ± 0.006); one line, exact KL / certificate ×1: 0.181, ×10: 0.175, ×100: 0.122, ×1000: none; line floor 0.9927; σ(r) slope over 40–70 m -0.58, σ(67.5)/σ(42.5) 0.78

## strong: the score test on synthesised images (presence, room; 24 grounds, a cyclic patch 45 × 28 m)

| amplification | AUC empirical | AUC predicted | shift / sd | predicted deflection | found at 5%, empirical | predicted | patch certificate |
|---|---|---|---|---|---|---|---|
| ×1 | 0.521 | 0.501 | 0.00 | 0.003 | 0.08 | 0.05 | 0.004 |
| ×800 | 0.960 | 0.930 | 2.40 | 2.083 | 0.88 | 0.67 | 1 |

## strong: over the FTA band (every 2 Hz)

| f (Hz) | presence (room) | presence (L tunnel) | location (6 m) | shape (room or L tunnel) |
|---|---|---|---|---|
| 8 | 0.006 / 0.086 | 0.021 / 0.284 | 0.004 / 0.062 | 0.017 / 0.237 |
| 10 | 0.006 / 0.082 | 0.017 / 0.239 | 0.004 / 0.055 | 0.014 / 0.191 |
| 12 | 0.006 / 0.084 | 0.015 / 0.215 | 0.004 / 0.053 | 0.011 / 0.163 |
| 14 | 0.006 / 0.090 | 0.014 / 0.206 | 0.004 / 0.056 | 0.010 / 0.144 |
| 16 | 0.007 / 0.099 | 0.014 / 0.207 | 0.004 / 0.062 | 0.009 / 0.132 |
| 18 | 0.008 / 0.111 | 0.015 / 0.212 | 0.005 / 0.072 | 0.009 / 0.125 |
| 20 | 0.009 / 0.125 | 0.015 / 0.219 | 0.006 / 0.086 | 0.008 / 0.122 |
| 22 | 0.010 / 0.142 | 0.016 / 0.226 | 0.007 / 0.103 | 0.008 / 0.122 |
| 24 | 0.011 / 0.161 | 0.016 / 0.234 | 0.009 / 0.125 | 0.009 / 0.126 |
| 26 | 0.013 / 0.185 | 0.017 / 0.242 | 0.010 / 0.151 | 0.009 / 0.135 |
| 28 | 0.014 / 0.211 | 0.017 / 0.248 | 0.012 / 0.182 | 0.010 / 0.152 |
| 30 | 0.016 / 0.238 | 0.018 / 0.254 | 0.015 / 0.218 | 0.012 / 0.179 |
| 32 | 0.018 / 0.266 | 0.018 / 0.261 | 0.018 / 0.257 | 0.015 / 0.214 |
| 34 | 0.020 / 0.292 | 0.019 / 0.269 | 0.021 / 0.298 | 0.018 / 0.258 |
| 36 | 0.022 / 0.317 | 0.019 / 0.280 | 0.024 / 0.338 | 0.021 / 0.308 |
| 38 | 0.024 / 0.341 | 0.021 / 0.298 | 0.027 / 0.379 | 0.025 / 0.363 |
| 40 | 0.026 / 0.364 | 0.023 / 0.327 | 0.030 / 0.418 | 0.030 / 0.420 |
| 42 | 0.027 / 0.385 | 0.026 / 0.366 | 0.033 / 0.456 | 0.035 / 0.478 |
| 44 | 0.029 / 0.405 | 0.029 / 0.412 | 0.036 / 0.491 | 0.039 / 0.534 |
| 46 | 0.030 / 0.422 | 0.033 / 0.462 | 0.038 / 0.521 | 0.044 / 0.587 |
| 48 | 0.031 / 0.436 | 0.038 / 0.515 | 0.040 / 0.547 | 0.049 / 0.633 |
| 50 | 0.032 / 0.447 | 0.042 / 0.565 | 0.042 / 0.568 | 0.053 / 0.673 |
| 52 | 0.033 / 0.454 | 0.046 / 0.611 | 0.044 / 0.584 | 0.056 / 0.705 |
| 54 | 0.033 / 0.458 | 0.050 / 0.650 | 0.045 / 0.594 | 0.059 / 0.728 |
| 56 | 0.033 / 0.459 | 0.054 / 0.683 | 0.045 / 0.598 | 0.061 / 0.744 |
| 58 | 0.033 / 0.458 | 0.057 / 0.707 | 0.045 / 0.598 | 0.062 / 0.751 |
| 60 | 0.033 / 0.455 | 0.059 / 0.725 | 0.045 / 0.594 | 0.062 / 0.751 |
| 62 | 0.032 / 0.452 | 0.060 / 0.737 | 0.044 / 0.588 | 0.061 / 0.746 |
| 64 | 0.032 / 0.451 | 0.061 / 0.745 | 0.044 / 0.588 | 0.060 / 0.734 |
| 66 | 0.032 / 0.454 | 0.062 / 0.749 | 0.045 / 0.592 | 0.058 / 0.721 |
| 68 | 0.033 / 0.462 | 0.062 / 0.750 | 0.045 / 0.592 | 0.057 / 0.708 |
| 70 | 0.034 / 0.469 | 0.062 / 0.751 | 0.044 / 0.587 | 0.055 / 0.695 |
| 72 | 0.034 / 0.477 | 0.062 / 0.751 | 0.043 / 0.577 | 0.054 / 0.685 |
| 74 | 0.035 / 0.483 | 0.062 / 0.754 | 0.042 / 0.562 | 0.053 / 0.678 |
| 76 | 0.035 / 0.488 | 0.063 / 0.757 | 0.040 / 0.542 | 0.053 / 0.675 |
| 78 | 0.036 / 0.492 | 0.063 / 0.761 | 0.038 / 0.519 | 0.053 / 0.676 |
| 80 | 0.036 / 0.496 | 0.064 / 0.767 | 0.036 / 0.494 | 0.054 / 0.680 |
| 82 | 0.036 / 0.501 | 0.065 / 0.774 | 0.034 / 0.470 | 0.055 / 0.690 |
| 84 | 0.037 / 0.507 | 0.066 / 0.781 | 0.032 / 0.445 | 0.056 / 0.702 |
| 86 | 0.037 / 0.513 | 0.067 / 0.787 | 0.030 / 0.422 | 0.057 / 0.713 |
| 88 | 0.038 / 0.519 | 0.068 / 0.793 | 0.028 / 0.401 | 0.058 / 0.722 |
| 90 | 0.038 / 0.525 | 0.069 / 0.798 | 0.027 / 0.384 | 0.059 / 0.729 |
| 92 | 0.039 / 0.531 | 0.069 / 0.803 | 0.026 / 0.370 | 0.060 / 0.734 |
| 94 | 0.039 / 0.536 | 0.070 / 0.807 | 0.025 / 0.359 | 0.060 / 0.738 |
| 96 | 0.040 / 0.542 | 0.071 / 0.812 | 0.024 / 0.349 | 0.061 / 0.742 |
| 98 | 0.040 / 0.548 | 0.071 / 0.816 | 0.024 / 0.341 | 0.061 / 0.746 |
| 100 | 0.041 / 0.553 | 0.072 / 0.820 | 0.023 / 0.336 | 0.062 / 0.750 |

Each cell of the band table: L1 / oracle at the nominal SNR.
