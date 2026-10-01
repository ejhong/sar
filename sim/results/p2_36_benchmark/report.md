# P2-36 benchmark (revision 3): generated report

Frozen parameters hash `2684dc43b817` (worlds `dee606bbff0f`, unchanged since revision 2); frozen frequency 69 Hz (the presence (room) pair's oracle on the measured disc, largest). SNR conditional: the site's median sigma0 measured in the image over ICEYE's specified noise floor; headline 18.36 dB (Dwell Fine's best, documentation 6.0.8), 23.7 dB for ground at 0 dB, every other specified value and the measured lower bound in the sweep. States: near chance (TV < 0.05), 95% at 5% excluded (TV < 0.9), unresolved. Each entry bounds detection rate minus false-alarm rate; 'predicted' is a weak-signal calculation, 'achieved' a simulated detector.

## quiet

| pair | L1 whole image | state | L1 within 70 m | growth to 0.9 | oracle | state | oracle bound stops excluding 95/5 at | oracle within 70 m | predicted (score test, within 70 m) | tail share | L1 averaged over realisations (linear bound) | worst realisation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 2.3e-07 | near chance | 2.2e-07 | 2.3e+06 | 8.6e-07 | near chance | 142 dB | 8.2e-07 | 2.8e-09 | 9% | 2.3e-07 | 1e-06 |
| presence (L tunnel) | 2.6e-07 | near chance | 2.5e-07 | 1.9e+06 | 9.8e-07 | near chance | 141 dB | 9.4e-07 | 3.1e-09 | 9% | 2.6e-07 | 1.2e-06 |
| location (6 m) | 1.9e-07 | near chance | 1.9e-07 | 2.1e+06 | 7.3e-07 | near chance | 143 dB | 7.3e-07 | 3e-09 | 1% | 1.9e-07 | 8.9e-07 |
| shape (room or L tunnel) | 1.9e-07 | near chance | 1.9e-07 | 2.3e+06 | 7.1e-07 | near chance | 144 dB | 7e-07 | 2.9e-09 | 1% | 1.9e-07 | 8.4e-07 |
| presence (10 m room) | 7.9e-07 | near chance | 7.3e-07 | 7.1e+05 | 3e-06 | near chance | 131 dB | 2.8e-06 | 7.1e-09 | 15% | 7.9e-07 | 3.5e-06 |

SNR sweep (L1 / oracle):

| pair | 5.88 dB | 7.26 dB | 9.66 dB | 12.66 dB | 18.36 dB | 23.7 dB | 30 dB | 40 dB |
|---|---|---|---|---|---|---|---|---|
| presence (room) | 2e-07 / 2e-07 | 2.1e-07 / 2.4e-07 | 2.2e-07 / 3.2e-07 | 2.2e-07 / 4.5e-07 | 2.3e-07 / 8.6e-07 | 2.3e-07 / 1.6e-06 | 2.3e-07 / 3.3e-06 | 2.3e-07 / 1e-05 |
| presence (L tunnel) | 2.3e-07 / 2.3e-07 | 2.4e-07 / 2.7e-07 | 2.5e-07 / 3.6e-07 | 2.5e-07 / 5.1e-07 | 2.6e-07 / 9.8e-07 | 2.6e-07 / 1.8e-06 | 2.6e-07 / 3.8e-06 | 2.6e-07 / 1.2e-05 |
| location (6 m) | 1.7e-07 / 1.7e-07 | 1.8e-07 / 2e-07 | 1.8e-07 / 2.7e-07 | 1.9e-07 / 3.8e-07 | 1.9e-07 / 7.3e-07 | 2e-07 / 1.4e-06 | 2e-07 / 2.8e-06 | 2e-07 / 8.9e-06 |
| shape (room or L tunnel) | 1.7e-07 / 1.7e-07 | 1.7e-07 / 2e-07 | 1.8e-07 / 2.6e-07 | 1.8e-07 / 3.7e-07 | 1.9e-07 / 7.1e-07 | 1.9e-07 / 1.3e-06 | 1.9e-07 / 2.7e-06 | 1.9e-07 / 8.6e-06 |
| presence (10 m room) | 7e-07 / 7.1e-07 | 7.3e-07 / 8.3e-07 | 7.5e-07 / 1.1e-06 | 7.7e-07 / 1.6e-06 | 7.9e-07 / 3e-06 | 7.9e-07 / 5.5e-06 | 8e-07 / 1.1e-05 | 8e-07 / 3.6e-05 |

## strong

| pair | L1 whole image | state | L1 within 70 m | growth to 0.9 | oracle | state | oracle bound stops excluding 95/5 at | oracle within 70 m | predicted (score test, within 70 m) | tail share |
|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 0.034 | near chance | 0.006 | 26.1 | 0.129 | 95% at 5% excluded | 38.4 dB | 0.023 | 0.002 | 97% |
| presence (L tunnel) | 0.064 | 95% at 5% excluded | 0.010 | 14.1 | 0.237 | 95% at 5% excluded | 33.1 dB | 0.036 | 0.004 | 98% |
| location (6 m) | 0.046 | near chance | 0.008 | 19.6 | 0.171 | 95% at 5% excluded | 36 dB | 0.031 | 0.003 | 97% |
| shape (room or L tunnel) | 0.058 | 95% at 5% excluded | 0.010 | 15.6 | 0.215 | 95% at 5% excluded | 34 dB | 0.036 | 0.004 | 97% |
| presence (10 m room) | 0.116 | 95% at 5% excluded | 0.017 | 7.77 | 0.415 | 95% at 5% excluded | 27.9 dB | 0.064 | 0.007 | 98% |

SNR sweep (L1 / oracle):

| pair | 5.88 dB | 7.26 dB | 9.66 dB | 12.66 dB | 18.36 dB | 23.7 dB | 30 dB | 40 dB |
|---|---|---|---|---|---|---|---|---|
| presence (room) | 0.031 / 0.031 | 0.032 / 0.036 | 0.033 / 0.048 | 0.034 / 0.067 | 0.034 / 0.129 | 0.035 / 0.237 | 0.035 / 0.466 | 0.035 / 0.951 |
| presence (L tunnel) | 0.057 / 0.057 | 0.059 / 0.067 | 0.061 / 0.088 | 0.062 / 0.124 | 0.064 / 0.237 | 0.064 / 0.423 | 0.064 / 0.750 | 0.064 / 1.000 |
| location (6 m) | 0.041 / 0.041 | 0.042 / 0.048 | 0.044 / 0.063 | 0.045 / 0.089 | 0.046 / 0.171 | 0.046 / 0.310 | 0.046 / 0.590 | 0.046 / 0.991 |
| shape (room or L tunnel) | 0.051 / 0.052 | 0.053 / 0.060 | 0.055 / 0.080 | 0.056 / 0.112 | 0.058 / 0.215 | 0.058 / 0.385 | 0.058 / 0.702 | 0.058 / 0.999 |
| presence (10 m room) | 0.103 / 0.103 | 0.106 / 0.121 | 0.110 / 0.159 | 0.113 / 0.223 | 0.116 / 0.415 | 0.116 / 0.687 | 0.117 / 0.963 | 0.117 / 1.000 |

## strong: allowances and sensitivity (the envelope is a declared assumption; these are robustness tests, not its verification)

| pair | case | L1 | oracle |
|---|---|---|---|
| presence (room) | grid allowance (x1.11 amplitude, empirical) | 0.038 (near chance) | 0.143 (95% at 5% excluded) |
| presence (room) | tail envelope x2 | 0.048 (near chance) | 0.181 (95% at 5% excluded) |
| presence (room) | tail envelope x5 | 0.076 (95% at 5% excluded) | 0.281 (95% at 5% excluded) |
| presence (room) | attenuation Q = 50 (illustrative) | 0.011 (near chance) | 0.041 (near chance) |
| presence (room) | attenuation Q = 20 (illustrative) | 0.008 (near chance) | 0.031 (near chance) |
| presence (L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.071 (95% at 5% excluded) | 0.262 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x2 | 0.090 (95% at 5% excluded) | 0.328 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x5 | 0.142 (95% at 5% excluded) | 0.496 (95% at 5% excluded) |
| presence (L tunnel) | attenuation Q = 50 (illustrative) | 0.019 (near chance) | 0.073 (95% at 5% excluded) |
| presence (L tunnel) | attenuation Q = 20 (illustrative) | 0.014 (near chance) | 0.054 (95% at 5% excluded) |
| location (6 m) | grid allowance (x1.11 amplitude, empirical) | 0.051 (95% at 5% excluded) | 0.189 (95% at 5% excluded) |
| location (6 m) | tail envelope x2 | 0.064 (95% at 5% excluded) | 0.238 (95% at 5% excluded) |
| location (6 m) | tail envelope x5 | 0.101 (95% at 5% excluded) | 0.366 (95% at 5% excluded) |
| location (6 m) | attenuation Q = 50 (illustrative) | 0.015 (near chance) | 0.055 (95% at 5% excluded) |
| location (6 m) | attenuation Q = 20 (illustrative) | 0.011 (near chance) | 0.042 (near chance) |
| shape (room or L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.064 (95% at 5% excluded) | 0.238 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x2 | 0.081 (95% at 5% excluded) | 0.298 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x5 | 0.128 (95% at 5% excluded) | 0.453 (95% at 5% excluded) |
| shape (room or L tunnel) | attenuation Q = 50 (illustrative) | 0.018 (near chance) | 0.067 (95% at 5% excluded) |
| shape (room or L tunnel) | attenuation Q = 20 (illustrative) | 0.014 (near chance) | 0.051 (95% at 5% excluded) |
| presence (10 m room) | grid allowance (x1.11 amplitude, empirical) | 0.128 (95% at 5% excluded) | 0.456 (95% at 5% excluded) |
| presence (10 m room) | tail envelope x2 | 0.163 (95% at 5% excluded) | 0.558 (95% at 5% excluded) |
| presence (10 m room) | tail envelope x5 | 0.258 (95% at 5% excluded) | 0.774 (95% at 5% excluded) |
| presence (10 m room) | attenuation Q = 50 (illustrative) | 0.035 (near chance) | 0.131 (95% at 5% excluded) |
| presence (10 m room) | attenuation Q = 20 (illustrative) | 0.026 (near chance) | 0.096 (95% at 5% excluded) |

## strong: checks

- presence (room): fixed scenes exact/formula 0.994 ± 0.006 (common wave ×10: 0.994 ± 0.006); one line, exact KL / certificate ×1: 0.155, ×10: 0.151, ×100: 0.114, ×1000: none; line floor 1.0088; σ(r) slope over 40–70 m -0.43, σ(67.5)/σ(42.5) 0.83
- presence (L tunnel): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.166, ×10: 0.162, ×100: 0.119, ×1000: none; line floor 1.0047; σ(r) slope over 40–70 m -0.34, σ(67.5)/σ(42.5) 0.86
- location (6 m): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.123, ×10: 0.119, ×100: 0.082, ×1000: none; line floor 1.0077; σ(r) slope over 40–70 m -0.36, σ(67.5)/σ(42.5) 0.85
- shape (room or L tunnel): fixed scenes exact/formula 0.996 ± 0.006 (common wave ×10: 0.996 ± 0.006); one line, exact KL / certificate ×1: 0.178, ×10: 0.173, ×100: 0.121, ×1000: none; line floor 1.0063; σ(r) slope over 40–70 m -0.58, σ(67.5)/σ(42.5) 0.78
- presence (10 m room): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.140, ×10: 0.132, ×100: 0.056, ×1000: none; line floor 1.0090; σ(r) slope over 40–70 m -0.37, σ(67.5)/σ(42.5) 0.85

## strong: the score test achieved on synthesised images (presence, room; a cyclic patch 45 × 28 m; threshold from 400 calibration grounds, evaluation on 100 + 100 grounds disjoint from them, seed streams {'noise_law': [35900, 35901], 'sign': [35950, 35951], 'calibration': [36000, 36400], 'null_evaluation': [37000, 37100], 'cavity_evaluation': [38000, 38100]})

| amplification | AUC achieved ± se (p) | AUC predicted | found, achieved [95%] | false alarms, achieved [95%] | found > false alarms, p | found at 5%, predicted | patch certificate |
|---|---|---|---|---|---|---|---|
| ×1 | 0.469 ± 0.041 (0.78) | 0.501 | 7/100 [0.03, 0.14] | 15/100 [0.09, 0.24] | 0.98 | 0.05 | 0.004 |
| ×800 | 0.913 ± 0.021 (3.2e-24) | 0.930 | 73/100 [0.63, 0.81] | 15/100 [0.09, 0.24] | 2.4e-17 | 0.67 | 1 |

## strong: over the sampled FTA band (every 2 Hz; not a proved maximum between samples)

| f (Hz) | presence (room) | presence (L tunnel) | location (6 m) | shape (room or L tunnel) | presence (10 m room) |
|---|---|---|---|---|---|
| 8 | 0.006 / 0.023 | 0.021 / 0.076 | 0.004 / 0.016 | 0.018 / 0.063 | 0.030 / 0.109 |
| 10 | 0.006 / 0.021 | 0.017 / 0.063 | 0.004 / 0.014 | 0.014 / 0.051 | 0.030 / 0.109 |
| 12 | 0.006 / 0.022 | 0.016 / 0.057 | 0.004 / 0.014 | 0.012 / 0.043 | 0.031 / 0.116 |
| 14 | 0.006 / 0.024 | 0.015 / 0.055 | 0.004 / 0.015 | 0.010 / 0.038 | 0.034 / 0.127 |
| 16 | 0.007 / 0.026 | 0.015 / 0.055 | 0.004 / 0.016 | 0.009 / 0.035 | 0.038 / 0.141 |
| 18 | 0.008 / 0.029 | 0.015 / 0.056 | 0.005 / 0.019 | 0.009 / 0.033 | 0.042 / 0.157 |
| 20 | 0.009 / 0.033 | 0.016 / 0.058 | 0.006 / 0.022 | 0.009 / 0.032 | 0.047 / 0.176 |
| 22 | 0.010 / 0.037 | 0.016 / 0.060 | 0.007 / 0.027 | 0.009 / 0.032 | 0.053 / 0.197 |
| 24 | 0.011 / 0.042 | 0.017 / 0.062 | 0.009 / 0.033 | 0.009 / 0.033 | 0.060 / 0.223 |
| 26 | 0.013 / 0.049 | 0.017 / 0.064 | 0.011 / 0.040 | 0.009 / 0.036 | 0.068 / 0.253 |
| 28 | 0.015 / 0.056 | 0.018 / 0.066 | 0.013 / 0.048 | 0.011 / 0.040 | 0.077 / 0.284 |
| 30 | 0.017 / 0.063 | 0.018 / 0.068 | 0.015 / 0.058 | 0.013 / 0.047 | 0.086 / 0.315 |
| 32 | 0.019 / 0.071 | 0.019 / 0.069 | 0.018 / 0.068 | 0.015 / 0.057 | 0.094 / 0.344 |
| 34 | 0.021 / 0.078 | 0.019 / 0.072 | 0.021 / 0.080 | 0.018 / 0.069 | 0.102 / 0.369 |
| 36 | 0.023 / 0.085 | 0.020 / 0.075 | 0.024 / 0.091 | 0.022 / 0.083 | 0.108 / 0.390 |
| 38 | 0.024 / 0.092 | 0.021 / 0.080 | 0.027 / 0.103 | 0.026 / 0.098 | 0.113 / 0.405 |
| 40 | 0.026 / 0.099 | 0.023 / 0.088 | 0.031 / 0.115 | 0.031 / 0.115 | 0.116 / 0.415 |
| 42 | 0.028 / 0.105 | 0.026 / 0.099 | 0.034 / 0.126 | 0.036 / 0.133 | 0.117 / 0.418 |
| 44 | 0.029 / 0.111 | 0.030 / 0.113 | 0.037 / 0.137 | 0.040 / 0.152 | 0.116 / 0.415 |
| 46 | 0.031 / 0.116 | 0.034 / 0.128 | 0.039 / 0.147 | 0.045 / 0.170 | 0.113 / 0.406 |
| 48 | 0.032 / 0.120 | 0.039 / 0.145 | 0.042 / 0.156 | 0.050 / 0.187 | 0.109 / 0.393 |
| 50 | 0.033 / 0.123 | 0.043 / 0.162 | 0.044 / 0.163 | 0.054 / 0.202 | 0.105 / 0.380 |
| 52 | 0.033 / 0.126 | 0.048 / 0.178 | 0.045 / 0.169 | 0.058 / 0.216 | 0.102 / 0.370 |
| 54 | 0.034 / 0.127 | 0.052 / 0.193 | 0.046 / 0.172 | 0.061 / 0.226 | 0.102 / 0.368 |
| 56 | 0.034 / 0.127 | 0.055 / 0.206 | 0.046 / 0.174 | 0.063 / 0.234 | 0.103 / 0.373 |
| 58 | 0.034 / 0.127 | 0.058 / 0.217 | 0.046 / 0.174 | 0.064 / 0.237 | 0.106 / 0.383 |
| 60 | 0.034 / 0.126 | 0.061 / 0.225 | 0.046 / 0.172 | 0.064 / 0.237 | 0.109 / 0.393 |
| 62 | 0.033 / 0.125 | 0.062 / 0.231 | 0.046 / 0.170 | 0.063 / 0.235 | 0.112 / 0.402 |
| 64 | 0.033 / 0.125 | 0.063 / 0.234 | 0.045 / 0.170 | 0.062 / 0.229 | 0.114 / 0.409 |
| 66 | 0.033 / 0.126 | 0.064 / 0.236 | 0.046 / 0.171 | 0.060 / 0.223 | 0.115 / 0.413 |
| 68 | 0.034 / 0.128 | 0.064 / 0.237 | 0.046 / 0.171 | 0.058 / 0.217 | 0.116 / 0.415 |
| 70 | 0.035 / 0.130 | 0.064 / 0.237 | 0.045 / 0.170 | 0.057 / 0.212 | 0.116 / 0.415 |
| 72 | 0.035 / 0.133 | 0.064 / 0.237 | 0.044 / 0.166 | 0.056 / 0.207 | 0.115 / 0.413 |
| 74 | 0.036 / 0.135 | 0.064 / 0.239 | 0.043 / 0.161 | 0.055 / 0.205 | 0.114 / 0.410 |
| 76 | 0.036 / 0.136 | 0.065 / 0.240 | 0.041 / 0.154 | 0.055 / 0.203 | 0.113 / 0.405 |
| 78 | 0.037 / 0.138 | 0.065 / 0.242 | 0.039 / 0.146 | 0.055 / 0.204 | 0.112 / 0.402 |
| 80 | 0.037 / 0.139 | 0.066 / 0.245 | 0.037 / 0.138 | 0.055 / 0.206 | 0.111 / 0.399 |
| 82 | 0.037 / 0.141 | 0.067 / 0.249 | 0.035 / 0.130 | 0.056 / 0.210 | 0.110 / 0.395 |
| 84 | 0.038 / 0.142 | 0.068 / 0.252 | 0.033 / 0.123 | 0.058 / 0.215 | 0.108 / 0.390 |
| 86 | 0.038 / 0.144 | 0.069 / 0.256 | 0.031 / 0.116 | 0.059 / 0.220 | 0.106 / 0.384 |
| 88 | 0.039 / 0.146 | 0.070 / 0.259 | 0.029 / 0.110 | 0.060 / 0.224 | 0.104 / 0.378 |
| 90 | 0.040 / 0.148 | 0.071 / 0.262 | 0.028 / 0.104 | 0.061 / 0.227 | 0.102 / 0.370 |
| 92 | 0.040 / 0.150 | 0.071 / 0.264 | 0.027 / 0.100 | 0.062 / 0.229 | 0.100 / 0.363 |
| 94 | 0.041 / 0.152 | 0.072 / 0.267 | 0.026 / 0.097 | 0.062 / 0.231 | 0.098 / 0.356 |
| 96 | 0.041 / 0.154 | 0.073 / 0.270 | 0.025 / 0.094 | 0.063 / 0.233 | 0.096 / 0.349 |
| 98 | 0.042 / 0.156 | 0.074 / 0.272 | 0.024 / 0.092 | 0.063 / 0.235 | 0.094 / 0.343 |
| 100 | 0.042 / 0.158 | 0.074 / 0.275 | 0.024 / 0.090 | 0.064 / 0.237 | 0.092 / 0.337 |

Each cell of the band table: L1 / oracle at the headline SNR.
