# P2-36 benchmark (revision 2): generated report

Frozen parameters hash `4bc4c58f8948`; frozen frequency 69 Hz (the presence (room) pair's oracle on the measured disc, largest). SNR 12.66 dB calibrated (the site's median sigma0 over ICEYE's best specified NESZ; 18 dB generous). States: near chance (TV < 0.05), 95% at 5% excluded (TV < 0.9), unresolved. Each entry bounds detection rate minus false-alarm rate; 'predicted' is a weak-signal calculation, 'achieved' a simulated detector.

## quiet

| pair | L1 whole image | state | L1 within 70 m | growth to 0.9 | oracle | state | oracle bound stops excluding 95/5 at | oracle within 70 m | predicted (score test, within 70 m) | tail share | L1 averaged over realisations (linear bound) | worst realisation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 2.2e-07 | near chance | 2.1e-07 | 3.6e+06 | 4.5e-07 | near chance | 142 dB | 4.3e-07 | 2.8e-09 | 9% | 2.2e-07 | 9.6e-07 |
| presence (L tunnel) | 2.5e-07 | near chance | 2.4e-07 | 3.1e+06 | 5.1e-07 | near chance | 141 dB | 4.9e-07 | 3.1e-09 | 9% | 2.5e-07 | 1.1e-06 |
| location (6 m) | 1.9e-07 | near chance | 1.9e-07 | 4e+06 | 3.8e-07 | near chance | 143 dB | 3.8e-07 | 3e-09 | 1% | 1.9e-07 | 8.6e-07 |
| shape (room or L tunnel) | 1.8e-07 | near chance | 1.8e-07 | 4.1e+06 | 3.7e-07 | near chance | 144 dB | 3.6e-07 | 2.9e-09 | 1% | 1.8e-07 | 8.1e-07 |
| presence (10 m room) | 7.7e-07 | near chance | 7.1e-07 | 1.1e+06 | 1.6e-06 | near chance | 131 dB | 1.4e-06 | 7.1e-09 | 15% | 7.7e-07 | 3.4e-06 |

SNR sweep (L1 / oracle):

| pair | 9.66 dB | 12.66 dB | 18 dB | 30 dB | 40 dB |
|---|---|---|---|---|---|
| presence (room) | 2.2e-07 / 3.2e-07 | 2.2e-07 / 4.5e-07 | 2.3e-07 / 8.3e-07 | 2.3e-07 / 3.3e-06 | 2.3e-07 / 1e-05 |
| presence (L tunnel) | 2.5e-07 / 3.6e-07 | 2.5e-07 / 5.1e-07 | 2.6e-07 / 9.4e-07 | 2.6e-07 / 3.8e-06 | 2.6e-07 / 1.2e-05 |
| location (6 m) | 1.8e-07 / 2.7e-07 | 1.9e-07 / 3.8e-07 | 1.9e-07 / 7e-07 | 2e-07 / 2.8e-06 | 2e-07 / 8.9e-06 |
| shape (room or L tunnel) | 1.8e-07 / 2.6e-07 | 1.8e-07 / 3.7e-07 | 1.9e-07 / 6.8e-07 | 1.9e-07 / 2.7e-06 | 1.9e-07 / 8.6e-06 |
| presence (10 m room) | 7.5e-07 / 1.1e-06 | 7.7e-07 / 1.6e-06 | 7.9e-07 / 2.9e-06 | 8e-07 / 1.1e-05 | 8e-07 / 3.6e-05 |

## strong

| pair | L1 whole image | state | L1 within 70 m | growth to 0.9 | oracle | state | oracle bound stops excluding 95/5 at | oracle within 70 m | predicted (score test, within 70 m) | tail share |
|---|---|---|---|---|---|---|---|---|---|---|
| presence (room) | 0.034 | near chance | 0.006 | 26.7 | 0.067 | 95% at 5% excluded | 38.4 dB | 0.012 | 0.002 | 97% |
| presence (L tunnel) | 0.062 | 95% at 5% excluded | 0.009 | 14.4 | 0.124 | 95% at 5% excluded | 33.1 dB | 0.019 | 0.004 | 98% |
| location (6 m) | 0.045 | near chance | 0.008 | 20.1 | 0.089 | 95% at 5% excluded | 36 dB | 0.016 | 0.003 | 97% |
| shape (room or L tunnel) | 0.056 | 95% at 5% excluded | 0.009 | 15.9 | 0.112 | 95% at 5% excluded | 34 dB | 0.019 | 0.004 | 97% |
| presence (10 m room) | 0.113 | 95% at 5% excluded | 0.017 | 7.93 | 0.223 | 95% at 5% excluded | 27.9 dB | 0.033 | 0.007 | 98% |

SNR sweep (L1 / oracle):

| pair | 9.66 dB | 12.66 dB | 18 dB | 30 dB | 40 dB |
|---|---|---|---|---|---|
| presence (room) | 0.033 / 0.048 | 0.034 / 0.067 | 0.034 / 0.124 | 0.035 / 0.466 | 0.035 / 0.951 |
| presence (L tunnel) | 0.061 / 0.088 | 0.062 / 0.124 | 0.064 / 0.228 | 0.064 / 0.750 | 0.064 / 1.000 |
| location (6 m) | 0.044 / 0.063 | 0.045 / 0.089 | 0.046 / 0.164 | 0.046 / 0.590 | 0.046 / 0.991 |
| shape (room or L tunnel) | 0.055 / 0.080 | 0.056 / 0.112 | 0.058 / 0.206 | 0.058 / 0.702 | 0.058 / 0.999 |
| presence (10 m room) | 0.110 / 0.159 | 0.113 / 0.223 | 0.116 / 0.400 | 0.117 / 0.963 | 0.117 / 1.000 |

## strong: allowances and sensitivity (the envelope is a declared assumption; these are robustness tests, not its verification)

| pair | case | L1 | oracle |
|---|---|---|---|
| presence (room) | grid allowance (x1.11 amplitude, empirical) | 0.037 (near chance) | 0.075 (95% at 5% excluded) |
| presence (room) | tail envelope x2 | 0.047 (near chance) | 0.094 (95% at 5% excluded) |
| presence (room) | tail envelope x5 | 0.075 (95% at 5% excluded) | 0.148 (95% at 5% excluded) |
| presence (room) | attenuation Q = 50 (illustrative) | 0.011 (near chance) | 0.021 (near chance) |
| presence (room) | attenuation Q = 20 (illustrative) | 0.008 (near chance) | 0.016 (near chance) |
| presence (L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.069 (95% at 5% excluded) | 0.138 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x2 | 0.088 (95% at 5% excluded) | 0.174 (95% at 5% excluded) |
| presence (L tunnel) | tail envelope x5 | 0.139 (95% at 5% excluded) | 0.271 (95% at 5% excluded) |
| presence (L tunnel) | attenuation Q = 50 (illustrative) | 0.019 (near chance) | 0.038 (near chance) |
| presence (L tunnel) | attenuation Q = 20 (illustrative) | 0.014 (near chance) | 0.028 (near chance) |
| location (6 m) | grid allowance (x1.11 amplitude, empirical) | 0.050 (near chance) | 0.099 (95% at 5% excluded) |
| location (6 m) | tail envelope x2 | 0.063 (95% at 5% excluded) | 0.125 (95% at 5% excluded) |
| location (6 m) | tail envelope x5 | 0.099 (95% at 5% excluded) | 0.195 (95% at 5% excluded) |
| location (6 m) | attenuation Q = 50 (illustrative) | 0.014 (near chance) | 0.028 (near chance) |
| location (6 m) | attenuation Q = 20 (illustrative) | 0.011 (near chance) | 0.022 (near chance) |
| shape (room or L tunnel) | grid allowance (x1.11 amplitude, empirical) | 0.063 (95% at 5% excluded) | 0.125 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x2 | 0.079 (95% at 5% excluded) | 0.157 (95% at 5% excluded) |
| shape (room or L tunnel) | tail envelope x5 | 0.125 (95% at 5% excluded) | 0.245 (95% at 5% excluded) |
| shape (room or L tunnel) | attenuation Q = 50 (illustrative) | 0.017 (near chance) | 0.035 (near chance) |
| shape (room or L tunnel) | attenuation Q = 20 (illustrative) | 0.013 (near chance) | 0.026 (near chance) |
| presence (10 m room) | grid allowance (x1.11 amplitude, empirical) | 0.126 (95% at 5% excluded) | 0.247 (95% at 5% excluded) |
| presence (10 m room) | tail envelope x2 | 0.160 (95% at 5% excluded) | 0.310 (95% at 5% excluded) |
| presence (10 m room) | tail envelope x5 | 0.253 (95% at 5% excluded) | 0.470 (95% at 5% excluded) |
| presence (10 m room) | attenuation Q = 50 (illustrative) | 0.034 (near chance) | 0.068 (95% at 5% excluded) |
| presence (10 m room) | attenuation Q = 20 (illustrative) | 0.025 (near chance) | 0.050 (95% at 5% excluded) |

## strong: checks

- presence (room): fixed scenes exact/formula 0.994 ± 0.006 (common wave ×10: 0.994 ± 0.006); one line, exact KL / certificate ×1: 0.157, ×10: 0.153, ×100: 0.115, ×1000: none; line floor 0.9952; σ(r) slope over 40–70 m -0.43, σ(67.5)/σ(42.5) 0.83
- presence (L tunnel): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.169, ×10: 0.164, ×100: 0.120, ×1000: none; line floor 0.9911; σ(r) slope over 40–70 m -0.34, σ(67.5)/σ(42.5) 0.86
- location (6 m): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.125, ×10: 0.121, ×100: 0.083, ×1000: none; line floor 0.9941; σ(r) slope over 40–70 m -0.36, σ(67.5)/σ(42.5) 0.85
- shape (room or L tunnel): fixed scenes exact/formula 0.996 ± 0.006 (common wave ×10: 0.996 ± 0.006); one line, exact KL / certificate ×1: 0.181, ×10: 0.175, ×100: 0.122, ×1000: none; line floor 0.9927; σ(r) slope over 40–70 m -0.58, σ(67.5)/σ(42.5) 0.78
- presence (10 m room): fixed scenes exact/formula 0.995 ± 0.006 (common wave ×10: 0.995 ± 0.006); one line, exact KL / certificate ×1: 0.135, ×10: 0.127, ×100: 0.055, ×1000: none; line floor 1.0486; σ(r) slope over 40–70 m -0.37, σ(67.5)/σ(42.5) 0.85

## strong: the score test achieved on synthesised images (presence, room; a cyclic patch 45 × 28 m; threshold from 24 calibration grounds, evaluation on 24 + 24 disjoint grounds)

| amplification | AUC achieved | AUC predicted | found, achieved | false alarms, achieved | found at 5%, predicted | patch certificate |
|---|---|---|---|---|---|---|
| ×1 | 0.530 | 0.501 | 0.08 | 0.21 | 0.05 | 0.004 |
| ×800 | 0.941 | 0.930 | 0.88 | 0.21 | 0.67 | 1 |

## strong: over the sampled FTA band (every 2 Hz; not a proved maximum between samples)

| f (Hz) | presence (room) | presence (L tunnel) | location (6 m) | shape (room or L tunnel) | presence (10 m room) |
|---|---|---|---|---|---|
| 8 | 0.006 / 0.012 | 0.021 / 0.039 | 0.004 / 0.008 | 0.017 / 0.033 | 0.029 / 0.057 |
| 10 | 0.006 / 0.011 | 0.017 / 0.033 | 0.004 / 0.007 | 0.014 / 0.026 | 0.029 / 0.057 |
| 12 | 0.006 / 0.011 | 0.015 / 0.030 | 0.004 / 0.007 | 0.012 / 0.022 | 0.031 / 0.060 |
| 14 | 0.006 / 0.012 | 0.014 / 0.028 | 0.004 / 0.008 | 0.010 / 0.020 | 0.033 / 0.066 |
| 16 | 0.007 / 0.014 | 0.014 / 0.028 | 0.004 / 0.008 | 0.009 / 0.018 | 0.037 / 0.073 |
| 18 | 0.008 / 0.015 | 0.015 / 0.029 | 0.005 / 0.010 | 0.009 / 0.017 | 0.041 / 0.082 |
| 20 | 0.009 / 0.017 | 0.015 / 0.030 | 0.006 / 0.012 | 0.008 / 0.017 | 0.046 / 0.092 |
| 22 | 0.010 / 0.019 | 0.016 / 0.031 | 0.007 / 0.014 | 0.008 / 0.017 | 0.052 / 0.103 |
| 24 | 0.011 / 0.022 | 0.016 / 0.032 | 0.009 / 0.017 | 0.009 / 0.017 | 0.059 / 0.117 |
| 26 | 0.013 / 0.025 | 0.017 / 0.033 | 0.010 / 0.021 | 0.009 / 0.018 | 0.067 / 0.133 |
| 28 | 0.015 / 0.029 | 0.017 / 0.034 | 0.013 / 0.025 | 0.010 / 0.021 | 0.076 / 0.149 |
| 30 | 0.017 / 0.033 | 0.018 / 0.035 | 0.015 / 0.030 | 0.012 / 0.024 | 0.084 / 0.166 |
| 32 | 0.018 / 0.037 | 0.018 / 0.036 | 0.018 / 0.036 | 0.015 / 0.029 | 0.093 / 0.183 |
| 34 | 0.020 / 0.041 | 0.019 / 0.037 | 0.021 / 0.041 | 0.018 / 0.036 | 0.100 / 0.197 |
| 36 | 0.022 / 0.044 | 0.020 / 0.039 | 0.024 / 0.047 | 0.022 / 0.043 | 0.106 / 0.209 |
| 38 | 0.024 / 0.048 | 0.021 / 0.041 | 0.027 / 0.053 | 0.026 / 0.051 | 0.111 / 0.218 |
| 40 | 0.026 / 0.051 | 0.023 / 0.046 | 0.030 / 0.060 | 0.030 / 0.060 | 0.113 / 0.223 |
| 42 | 0.027 / 0.055 | 0.026 / 0.051 | 0.033 / 0.066 | 0.035 / 0.069 | 0.114 / 0.225 |
| 44 | 0.029 / 0.058 | 0.029 / 0.059 | 0.036 / 0.071 | 0.040 / 0.079 | 0.113 / 0.223 |
| 46 | 0.030 / 0.060 | 0.033 / 0.067 | 0.038 / 0.077 | 0.044 / 0.088 | 0.111 / 0.218 |
| 48 | 0.031 / 0.062 | 0.038 / 0.076 | 0.041 / 0.081 | 0.049 / 0.098 | 0.107 / 0.211 |
| 50 | 0.032 / 0.064 | 0.042 / 0.084 | 0.043 / 0.085 | 0.053 / 0.106 | 0.103 / 0.203 |
| 52 | 0.033 / 0.065 | 0.047 / 0.093 | 0.044 / 0.088 | 0.057 / 0.113 | 0.100 / 0.198 |
| 54 | 0.033 / 0.066 | 0.051 / 0.101 | 0.045 / 0.090 | 0.060 / 0.119 | 0.099 / 0.196 |
| 56 | 0.033 / 0.066 | 0.054 / 0.108 | 0.046 / 0.091 | 0.062 / 0.122 | 0.101 / 0.199 |
| 58 | 0.033 / 0.066 | 0.057 / 0.114 | 0.046 / 0.091 | 0.063 / 0.124 | 0.104 / 0.205 |
| 60 | 0.033 / 0.066 | 0.059 / 0.118 | 0.045 / 0.090 | 0.063 / 0.124 | 0.107 / 0.211 |
| 62 | 0.033 / 0.065 | 0.061 / 0.121 | 0.045 / 0.089 | 0.062 / 0.123 | 0.109 / 0.216 |
| 64 | 0.032 / 0.065 | 0.062 / 0.123 | 0.044 / 0.089 | 0.060 / 0.120 | 0.111 / 0.220 |
| 66 | 0.033 / 0.065 | 0.062 / 0.124 | 0.045 / 0.089 | 0.059 / 0.117 | 0.113 / 0.222 |
| 68 | 0.033 / 0.067 | 0.062 / 0.124 | 0.045 / 0.089 | 0.057 / 0.114 | 0.113 / 0.223 |
| 70 | 0.034 / 0.068 | 0.062 / 0.124 | 0.044 / 0.089 | 0.056 / 0.111 | 0.113 / 0.223 |
| 72 | 0.035 / 0.069 | 0.063 / 0.125 | 0.044 / 0.087 | 0.054 / 0.108 | 0.113 / 0.222 |
| 74 | 0.035 / 0.070 | 0.063 / 0.125 | 0.042 / 0.084 | 0.054 / 0.107 | 0.112 / 0.220 |
| 76 | 0.036 / 0.071 | 0.063 / 0.126 | 0.040 / 0.080 | 0.053 / 0.106 | 0.110 / 0.217 |
| 78 | 0.036 / 0.072 | 0.064 / 0.127 | 0.038 / 0.076 | 0.054 / 0.107 | 0.109 / 0.216 |
| 80 | 0.036 / 0.072 | 0.065 / 0.129 | 0.036 / 0.072 | 0.054 / 0.108 | 0.109 / 0.214 |
| 82 | 0.037 / 0.073 | 0.066 / 0.131 | 0.034 / 0.068 | 0.055 / 0.110 | 0.107 / 0.212 |
| 84 | 0.037 / 0.074 | 0.067 / 0.133 | 0.032 / 0.064 | 0.056 / 0.112 | 0.106 / 0.209 |
| 86 | 0.038 / 0.075 | 0.068 / 0.134 | 0.030 / 0.060 | 0.058 / 0.115 | 0.104 / 0.206 |
| 88 | 0.038 / 0.076 | 0.068 / 0.136 | 0.029 / 0.057 | 0.059 / 0.117 | 0.102 / 0.202 |
| 90 | 0.039 / 0.077 | 0.069 / 0.138 | 0.027 / 0.054 | 0.060 / 0.119 | 0.100 / 0.198 |
| 92 | 0.039 / 0.078 | 0.070 / 0.139 | 0.026 / 0.052 | 0.060 / 0.120 | 0.098 / 0.193 |
| 94 | 0.040 / 0.079 | 0.071 / 0.140 | 0.025 / 0.050 | 0.061 / 0.121 | 0.096 / 0.189 |
| 96 | 0.040 / 0.080 | 0.071 / 0.142 | 0.025 / 0.049 | 0.061 / 0.122 | 0.094 / 0.185 |
| 98 | 0.041 / 0.081 | 0.072 / 0.143 | 0.024 / 0.048 | 0.062 / 0.123 | 0.092 / 0.182 |
| 100 | 0.041 / 0.082 | 0.073 / 0.145 | 0.024 / 0.047 | 0.062 / 0.124 | 0.090 / 0.179 |

Each cell of the band table: L1 / oracle at the calibrated SNR.
