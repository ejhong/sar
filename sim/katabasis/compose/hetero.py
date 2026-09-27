"""Von Karman random media: small-scale heterogeneity inside a stratum."""
from __future__ import annotations

import numpy as np


def von_karman(shape: tuple[int, int, int], spacing: float, corr_m: tuple[float, float, float],
               hurst: float, sigma: float, seed: int) -> np.ndarray:
    """Zero-mean field with standard deviation sigma and a von Karman spectrum.

    Power spectrum P(k) ∝ (1 + k²a²)^-(H + 3/2) in three dimensions, with
    anisotropic correlation lengths a = corr_m along (x, y, z).
    """
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(shape)
    ks = [np.fft.fftfreq(n, spacing) * 2 * np.pi for n in shape]
    KX, KY, KZ = np.meshgrid(*ks, indexing='ij', sparse=True)
    ka2 = (KX * corr_m[0]) ** 2 + (KY * corr_m[1]) ** 2 + (KZ * corr_m[2]) ** 2
    amp = (1 + ka2) ** (-(hurst + 1.5) / 2)
    field = np.real(np.fft.ifftn(np.fft.fftn(noise) * amp))
    field -= field.mean()
    sd = field.std()
    return field * (sigma / sd) if sd > 0 else field
