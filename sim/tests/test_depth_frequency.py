"""P2-06: the published focusing puts motion at a depth set by its frequency."""
import numpy as np
import pytest

from sarsim.acquisition import DwellGeometry
from sarsim.subap import SubapBank
from sarsim.tomo import focus_paper, kz_for_bank, nyquist_depth

G = DwellGeometry.from_record('giza-20250827')


def test_depth_is_frequency_times_a_constant_of_the_geometry():
    bank = SubapBank()
    t = G.nu_to_time(bank.bands(G, 113692)['nu_c'])
    kz = kz_for_bank(bank, G, 0.48, 113692)
    m_per_hz = 0.48 * G.R0 * np.sin(G.theta) / (2 * G.V_platform)
    assert 2 * np.pi / abs(np.polyfit(t, kz, 1)[0]) == pytest.approx(m_per_hz, rel=1e-4)
    z = np.linspace(0.02, 0.97 * nyquist_depth(kz), 1200)
    for f in (0.3, 1.0, 1.5):
        for amp in (1.0, 1e-9):                     # a chamber's imprint only scales the trajectory
            y = amp * np.cos(2 * np.pi * f * t + 0.3)[None, :].astype(complex)
            assert z[np.argmax(focus_paper(y, kz, z)[0])] == pytest.approx(m_per_hz * f, abs=0.05)
