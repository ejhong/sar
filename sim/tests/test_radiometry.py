"""The SNR per cell is read from each acquisition's own radiometry, and traced to its sources."""
import unittest

from sarsim.radiometry import snr_per_cell


class RadiometryTests(unittest.TestCase):
    def test_giza_passes(self):
        for name in ('giza-20250827', 'giza-20220715'):
            s = snr_per_cell(name)
            self.assertFalse(s['nesz_in_product'])
            self.assertEqual(s['nesz_specified_db'], [-18.0, -15.0])
            self.assertIn('Table 2-11', s['nesz_source'])
            # the nominal is the median ground over the best specified floor; the measured ceiling only bounds from below
            self.assertAlmostEqual(s['nominal_db'], s['ground_sigma0_median_db'] + 18.0)
            self.assertLess(s['measured_floor_lower_db'], s['specified_worst_db'])
            self.assertLess(s['specified_worst_db'], s['nominal_db'])
            self.assertLess(s['nominal_db'], s['generous_db'])
            self.assertLess(s['nominal_db'], 20.0)      # far below the withdrawn 30 dB


if __name__ == '__main__':
    unittest.main()
