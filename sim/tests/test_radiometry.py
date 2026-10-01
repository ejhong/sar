"""The SNR per cell: the backscatter measured in each acquisition over ICEYE's specified noise, version by version."""
import unittest

from sarsim.radiometry import SPECIFICATION, snr_per_cell


class RadiometryTests(unittest.TestCase):
    def test_giza_passes(self):
        for name, mode, best in (('giza-20250827', 'SpotlightDwellFine', -23.7), ('giza-20220715', 'SpotlightDwell', -26.7)):
            s = snr_per_cell(name)
            self.assertFalse(s['nesz_in_product'])
            self.assertEqual(s['mode'], mode)
            self.assertEqual(s['nesz_headline_db'], best)
            self.assertIn('Table 2-11', s['nesz_headline_source'])
            self.assertIn('conditional', s['status'])
            # the headline is the median ground over the newest documentation's best for the mode: the most favourable
            # specified, so every other scenario lies below it; the measured ceiling bounds the ratio only from below
            self.assertAlmostEqual(s['headline_db'], s['ground_sigma0_median_db'] - best)
            self.assertTrue(all(sc['snr_db'] <= s['headline_db'] + 1e-12 for sc in s['scenarios']))
            self.assertLess(s['measured_floor_lower_db'], min(sc['snr_db'] for sc in s['scenarios']))
            self.assertLess(s['headline_db'], s['bright_db'])
            self.assertLessEqual(s['headline_db'], s['any_dwell_mode_best_db'])

    def test_versions_as_published(self):
        self.assertEqual(SPECIFICATION['6.0.0']['nesz_db']['SpotlightDwellFine'], (-18.0, -15.0))
        self.assertEqual(SPECIFICATION['6.0.8']['nesz_db']['SpotlightDwell'], (-26.7, -15.6))
        self.assertEqual(SPECIFICATION['6.0.8']['nesz_db']['SpotlightDwellFine'], (-23.7, -12.6))


if __name__ == '__main__':
    unittest.main()
