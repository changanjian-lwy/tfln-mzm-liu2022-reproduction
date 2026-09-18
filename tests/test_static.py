"""Physical invariants and an independent hand-calculated lecture benchmark."""
import unittest
from dataclasses import replace
import numpy as np
from tfln_mzm.static import StaticParameters, delta_n, phase_difference, half_wave_voltage, transmission, dc_gap_voltage


class TestStaticPhysics(unittest.TestCase):
    def setUp(self):
        self.p = StaticParameters()
        self.vpi = half_wave_voltage(self.p)

    def test_lecture_benchmark(self):
        self.assertAlmostEqual(self.vpi, 2.426120711244678, places=10)
        self.assertAlmostEqual(float(delta_n(self.vpi, self.p)), -3.875e-5, places=12)
        self.assertAlmostEqual(float(phase_difference(self.vpi, self.p)), -np.pi, places=12)

    def test_switching_and_period(self):
        np.testing.assert_allclose(transmission([0, self.vpi, 2*self.vpi], self.p), [1, 0, 1], atol=1e-14)
        v = np.linspace(-10, 10, 501)
        np.testing.assert_allclose(transmission(v, self.p), transmission(v+2*self.vpi, self.p), atol=1e-14)

    def test_field_interference_and_power_conservation(self):
        phi = phase_difference(np.linspace(-10, 10, 501), self.p)
        field_bright = (1 + np.exp(1j*phi))/2
        field_dark = (1 - np.exp(1j*phi))/2
        np.testing.assert_allclose(abs(field_bright)**2, transmission(np.linspace(-10, 10, 501), self.p), atol=1e-14)
        np.testing.assert_allclose(abs(field_bright)**2+abs(field_dark)**2, 1, atol=1e-14)

    def test_quadrature_slope(self):
        h = self.vpi*1e-5
        slope = (transmission(h, self.p, bias_rad=np.pi/2)-transmission(-h, self.p, bias_rad=np.pi/2))/(2*h)
        self.assertAlmostEqual(float(slope), np.pi/(2*self.vpi), places=8)

    def test_geometry_and_drive_scaling(self):
        self.assertAlmostEqual(half_wave_voltage(replace(self.p, length_m=0.04)), self.vpi/2)
        self.assertAlmostEqual(half_wave_voltage(replace(self.p, gap_m=20e-6)), 2*self.vpi)
        self.assertAlmostEqual(half_wave_voltage(replace(self.p, overlap=0.5)), 2*self.vpi)
        self.assertAlmostEqual(half_wave_voltage(self.p, "push_pull"), self.vpi/2)

    def test_dc_voltage_divider(self):
        np.testing.assert_allclose(dc_gap_voltage(1, [0, 20, 40, 50, 80, np.inf]), [0, 2/7, 4/9, 0.5, 8/13, 1])

    def test_invalid_parameters(self):
        for kw in ({"length_m":0}, {"overlap":1.1}, {"n_e":float('nan')}):
            with self.assertRaises(ValueError):
                StaticParameters(**kw)
        with self.assertRaises(ValueError):
            half_wave_voltage(self.p, "unspecified")


if __name__ == "__main__":
    unittest.main()
