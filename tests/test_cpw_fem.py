import importlib.util
import unittest

HAVE_FEMWELL = importlib.util.find_spec("femwell") is not None

if HAVE_FEMWELL:
    from tfln_mzm.cpw_fem import CPW, build_mesh, half_geometry, resolutions, solve, tutorial_domain

TUNCER = dict(w_sig=7.0, gap=10.0, w_gnd=100.0, t_metal=0.8, sigma=6e7, eps_sub=13.0)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestCPWGeometry(unittest.TestCase):
    def test_domain_and_ground_clipping(self):
        c = CPW(**TUNCER)
        self.assertEqual(tutorial_domain(c, 1), (40.5, 27.0, 27.0))
        shapes, off = half_geometry(c, 1)
        self.assertEqual(shapes["metal_gnd"].bounds, (13.5, 0.0, 40.5, 0.8))   # clipped to the domain
        self.assertEqual(half_geometry(c, 4)[0]["metal_gnd"].bounds, (13.5, 0.0, 113.5, 0.8))  # full 100 um
        self.assertAlmostEqual(off, 0.04)

    def test_contour_is_open_on_the_symmetry_plane(self):
        shapes, off = half_geometry(CPW(**TUNCER), 1)
        line = shapes["sig_contour"]
        self.assertEqual(line.geom_type, "LineString")
        ends = [line.coords[0], line.coords[-1]]
        self.assertTrue(all(abs(x) < 1e-12 for x, _ in ends))
        self.assertAlmostEqual(max(x for x, _ in line.coords), 3.5 + off)

    def test_mesh_scale_multiplies_sizes_only(self):
        a, b = resolutions(1.0), resolutions(0.5)
        for k in a:
            self.assertEqual(a[k]["distance"], b[k]["distance"])
            self.assertAlmostEqual(b[k]["resolution"], a[k]["resolution"] / 2)
            if "SizeMax" in a[k]:
                self.assertAlmostEqual(b[k]["SizeMax"], a[k]["SizeMax"] / 2)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestCPWLowFrequencyResistance(unittest.TestCase):
    def test_series_resistance_approaches_dc_value(self):
        """At 50 MHz the thin-sheet skin length (~100 um) exceeds every conductor width, so R ~ R_dc.

        With the tutorial's contour current this test failed by 24% on this mesh (A17 boundary addendum).
        """
        c = CPW(**TUNCER)
        r = solve(c, build_mesh(c, 4.0, 1), 5e7)
        g = 40.5 - 13.5   # ground width inside the d=1 half-domain
        r_dc = 1 / (c.sigma * c.w_sig * 1e-6 * c.t_metal * 1e-6) + 1 / (c.sigma * 2 * g * 1e-6 * c.t_metal * 1e-6)
        self.assertLess(abs(r["R"] / r_dc - 1), 0.02)
        self.assertLess(abs(r["gamma_rlgc"] / r["gamma_fem"] - 1), 0.05)
        self.assertGreater(r["Z0_pi"].real, 0)

    def test_second_order_series_resistance_approaches_dc_value(self):
        """Same check with second-order elements (A17b G1)."""
        c = CPW(**TUNCER)
        r = solve(c, build_mesh(c, 4.0, 1), 5e7, order=2)
        g = 40.5 - 13.5
        r_dc = 1 / (c.sigma * c.w_sig * 1e-6 * c.t_metal * 1e-6) + 1 / (c.sigma * 2 * g * 1e-6 * c.t_metal * 1e-6)
        self.assertLess(abs(r["R"] / r_dc - 1), 0.02)
        self.assertLess(abs(r["gamma_rlgc"] / r["gamma_fem"] - 1), 0.05)
        self.assertGreater(r["n_dofs"], 5 * r["n_elements"])   # second-order space, not the default


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestCPWModeTracking(unittest.TestCase):
    def test_complex_n_guess_returns_the_same_mode(self):
        """femwell accepts a complex n_guess and, where the eigenproblem is well conditioned, keeps the mode.

        At 50 MHz this check failed with second-order elements (neff moved 6e-3 with the shift): the E-field
        formulation breaks down at low frequency, so the answer drifts with round-off (A17b BOUNDARY section 8).
        """
        c = CPW(**TUNCER)
        mesh = build_mesh(c, 4.0, 1)
        default = solve(c, mesh, 2e9, order=2)
        self.assertGreater(abs(default["neff"].imag), 1e-3 * abs(default["neff"]))   # lossy, so genuinely complex
        tracked = solve(c, mesh, 2e9, order=2, n_guess=default["neff"])
        self.assertLess(abs(tracked["neff"] / default["neff"] - 1), 1e-5)
        self.assertLess(abs(tracked["Z0_pi"] / default["Z0_pi"] - 1), 1e-4)


if __name__ == "__main__":
    unittest.main()
