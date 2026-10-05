import importlib.util
import unittest
from dataclasses import replace

import numpy as np

HAVE_FEMWELL = importlib.util.find_spec("femwell") is not None

if HAVE_FEMWELL:
    from tfln_mzm.optical_fem import (Section, build_slab_mesh, geometry, group_index, index, permittivity,
                                      slab_exact, solve_slab, LN_FILES, SILICA, GOLD)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestMaterials(unittest.TestCase):
    def test_sellmeier_values(self):
        """Values from a separate implementation of the published formulas (O01 boundary, cheaper check)."""
        self.assertAlmostEqual(index(LN_FILES["CLN"][0], 1.55).real, 2.13756, places=5)
        self.assertAlmostEqual(index(LN_FILES["CLN"][1], 1.55).real, 2.21111, places=5)
        self.assertAlmostEqual(index(SILICA, 1.55).real, 1.44402, places=5)
        self.assertAlmostEqual(index("SiO2-quartz_Ghosh-o.yml", 1.55).real, 1.52770, places=5)
        self.assertLess(index(LN_FILES["MgO"][0], 1.55).real, index(LN_FILES["MgO"][1], 1.55).real)  # negative uniaxial

    def test_gold_table_and_sign(self):
        n = index(GOLD, 1.55)   # Olmon-ev table row 1.550 um: 0.3226, 10.62
        self.assertAlmostEqual(n.real, 0.3226, places=6)
        self.assertAlmostEqual(n.imag, -10.62, places=6)
        self.assertLess((n * n).imag, 0)   # femwell lossy convention, as in cpw_fem

    def test_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            index(LN_FILES["CLN"][0], 6.0)

    def test_ln_tensor_order(self):
        eps = permittivity(Section(), 1.55)
        exx, eyy, ezz = eps["ridge"]
        self.assertAlmostEqual(exx, 2.13756**2, places=3)   # x = crystal Z: extraordinary
        self.assertAlmostEqual(eyy, ezz)
        self.assertGreater(eyy, exx)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestGeometry(unittest.TestCase):
    def test_widths(self):
        s = Section(theta_deg=60)
        top, bot = s.widths()
        self.assertAlmostEqual(top, 1.0)
        self.assertAlmostEqual(bot, 1.0 + 2 * 0.3 / np.tan(np.radians(60)))
        top, bot = replace(s, width_ref="bottom").widths()
        self.assertAlmostEqual(bot, 1.0)
        self.assertEqual(Section(theta_deg=90).widths(), (1.0, 1.0))

    def test_cladding_is_conformal(self):
        g = geometry(Section(theta_deg=90, rail_t=None))
        # vertical ridge: 7 um of slab-top band, a 1.2 um top cover (mitred corners included), two 0.2 um sides between them
        self.assertAlmostEqual(g["clad"].area, 0.1 * (8 - 1) + 0.1 * 1.2 + 2 * 0.1 * 0.2, places=9)
        self.assertAlmostEqual(g["clad"].bounds[1], 0.3)
        self.assertAlmostEqual(g["clad"].bounds[3], 0.7)

    def test_rails(self):
        g = geometry(Section(rail_t=0.5))
        self.assertEqual(g["rail_r"].bounds, (1.5, 0.4, 3.5, 0.9))
        self.assertEqual(g["rail_l"].bounds, (-3.5, 0.4, -1.5, 0.9))
        self.assertNotIn("rail_l", geometry(Section(rail_t=None)))

    def test_half_domain(self):
        full, half = geometry(Section()), geometry(Section(), half=True)
        self.assertNotIn("rail_l", half)
        self.assertNotIn("rail_face_l", half)
        for k, v in half.items():
            self.assertGreaterEqual(v.bounds[0], 0.0, k)
            if v.geom_type != "LineString":
                self.assertAlmostEqual(v.area, full[k].area / 2 if k not in ("rail_r",) else full[k].area, places=9)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestSlabExact(unittest.TestCase):
    def test_residual_and_fundamental(self):
        lam, d, nc = 1.55, 0.6, 1.444
        eps = (2.1376**2, 2.2111**2, 2.2111**2)
        k0 = 2 * np.pi / lam
        n = slab_exact(eps, nc, d, lam, "TE")
        kap, gam = k0 * np.sqrt(eps[0] - n * n), k0 * np.sqrt(n * n - nc * nc)
        self.assertLess(abs(kap * np.tan(kap * d / 2) - gam) / gam, 1e-9)
        self.assertLess(kap * d / 2, np.pi / 2)
        n = slab_exact(eps, nc, d, lam, "TM")
        kap = k0 * np.sqrt(eps[2] / eps[1]) * np.sqrt(eps[1] - n * n)
        gam = k0 * np.sqrt(n * n - nc * nc)
        self.assertLess(abs(kap / eps[2] * np.tan(kap * d / 2) - gam / nc**2) / (gam / nc**2), 1e-9)

    def test_isotropic_te_above_tm_and_thick_limit(self):
        e = (2.2**2,) * 3
        self.assertGreater(slab_exact(e, 1.444, 0.6, 1.55, "TE"), slab_exact(e, 1.444, 0.6, 1.55, "TM"))
        self.assertGreater(slab_exact(e, 1.444, 20.0, 1.55, "TE"), 2.199)

    def test_group_index_linear(self):
        # n(lam) = a + b lam  ->  ng = a exactly
        self.assertAlmostEqual(group_index(2 + 0.1 * 1.545, 2 + 0.1 * 1.55, 2 + 0.1 * 1.555, 1.55, 0.005), 2.0, places=12)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestSlabPlumbing(unittest.TestCase):
    """Coarse-mesh plumbing check only (m = 4, loose tolerance); the G1 gate itself runs in run_o01.py."""

    def test_te_and_tm_track_exact(self):
        mesh = build_slab_mesh(4.0)
        for pol, te_ok in (("TE", lambda t: t > 0.99), ("TM", lambda t: t < 0.01)):
            r = solve_slab(mesh, 1.55, pol)
            self.assertIsNotNone(r["neff"])
            self.assertLess(abs(r["neff"].real - r["exact"]), 2e-3, pol)
            self.assertTrue(te_ok(r["te_fraction"]), pol)


if __name__ == "__main__":
    unittest.main()
