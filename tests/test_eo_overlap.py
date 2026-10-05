import importlib.util
import unittest
from collections import OrderedDict

import numpy as np

HAVE_FEMWELL = importlib.util.find_spec("femwell") is not None

if HAVE_FEMWELL:
    from shapely.geometry import box
    from tfln_mzm.eo_overlap import (FieldX, deps_xx, dc_permittivity, element_ex, gamma_a, solve_potential, vpil_a, vpil_b)
    from tfln_mzm.optical_fem import Section, solve_modes


def _mesh(shapes, res):
    from femwell.mesh import mesh_from_OrderedDict
    from skfem.io.meshio import from_meshio
    return from_meshio(mesh_from_OrderedDict(shapes, {}, default_resolution_max=res))


def two_layer(axis, a, d, e1, e2, V=1.0, res=0.25):
    """Plates normal to `axis` at 0 and d, layer 1 on [0, a]; returns (FEM fields, exact fields) per layer."""
    from skfem import ElementTriP0
    if axis == 0:
        shapes = OrderedDict(l1=box(0, 0, a, 2.0), l2=box(a, 0, d, 2.0))
    else:
        shapes = OrderedDict(l1=box(0, 0, 2.0, a), l2=box(0, a, 2.0, d))
    mesh = _mesh(shapes, res)
    f_lo = mesh.facets_satisfying(lambda p: np.abs(p[axis]) < 1e-9, boundaries_only=True)
    f_hi = mesh.facets_satisfying(lambda p: np.abs(p[axis] - d) < 1e-9, boundaries_only=True)
    from skfem import Basis, ElementTriP2
    b = Basis(mesh, ElementTriP2())
    basis, phi = solve_potential(mesh, {"l1": e1, "l2": e2},
                                 [(b.get_dofs(facets=f_lo).all(), 0.0), (b.get_dofs(facets=f_hi).all(), V)])
    b0 = basis.with_element(ElementTriP0())
    grad = basis.interpolate(phi).grad[axis]
    fem = [float(np.mean(grad[mesh.subdomains[n]])) for n in ("l1", "l2")]
    k1, k2 = e1[axis], e2[axis]
    E1 = V / (a + (d - a) * k1 / k2)            # D continuous across the interface
    return fem, [E1, E1 * k1 / k2], b0


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestPotential(unittest.TestCase):
    def test_tensor_orientation(self):
        eps = dc_permittivity(Section(), "clamped", 3.9)
        self.assertEqual(eps["ridge"], (27.9, 44.0))      # x = crystal Z (eps33), y = crystal X (eps11)
        self.assertEqual(dc_permittivity(Section(), "unclamped", 3.9)["slab"], (28.7, 85.0))
        self.assertEqual(eps["substrate"], (4.5, 4.5))
        self.assertEqual(dc_permittivity(Section(substrate="fused"), "clamped", 3.9)["substrate"], (3.8, 3.8))

    def test_two_layer_exact_both_axes(self):
        # anisotropic layers: only the component along the stacking axis may matter
        e1, e2 = (27.9, 44.0), (3.9, 3.9)
        for axis in (0, 1):
            fem, exact, _ = two_layer(axis, 0.7, 2.0, e1, e2)
            for f, x in zip(fem, exact):
                self.assertLess(abs(f / x - 1), 1e-8, axis)
        # swapping the orientation changes the answer, so the test can see a transposed tensor
        fx, _, _ = two_layer(0, 0.7, 2.0, e1, e2)
        fy, _, _ = two_layer(1, 0.7, 2.0, e1, e2)
        self.assertGreater(abs(fx[0] / fy[0] - 1), 0.1)


@unittest.skipUnless(HAVE_FEMWELL, "femwell not installed (see requirements-fem-lock.txt)")
class TestOverlap(unittest.TestCase):
    def test_field_transfer_between_meshes(self):
        """E_x of phi = x^2 / 2 (so E_x = -x) on a coarse mesh, probed at the points of a different mesh."""
        from skfem import Basis, ElementTriP2
        coarse = _mesh(OrderedDict(slab=box(0, 0, 2, 1)), 0.4)
        fine = _mesh(OrderedDict(slab=box(0.1, 0.1, 1.9, 0.9)), 0.13)
        b = Basis(coarse, ElementTriP2())
        f = FieldX(b, 0.5 * b.doflocs[0] ** 2)
        pts = fine.p
        self.assertTrue(np.allclose(f.at(pts), -pts[0], atol=1e-12))

    def test_uniform_field_gives_gamma_one(self):
        from skfem import Basis, ElementTriP2
        mesh = _mesh(OrderedDict(slab=box(-1.5, -1.5, 1.5, 1.5)), 0.3)
        modes, _ = solve_modes(mesh, {"slab": (4.0, 4.0, 4.0)}, 1.55, n_guess=2.0, num_modes=1, intorder=4)
        b = Basis(mesh, ElementTriP2())
        g, V = 3.0, 1.0
        phi = -(V / g) * b.doflocs[0]                     # E_x = V/g everywhere, exact in P2
        field = FieldX(b, phi)
        self.assertAlmostEqual(gamma_a(modes[0], field, g, V), 1.0, places=10)
        ex = element_ex(mesh, field)
        self.assertTrue(np.allclose(ex, V / g, rtol=1e-10))
        d = deps_xx(mesh, ex, 2.0, 31e-6, 10.0)
        self.assertTrue(np.allclose(d, -2.0**4 * 31e-6 * (V / g) * 10.0))

    def test_vpil_formulas(self):
        self.assertAlmostEqual(vpil_a(1.55, 3.0, 2.13756, 31e-6, 1.0), 0.7679, places=4)   # O02 cheaper check
        # convention B with dn/dV equal to convention A's 0.5 ne^3 r33 Gamma / g gives the same VpiL
        dn = 0.5 * 2.13756**3 * 31e-6 * 0.5 / 3.0
        self.assertAlmostEqual(vpil_b(1.55, dn), vpil_a(1.55, 3.0, 2.13756, 31e-6, 0.5), places=12)


if __name__ == "__main__":
    unittest.main()
