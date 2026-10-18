"""Synthetic mathematical checks; no material calculations or VASP required."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from fit_effective_mass import HBAR2_OVER_ME_EV_A2, analyze, quadratic_fit
from make_mass_kpoints import make_mesh


class EffectiveMassTests(unittest.TestCase):
    def setUp(self):
        self.lattice = np.array([[4., 0., 0.], [-2., 3.4641, 0.], [.3, .5, 20.]])
        _, _, self.q, _ = make_mesh(self.lattice, np.array([.1, .2, .5]), .005, 5)
        # Known principal masses rotated away from the Cartesian axes.
        self.masses = np.array([.25, .7, 1.6])
        raw = np.array([[1., 2., -.3], [.2, -.3, 2.], [2., 1., 1.]])
        self.rotation, _ = np.linalg.qr(raw)
        self.hessian = self.rotation @ np.diag(HBAR2_OVER_ME_EV_A2/self.masses) @ self.rotation.T

    def energies(self, hessian, gradient=None):
        energy = 6.1 + .5*np.einsum('ni,ij,nj->n', self.q, hessian, self.q)
        return energy if gradient is None else energy+self.q@gradient

    def test_rotated_anisotropic_electrons(self):
        result, _ = analyze(self.q, self.energies(self.hessian), 'electron')
        self.assertTrue(result['all_curvatures_have_required_sign'])
        self.assertEqual(result['warnings'], [])
        np.testing.assert_allclose(sorted(result['principal_masses_me']), self.masses, rtol=1e-10)
        axes = np.array(result['principal_axes_cartesian_columns'])
        eig = np.array(result['hessian_eigenvalues_eV_A2'])
        np.testing.assert_allclose(axes@np.diag(eig)@axes.T, self.hessian, atol=1e-9)

    def test_holes_do_not_take_absolute_electron_curvature_blindly(self):
        holes, _ = analyze(self.q, self.energies(-self.hessian), 'hole')
        np.testing.assert_allclose(sorted(holes['principal_masses_me']), self.masses, rtol=1e-10)
        wrong, _ = analyze(self.q, self.energies(-self.hessian), 'electron')
        self.assertFalse(wrong['all_curvatures_have_required_sign'])
        self.assertEqual(wrong['principal_masses_me'], [None]*3)
        self.assertTrue(any('WRONG_CURVATURE_SIGN' in x for x in wrong['warnings']))

    def test_detect_off_center_extremum(self):
        gradient = np.array([.1, -.03, .015])
        result, _ = analyze(self.q, self.energies(self.hessian, gradient), 'electron')
        np.testing.assert_allclose(result['gradient_eV_A'], gradient, atol=1e-11)
        np.testing.assert_allclose(result['quadratic_stationary_shift_Ainv'],
                                   -np.linalg.solve(self.hessian, gradient), atol=1e-11)
        self.assertTrue(any('NONZERO_GRADIENT' in x for x in result['warnings']))

    def test_nonparabolic_window_diagnostic(self):
        energy = self.energies(self.hessian) + 2e5*self.q[:, 0]**4
        result, _ = analyze(self.q, energy, 'electron')
        self.assertTrue(result['inner_window_diagnostic']['available'])
        self.assertTrue(any('WINDOW_DEPENDENT_CURVATURE' in x for x in result['warnings']))

    def test_reject_path_only_sampling(self):
        line = np.column_stack([np.linspace(-.02, .02, 21), np.zeros(21), np.zeros(21)])
        with self.assertRaisesRegex(ValueError, 'rank='):
            quadratic_fit(line, line[:, 0]**2)

    def test_cli_mapping_output_and_rejection(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            meta = tmp/'metadata.json'
            energy = tmp/'energies.csv'
            output = tmp/'fit.json'
            records = [{'kpoint_index_1based': i+1, 'offset_cartesian_Ainv': row.tolist()}
                       for i, row in enumerate(self.q)]
            meta.write_text(json.dumps({'center_fractional': [.1, .2, .5], 'points': records}))
            expected = self.energies(self.hessian)
            # Reverse rows deliberately: matching uses explicit indices.
            with energy.open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['kpoint_index_1based', 'energy_eV'])
                writer.writerows([(i+1, expected[i]) for i in reversed(range(len(expected)))])
            script = Path(__file__).with_name('fit_effective_mass.py')
            command = [sys.executable, str(script), '--metadata', str(meta), '--energies', str(energy),
                       '--carrier', 'electron', '--output', str(output)]
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            result = json.loads(output.read_text())
            np.testing.assert_allclose(sorted(result['principal_masses_me']), self.masses, rtol=1e-10)
            self.assertEqual(len(result['point_residuals']), 125)
            again = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)
            self.assertIn('already exists', again.stderr)
            # Missing indexed data must fail, not silently change the sample.
            energy.write_text('kpoint_index_1based,energy_eV\n1,6.1\n')
            bad = subprocess.run(command+['--force'], capture_output=True, text=True)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn('match exactly', bad.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
