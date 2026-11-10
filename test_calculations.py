"""Lightweight regression tests for PVT Lab (run: python -m unittest)."""
import math
import unittest

import app


class ThermodynamicRegressionTests(unittest.TestCase):
    def test_textbook_example_13_5_n_butane_rk(self):
        """Smith et al., Example 13.5: n-butane, RK, 500 K, 50 bar."""
        result = app.calculate({
            "gas": "n-butane", "eos": "rk", "mode": "volume", "phase": "vapor",
            "t_k": 500, "p_bar": 50, "t_ref": 298.15, "p_ref_bar": 1,
            "h_ref": 0, "s_ref": 0,
        })["properties"]
        # Textbook reports Z=0.6850, HR=-4505 J/mol, SR=-6.546 J/(mol K).
        # Small differences are expected from rounded table inputs/constants.
        self.assertAlmostEqual(result["Z"], 0.6850, delta=0.001)
        self.assertAlmostEqual(result["H_res_Jmol"], -4505, delta=10)
        self.assertAlmostEqual(result["S_res_JmolK"], -6.546, delta=0.02)

    def test_all_gas_eos_pairs_produce_finite_properties(self):
        for gas_key in app.GASES:
            for eos_key in app.EOS:
                with self.subTest(gas=gas_key, eos=eos_key):
                    props = app.calculate({
                        "gas": gas_key, "eos": eos_key, "mode": "volume", "phase": "stable",
                        "t_k": 350, "p_bar": 10, "t_ref": 298.15, "p_ref_bar": 1,
                        "h_ref": 0, "s_ref": 0,
                    })["properties"]
                    for key in ("Z", "V_Lmol", "phi", "fugacity_bar", "H_Jmol",
                                "S_JmolK", "H_res_Jmol", "S_res_JmolK"):
                        self.assertTrue(math.isfinite(props[key]), f"{key} is not finite")

    def test_three_calculation_modes(self):
        common = {"gas": "methane", "eos": "pr", "t_ref": 298.15,
                  "p_ref_bar": 1, "h_ref": 0, "s_ref": 0}
        by_volume = app.calculate({**common, "mode": "volume", "t_k": 300, "p_bar": 10})["properties"]
        by_pressure = app.calculate({**common, "mode": "pressure", "t_k": 300,
                                     "v_lmol": by_volume["V_Lmol"]})["properties"]
        by_temperature = app.calculate({**common, "mode": "temperature", "p_bar": 10,
                                        "v_lmol": by_volume["V_Lmol"]})["properties"]
        self.assertAlmostEqual(by_pressure["P_bar"], 10, delta=1e-7)
        self.assertAlmostEqual(by_temperature["T_K"], 300, delta=1e-5)

    def test_middle_cubic_root_is_not_selected_as_stable(self):
        gas = app.GASES["n-butane"]
        roots = app.compressibility_roots(300.0, 10.0e5, gas, "pr")
        self.assertEqual(len(roots), 3)
        selected, label = app._select_root(roots, 300.0, 10.0e5, gas, "pr", "stable")
        self.assertIn(selected, (roots[0], roots[-1]))
        self.assertNotEqual(selected, roots[1])
        self.assertIn("outer root", label)

    def test_cp_table_values_at_298_k(self):
        for gas_key, gas in app.GASES.items():
            with self.subTest(gas=gas_key):
                cp_r = app.cp_over_r(298.15, gas)
                self.assertAlmostEqual(cp_r, gas["cp298_R"], delta=0.002)


if __name__ == "__main__":
    unittest.main(verbosity=2)
