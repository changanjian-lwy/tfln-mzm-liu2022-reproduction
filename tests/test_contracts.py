import copy
import json
from pathlib import Path
import unittest
from tfln_mzm.contracts import preflight

BASE = json.loads((Path(__file__).resolve().parents[1]/"configs/paper_baseline.json").read_text())

class TestBaselineGuards(unittest.TestCase):
    def test_missing_fem_data_is_visible_and_not_a_pass(self):
        result = preflight(BASE)
        self.assertEqual(set(result["missing"]), {"characteristic_impedance", "rf_attenuation", "microwave_index"})
        self.assertEqual(result["status"], "BLOCKED_BY_MISSING_DATA")
        self.assertFalse(result["physics_validated"])

    def test_units_source_unknown_and_scope_guards(self):
        mutations = [("length","unit","mm"), ("source_voltage","source",""),
                     ("characteristic_impedance","value",50),
                     ("source_impedance","evidence","sensitivity_only"),
                     ("length","value",0), ("source_voltage","value",float('nan'))]
        for name,key,value in mutations:
            cfg=copy.deepcopy(BASE)
            cfg["parameters"][name][key]=value
            with self.subTest(name=name,key=key), self.assertRaises(ValueError):
                preflight(cfg)

    def test_assumption_provider_does_not_earn_physics_validation(self):
        cfg=copy.deepcopy(BASE)
        for name,value in [("characteristic_impedance",50),("rf_attenuation",0),("microwave_index",2.25)]:
            cfg["parameters"][name].update(value=value,evidence="assumption",source="Explicit ideal-limit test")
        result=preflight(cfg)
        self.assertEqual(result["status"], "INPUTS_READY_NOT_VALIDATED")
        self.assertFalse(result["physics_validated"])
        self.assertIn("rf_attenuation", result["nonpaper_inputs"])
