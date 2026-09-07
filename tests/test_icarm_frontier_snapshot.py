from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import math
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "baseline/verify_icarm_frontier_20260907.py"
spec = importlib.util.spec_from_file_location("icarm_frontier_snapshot", SCRIPT)
assert spec is not None and spec.loader is not None
MODULE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(MODULE)


class ICARMFrontierSnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(MODULE.INPUT.read_text())
        cls.helpers, cls.group = MODULE.load_helpers()
        cls.fresh = MODULE.compute_certificate(cls.data)

    def test_exact_certificate_replays(self):
        self.assertEqual(self.fresh, json.loads(MODULE.CERTIFICATE.read_text()))
        core = dict(self.fresh)
        digest = core.pop("certificate_sha256")
        self.assertEqual(digest, self.group.canonical_sha256(core))
        for path, expected in self.fresh["implementation_sha256"].items():
            self.assertEqual(hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), expected)

    def test_all_four_lower_bounds_and_no_upper_bound(self):
        expected = [(273, 30, 30), (302, 31, 34), (398, 30, 32), (582, 30, 32)]
        actual = []
        for record in self.fresh["curves"]:
            matrix = record["binary_matrix"]
            actual.append((record["id"], record["unconditional_lower_bound"], matrix["column_count"]))
            self.assertEqual(matrix["rank"], matrix["row_count"])
            self.assertEqual(len(matrix["pivot_columns_zero_based"]), matrix["rank"])
            self.assertTrue(record["point_membership_both_models"])
            self.assertFalse(record["exact_rank_upper_bound_proved"])
            first, second = record["torsion_reductions"]
            self.assertEqual(math.gcd(first["group_order"], second["group_order"]), 1)
            self.assertNotEqual(second["group_order"] % first["prime"], 0)
            self.assertNotEqual(first["group_order"] % second["prime"], 0)
        self.assertEqual(actual, expected)
        self.assertFalse(self.fresh["rank32_curve_found"])

    def test_rank31_matrix_agrees_with_original_certificate(self):
        old = json.loads((ROOT / "baseline/rank31_icarm302_mod2_certificate.json").read_text())
        new = next(c for c in self.fresh["curves"] if c["id"] == 302)
        self.assertEqual(old["binary_matrix"]["rows"], new["binary_matrix"]["rows"])
        self.assertEqual(old["torsion_certificate"]["reductions"], new["torsion_reductions"])

    def test_census_is_dated_and_j_invariants_are_distinct(self):
        self.assertEqual(self.fresh["source_curve_count"], 631)
        self.assertEqual(self.fresh["source_maximum_rank_lower_bound"], 31)
        self.assertTrue(self.fresh["all_j_invariants_pairwise_distinct"])
        self.assertIn("not nonexistence", self.fresh["claim_boundary"])
        self.assertIn("common parameter family", self.fresh["claim_boundary"])

    def test_altered_point_rejected(self):
        record = copy.deepcopy(self.data["curves"][2])
        record["points"][0][1] = str(int(record["points"][0][1]) + 1)
        with self.assertRaisesRegex(AssertionError, "not on the curve"):
            MODULE.certify_curve(record, self.helpers, self.group)

    def test_duplicate_or_negative_pair_rejected(self):
        record = copy.deepcopy(self.data["curves"][2])
        record["points"][1] = list(record["points"][0])
        with self.assertRaisesRegex(AssertionError, "equal x-coordinates"):
            MODULE.certify_curve(record, self.helpers, self.group)

    def test_insufficient_local_rank_is_not_a_certificate(self):
        # p=11,23 suffice for the torsion argument but not for rank 30.
        with patch.object(self.helpers, "search_primes", return_value=(11, 23)):
            with self.assertRaisesRegex(AssertionError, "local image rank 0, expected 30"):
                MODULE.certify_curve(self.data["curves"][2], self.helpers, self.group)

    def test_wrong_discriminant_and_point_count_rejected(self):
        for field, value in (("discriminant", "1"), ("rank_lower_bound", 32)):
            record = copy.deepcopy(self.data["curves"][2])
            record[field] = value
            with self.assertRaises(AssertionError):
                MODULE.certify_curve(record, self.helpers, self.group)

    def test_changed_rank31_baseline_rejected(self):
        data = copy.deepcopy(self.data)
        data["curves"][1]["points"][0][1] = "0"
        with self.assertRaisesRegex(AssertionError, "baseline changed"):
            MODULE.validate_snapshot(data)

    def test_missing_projection_and_wrong_census_rejected(self):
        data = copy.deepcopy(self.data)
        data["curves"].pop()
        with self.assertRaisesRegex(AssertionError, "curve IDs"):
            MODULE.validate_snapshot(data)
        data = copy.deepcopy(self.data)
        data["source_curve_count"] += 1
        with self.assertRaisesRegex(AssertionError, "census mismatch"):
            MODULE.validate_snapshot(data)

    def test_full_source_projection_roundtrip_and_tamper_gate(self):
        # A compact synthetic source tests projection plumbing, not a new census.
        data = copy.deepcopy(self.data)
        source = dict(data["source"], count=4, curves=data["curves"])
        raw = json.dumps(source).encode()
        data["raw_database_sha256"] = hashlib.sha256(raw).hexdigest()
        data["raw_database_bytes"] = len(raw)
        data["source_curve_count"] = 4
        data["source_rank_histogram"] = {"30": 3, "31": 1}
        MODULE.audit_source_database(data, raw)
        with self.assertRaisesRegex(AssertionError, "SHA-256 mismatch"):
            MODULE.audit_source_database(data, raw + b"\n")
        data["curves"][2]["commentary"] = "changed after download"
        with self.assertRaisesRegex(AssertionError, "projection changed"):
            MODULE.audit_source_database(data, raw)


if __name__ == "__main__":
    unittest.main()
