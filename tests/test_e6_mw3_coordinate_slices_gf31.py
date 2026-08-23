import json
import math
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research"))

from replay_e6_mw3_p2_split_cores import canonical_sha256
from replay_e6_mw3_degenerate_pivots_gf31 import (
    compute_certificate as compute_degenerate,
)
from summarize_e6_mw3_coordinate_slices_gf31 import compute as compute_summary


SLICE_FILES = [
    "e6_mw3_core_slice_r1_s1_x1_gf31.json",
    "e6_mw3_core_slice_r2_s1_x1_gf31.json",
    "e6_mw3_core_slice_r3_s1_x1_gf31.json",
    "e6_mw3_core_slice_r1_s2_x1_gf31.json",
    "e6_mw3_core_slice_r1_s1_x2_gf31.json",
]


def load_certificate(filename):
    return json.loads((ROOT / "certificates" / filename).read_text())


def assert_valid_hash(testcase, payload):
    stored = payload["certificate_sha256"]
    unhashed = dict(payload)
    del unhashed["certificate_sha256"]
    testcase.assertEqual(canonical_sha256(unhashed), stored)


class E6MW3CoordinateSliceTests(unittest.TestCase):
    def test_original_core_equation_exhaustion_counts(self):
        payload = load_certificate("e6_mw3_core_exhaustion_gf31.json")
        assert_valid_hash(self, payload)
        self.assertEqual(payload["counts"]["quadruples_examined"], 893730)
        self.assertEqual(payload["counts"]["simultaneous_core_solutions"], 1853)
        self.assertEqual(
            payload["counts"]["valid_surfaces_with_at_least_two_I2_roots"], 2
        )

    def test_five_slice_certificates(self):
        payloads = [load_certificate(filename) for filename in SLICE_FILES]
        for payload in payloads:
            assert_valid_hash(self, payload)
        self.assertEqual(
            sum(item["counts"]["quadruples_examined"] for item in payloads),
            4468650,
        )
        self.assertEqual(
            sum(
                item["counts"]["valid_surfaces_with_at_least_two_I2_roots"]
                for item in payloads
            ),
            12,
        )

    def test_relation_certificate_has_primitive_dependency(self):
        payload = load_certificate(
            "e6_mw3_candidate_r3_s1_x1_core_1_0_14_23_relations.json"
        )
        assert_valid_hash(self, payload)
        primitives = set()
        for relation in payload["relations"]:
            divisor = math.gcd(*(abs(value) for value in relation))
            primitive = tuple(value // divisor for value in relation)
            if next(value for value in primitive if value) < 0:
                primitive = tuple(-value for value in primitive)
            primitives.add(primitive)
        self.assertEqual(primitives, {(1, 0, 2)})

    def test_degenerate_pivot_planes_recompute_and_match_generic_counts(self):
        stored = load_certificate("e6_mw3_degenerate_pivots_gf31.json")
        assert_valid_hash(self, stored)
        self.assertEqual(compute_degenerate(), stored)
        self.assertEqual(
            stored["totals"]["degenerate_core_triples_examined"], 172980
        )
        self.assertEqual(stored["totals"]["simultaneous_core_solutions"], 354)
        self.assertEqual(stored["totals"]["valid_two_I2_surfaces"], 0)

        generic = [
            load_certificate("e6_mw3_core_exhaustion_gf31.json"),
            *(load_certificate(filename) for filename in SLICE_FILES),
        ]
        generic_counts = [
            item["counts"]["solutions_on_degenerate_a1_pivot_plane"]
            for item in generic
        ]
        degenerate_counts = [
            item["counts"]["simultaneous_core_solutions"]
            for item in stored["slice_records"]
        ]
        self.assertEqual(degenerate_counts, generic_counts)

    def test_summary_recomputes(self):
        stored = load_certificate("e6_mw3_coordinate_slices_gf31_summary.json")
        assert_valid_hash(self, stored)
        recomputed = compute_summary()
        self.assertEqual(recomputed, stored)
        self.assertEqual(
            stored["degenerate_pivot_totals"]["valid_two_I2_surfaces"], 0
        )
        self.assertEqual(
            stored["totals"]["independent_rank3_seeds_after_exact_relation_filter"],
            0,
        )


if __name__ == "__main__":
    unittest.main()
