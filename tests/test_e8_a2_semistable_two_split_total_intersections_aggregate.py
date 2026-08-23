from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "certificates"
    / "e8_a2_semistable_two_split_total_intersections_aggregate.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class E8A2TwoSplitTotalIntersectionAggregateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    def test_source_certificates_are_hash_locked(self) -> None:
        for source in self.payload["source_certificates"]:
            self.assertEqual(sha256(ROOT / source["path"]), source["sha256"])
        replay = self.payload["supplemental_boundary_replay"]
        self.assertEqual(
            sha256(ROOT / replay["source"]), replay["source_sha256"]
        )
        self.assertEqual(
            sha256(ROOT / replay["shared_arithmetic_source"]),
            replay["shared_arithmetic_source_sha256"],
        )

    def test_aggregate_is_derived_from_locked_sources(self) -> None:
        by_role = {
            source["role"]: json.loads(
                (ROOT / source["path"]).read_text(encoding="utf-8")
            )
            for source in self.payload["source_certificates"]
        }
        dense_original = by_role["dense_chart_original_enumeration"]
        dense_corrected = by_role[
            "dense_chart_corrected_total_intersections"
        ]
        boundaries = by_role[
            "three_intrinsic_boundary_charts_original_enumeration"
        ]
        supplement = self.payload["supplemental_boundary_replay"]["summary"]
        aggregate = self.payload["aggregate"]

        self.assertEqual(
            aggregate["raw_surface_parameter_tuples"],
            dense_original["totals"]["raw_parameter_tuples"]
            + boundaries["totals"]["raw_parameter_tuples"],
        )
        self.assertEqual(
            aggregate["rational_chart_tuples"],
            dense_original["totals"]["rational_chart_tuples"]
            + boundaries["totals"]["rational_chart_tuples"],
        )
        self.assertEqual(
            aggregate["kodaira_open_surfaces"],
            dense_original["totals"]["kodaira_surfaces"]
            + boundaries["totals"]["kodaira_surfaces"],
        )
        old_p3_tests = (
            dense_original["totals"]["P3_tests_after_eligible_pairs"]
            + boundaries["totals"]["P3_tests_after_target_pairs"]
        )
        self.assertEqual(
            aggregate["old_simple_infinity_P3_tests"], old_p3_tests
        )
        corrected_p3_tests = (
            dense_corrected["totals"]["p3_tests"]
            + boundaries["totals"]["P3_tests_after_target_pairs"]
            + supplement["p3_tests"]
        )
        self.assertEqual(
            aggregate["corrected_P3_tests"], corrected_p3_tests
        )
        p1_p2_tests = (
            dense_original["totals"]["section_tests"]
            - dense_original["totals"]["P3_tests_after_eligible_pairs"]
            + boundaries["totals"]["section_tests"]
            - boundaries["totals"]["P3_tests_after_target_pairs"]
        )
        self.assertEqual(aggregate["P1_P2_tests"], p1_p2_tests)
        self.assertEqual(
            aggregate["corrected_all_section_tests"],
            p1_p2_tests + corrected_p3_tests,
        )
        boundary_safe_pairs = boundaries["totals"][
            "P1_P2_gcd_degree_counts"
        ]["2"]
        self.assertEqual(
            aggregate["eligible_total_two_P1_P2_pairs"],
            dense_corrected["totals"]["eligible_total_two_pairs"]
            + boundary_safe_pairs
            + supplement["new_collision_pairs"],
        )
        self.assertEqual(
            aggregate["infinity_collision_pairs"],
            dense_corrected["totals"]["new_collision_pairs"]
            + supplement["new_collision_pairs"],
        )
        self.assertEqual(
            aggregate["retained_sections"],
            dense_original["totals"]["retained_sections"]
            + boundaries["totals"]["retained_sections"],
        )
        self.assertEqual(
            aggregate["net_P3_test_increase"],
            aggregate["new_tests_on_genuine_collision_pairs"]
            - aggregate["removed_tests_from_total_three_pairs"],
        )
        self.assertEqual(aggregate["P3_sections"], 0)
        self.assertEqual(aggregate["complete_target_triples"], 0)
        self.assertFalse(self.payload["claim_boundary"]["rank31_curve_found"])

    def test_four_boundary_records_are_exact_and_distinct(self) -> None:
        replay = self.payload["supplemental_boundary_replay"]
        records = replay["records"]
        self.assertEqual(len(records), 4)
        self.assertEqual(
            [(item["parameters"]["q"], item["parameters"]["W"])
             for item in records],
            [(2, 2), (4, 3), (7, 3), (9, 2)],
        )
        self.assertTrue(
            all(
                item["parameters"]["lambda"] == 9
                and item["affine_degree"] == 1
                and item["infinity_degree"] == 1
                for item in records
            )
        )
        self.assertEqual(
            len(
                {
                    (
                        tuple(item["P1"]["U"]), tuple(item["P1"]["V"]),
                        tuple(item["P2"]["U"]), tuple(item["P2"]["V"]),
                    )
                    for item in records
                }
            ),
            4,
        )

    @unittest.skipUnless(shutil.which("cc"), "C compiler unavailable")
    def test_strict_boundary_replay(self) -> None:
        replay = self.payload["supplemental_boundary_replay"]
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "replay"
            subprocess.run(
                [
                    "cc", "-std=c11", "-O3", "-Wall", "-Wextra",
                    "-Wconversion", "-Wshadow", "-pedantic", "-Werror",
                    "-I", str(ROOT / "research"), "-o", str(binary),
                    str(ROOT / replay["source"]),
                ],
                check=True,
                cwd=ROOT,
            )
            completed = subprocess.run(
                [str(binary)], check=True, cwd=ROOT,
                capture_output=True, text=True,
            )
        self.assertEqual(
            hashlib.sha256(completed.stdout.encode("utf-8")).hexdigest(),
            replay["stdout_sha256"],
        )
        rows = [json.loads(line) for line in completed.stdout.splitlines()]
        self.assertEqual(len(rows), 5)
        self.assertTrue(
            all(row["p3_sections"] == 0 for row in rows[:-1])
        )
        self.assertEqual(rows[-1], {
            "event": "summary",
            "prime": 11,
            "records_verified": 4,
            "new_collision_pairs": 4,
            "p3_tests": 644204,
            "p3_sections": 0,
            "total_gram_triples": 0,
        })


if __name__ == "__main__":
    unittest.main()
