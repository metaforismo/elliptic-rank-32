from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "search_e8_a2_semistable_two_split_boundaries.c"
SHARED_SOURCE = ROOT / "research" / "search_e8_a2_semistable_target.c"
CERTIFICATE = (
    ROOT
    / "certificates"
    / "e8_a2_semistable_two_split_boundary_search_small_primes.json"
)


class E8A2SemistableTwoSplitBoundarySearchTests(unittest.TestCase):
    def test_certificate_matches_sources_and_nested_totals(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(
            hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            payload["engine"]["source_sha256"],
        )
        self.assertEqual(
            hashlib.sha256(SHARED_SOURCE.read_bytes()).hexdigest(),
            payload["engine"]["shared_arithmetic_engine_sha256"],
        )
        charts = [chart for run in payload["runs"] for chart in run["charts"].values()]
        self.assertEqual(payload["totals"]["raw_parameter_tuples"], sum(item["raw"] for item in charts))
        self.assertEqual(
            payload["totals"]["rational_chart_tuples"],
            sum(item["rational"] for item in charts),
        )
        self.assertEqual(payload["totals"]["kodaira_surfaces"], sum(item["kodaira"] for item in charts))
        self.assertEqual(
            payload["totals"]["section_tests"],
            sum(sum(item["tests"].values()) for item in charts),
        )
        self.assertEqual(
            payload["totals"]["retained_sections"],
            sum(sum(item["sections"].values()) for item in charts),
        )
        self.assertEqual(
            payload["totals"]["P3_tests_after_target_pairs"],
            sum(item["tests"]["P3"] for item in charts),
        )
        for degree in ["minus_one", "0", "1", "2", "3", "4", "5", "other"]:
            self.assertEqual(
                payload["totals"]["P1_P2_gcd_degree_counts"][degree],
                sum(item["gcd_degrees"][degree] for item in charts),
            )
        self.assertEqual(len(payload["representative_target_pair_seeds_gf13"]), 3)
        self.assertFalse(payload["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(shutil.which("cc"), "C compiler unavailable")
    def test_strict_build_and_gf5_replay(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        expected = next(run for run in payload["runs"] if run["prime"] == 5)
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "search"
            subprocess.run(
                [
                    "cc", "-std=c11", "-O3", "-Wall", "-Wextra",
                    "-Wconversion", "-Wshadow", "-pedantic", "-Werror",
                    "-o", str(binary), str(SOURCE),
                ],
                check=True,
                cwd=ROOT,
            )
            completed = subprocess.run(
                [str(binary), "5", "0"],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertEqual(
            hashlib.sha256(completed.stdout.encode("utf-8")).hexdigest(),
            expected["stdout_sha256"],
        )
        summaries = [json.loads(line) for line in completed.stdout.splitlines()[1:]]
        self.assertEqual(
            [item["chart"] for item in summaries],
            ["s0_zero", "s_lambda_zero", "node1_beta_gamma_zero"],
        )
        self.assertTrue(all(item["complete_gram_triples"] == 0 for item in summaries))


if __name__ == "__main__":
    unittest.main()
