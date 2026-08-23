from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "research"
    / "search_e8_a2_semistable_two_split_total_intersections.c"
)
SHARED_SOURCE = ROOT / "research" / "search_e8_a2_semistable_target.c"
CERTIFICATE = (
    ROOT
    / "certificates"
    / "e8_a2_semistable_two_split_total_intersections_small_primes.json"
)


class E8A2TwoSplitTotalIntersectionTests(unittest.TestCase):
    def test_certificate_hashes_and_totals(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(
            hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            payload["engine"]["source_sha256"],
        )
        self.assertEqual(
            hashlib.sha256(SHARED_SOURCE.read_bytes()).hexdigest(),
            payload["engine"]["shared_arithmetic_engine_sha256"],
        )
        runs = payload["runs"]
        self.assertEqual(
            payload["totals"]["eligible_total_two_pairs"],
            sum(run["eligible_total_two_pairs"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["new_collision_pairs"],
            sum(run["new_collision_pairs"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["p3_tests"],
            sum(run["p3_tests"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["enumerated_p1_p2_pairs"],
            sum(sum(run["pair_strata"].values()) for run in runs),
        )
        self.assertEqual(payload["totals"]["new_collision_pairs"], 6)
        self.assertTrue(all(run["total_gram_triples"] == 0 for run in runs))
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
                [str(binary), "5"],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertEqual(
            hashlib.sha256(completed.stdout.encode("utf-8")).hexdigest(),
            expected["stdout_sha256"],
        )
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["eligible_total_two_pairs"], 0)
        self.assertEqual(summary["total_gram_triples"], 0)


if __name__ == "__main__":
    unittest.main()
