from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "search_e8_a2_semistable_two_split_pair_first.c"
SHARED_SOURCE = ROOT / "research" / "search_e8_a2_semistable_target.c"
CERTIFICATE = (
    ROOT
    / "certificates"
    / "e8_a2_semistable_two_split_pair_first_large_primes.json"
)


class E8A2SemistableTwoSplitPairFirstTests(unittest.TestCase):
    def test_certificate_matches_sources_sampling_and_totals(self) -> None:
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
        for run in runs:
            self.assertLess(run["sampled_incidences"], run["incidence_universe"])
            self.assertLess(run["offset"], run["incidence_universe"])
            self.assertEqual(math.gcd(run["step"], run["incidence_universe"]), 1)
            self.assertEqual(
                run["tests"]["P3"],
                run["eligible_pairs"] * run["prime"] ** 5,
            )
            self.assertEqual(run["sections"]["P3"], 0)
            self.assertEqual(run["complete_gram_triples"], 0)
            self.assertEqual(
                run["affine_degree_two_pairs"],
                run["P1_P2_gcd_degree_counts"]["2"],
            )
            self.assertEqual(
                run["eligible_pairs"],
                run["P1_P2_total_degree_counts"]["2"],
            )
            self.assertEqual(
                run["sections"]["P2"],
                sum(run["P1_P2_infinity_degree_counts"].values()),
            )
            self.assertEqual(run["total_degree_two_with_positive_infinity"], 0)
        self.assertEqual(
            payload["totals"]["sampled_incidences"],
            sum(run["sampled_incidences"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["rational_chart_incidences"],
            sum(run["rational_chart_incidences"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["kodaira_incidences"],
            sum(run["kodaira_incidences"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["section_tests"],
            sum(sum(run["tests"].values()) for run in runs),
        )
        self.assertEqual(
            payload["totals"]["P3_tests_after_eligible_pairs"],
            sum(run["tests"]["P3"] for run in runs),
        )
        for section in ["P1", "P2", "P3"]:
            self.assertEqual(
                payload["totals"]["sections"][section],
                sum(run["sections"][section] for run in runs),
            )
        for field in [
            "P1_P2_gcd_degree_counts",
            "P1_P2_infinity_degree_counts",
            "P1_P2_total_degree_counts",
        ]:
            for degree in [
                "minus_one", "0", "1", "2", "3", "4", "5", "other"
            ]:
                self.assertEqual(
                    payload["totals"][field][degree],
                    sum(run[field][degree] for run in runs),
                )
        self.assertEqual(
            payload["totals"]["affine_degree_two_pairs"],
            sum(run["affine_degree_two_pairs"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["eligible_pairs"],
            sum(run["eligible_pairs"] for run in runs),
        )
        self.assertEqual(
            payload["totals"]["total_degree_two_with_positive_infinity"],
            sum(
                run["total_degree_two_with_positive_infinity"] for run in runs
            ),
        )
        pairs = payload["previous_affine_degree_two_pairs"]
        self.assertEqual(len(pairs), 3)
        self.assertEqual([pair["prime"] for pair in pairs], [17, 17, 19])
        for pair in pairs:
            self.assertEqual(
                pair["intersection"],
                {"affine": 2, "at_infinity": 0, "total": 2},
            )
        self.assertFalse(payload["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(shutil.which("cc"), "C compiler unavailable")
    def test_strict_build_and_exhaustive_gf5_replay(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        replay = next(
            item for item in payload["validation_replays"] if item["prime"] == 5
        )
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
                [str(binary), "5", "24000", "20260824", "0"],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertEqual(
            hashlib.sha256(completed.stdout.encode("utf-8")).hexdigest(),
            replay["stdout_sha256"],
        )
        start, summary = [json.loads(line) for line in completed.stdout.splitlines()]
        self.assertTrue(start["exhaustive_first_stage"])
        self.assertEqual(
            start["intersection_mode"], "affine gcd plus infinity order"
        )
        self.assertEqual(summary["tests"]["P1"], 2000)
        self.assertEqual(summary["sections"]["P1"], 4)
        self.assertEqual(summary["tests"]["P2"], 2500)
        self.assertEqual(summary["eligible_pairs"], 0)
        self.assertEqual(summary["P1_P2_total_degree_counts"]["2"], 0)


if __name__ == "__main__":
    unittest.main()
