from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "search_e8_a2_semistable_target.c"
CERTIFICATE = (
    ROOT / "certificates" / "e8_a2_semistable_target_small_primes.json"
)


class E8A2SemistableTargetSearchTests(unittest.TestCase):
    def test_certificate_matches_source_and_run_totals(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        self.assertEqual(source_hash, payload["engine"]["source_sha256"])
        self.assertEqual(
            payload["totals"]["section_tests"],
            sum(sum(run["tests"].values()) for run in payload["runs"]),
        )
        self.assertEqual(
            payload["totals"]["retained_sections"],
            sum(sum(run["sections"].values()) for run in payload["runs"]),
        )
        self.assertEqual(
            payload["totals"]["kodaira_surfaces"],
            sum(run["kodaira_surfaces"] for run in payload["runs"]),
        )
        self.assertEqual(
            [run["complete_gram_triples"] for run in payload["runs"]],
            [0, 0, 0],
        )
        self.assertFalse(payload["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(shutil.which("cc"), "C compiler unavailable")
    def test_strict_build_and_gf5_replay(self) -> None:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        expected = next(run for run in payload["runs"] if run["prime"] == 5)
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "search"
            subprocess.run(
                [
                    "cc",
                    "-std=c11",
                    "-O3",
                    "-Wall",
                    "-Wextra",
                    "-Wconversion",
                    "-Wshadow",
                    "-pedantic",
                    "-Werror",
                    "-o",
                    str(binary),
                    str(SOURCE),
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
        summary = json.loads(completed.stdout.splitlines()[-1])
        self.assertEqual(summary["kodaira_surfaces"], 16)
        self.assertEqual(summary["complete_gram_triples"], 0)


if __name__ == "__main__":
    unittest.main()
