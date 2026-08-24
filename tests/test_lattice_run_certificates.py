from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATES = (
    ROOT / "certificates" / "rank17_exact_neighbor_chain_run_32673843229.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_32674002260.json",
)


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class LatticeRunCertificateTests(unittest.TestCase):
    def test_compact_run_certificates_have_valid_record_hashes(self) -> None:
        for path in CERTIFICATES:
            with self.subTest(path=path.name):
                record = json.loads(path.read_text(encoding="utf-8"))
                digest = record.pop("record_sha256")
                self.assertEqual(digest, canonical_sha256(record))

    def test_run_certificates_do_not_promote_a_rank_record(self) -> None:
        chain = json.loads(CERTIFICATES[0].read_text(encoding="utf-8"))
        search = json.loads(CERTIFICATES[1].read_text(encoding="utf-8"))
        self.assertIs(chain["claim_boundary"]["rank31_curve_found"], False)
        self.assertIs(search["claim_boundary"]["rank31_curve_found"], False)
        self.assertIs(search["claim_boundary"]["rank32_curve_found"], False)
        self.assertIs(search["claim_boundary"]["exact_lattice_bridge_found"], False)
        self.assertEqual(search["search"]["artifact_schema_generation"], "legacy")
        self.assertEqual(
            search["search"]["legacy_exact_p2_boundary_failures"],
            2,
        )
        self.assertEqual(
            search["search"]["total_classified_p2_boundary_failures"], 2
        )
        self.assertEqual(
            search["search"]["legacy_unclassified_construction_errors"], 0
        )
        self.assertTrue(search["search"]["source_checkout"]["git_head_verified"])
        self.assertIn("verification_implementation", search)
        self.assertIn(
            "Only the 4990 successfully constructed moves",
            search["claim_boundary"]["negative_result_scope"],
        )

    def test_frontier_points_to_both_audits_and_remains_unsolved(self) -> None:
        frontier = json.loads((ROOT / "STATUS.frontier.json").read_text(encoding="utf-8"))
        supporting = set(frontier["newest_result"]["supporting_certificates"])
        self.assertIn(
            "certificates/rank17_exact_neighbor_chain_run_32673843229.json",
            supporting,
        )
        self.assertIn(
            "certificates/e8_a2_target_neighbor_bridge_run_32674002260.json",
            supporting,
        )
        self.assertEqual(frontier["authoritative_lower_bound"], 31)
        self.assertIs(frontier["rank31_lower_bound_reproduced"], True)
        self.assertIs(frontier["rank31_exact_rank_claim_unconditional"], False)
        self.assertIs(frontier["rank32_record_claim"], False)
        self.assertIn("rank 32 remains unsolved", frontier["truth_status"])


if __name__ == "__main__":
    unittest.main()
