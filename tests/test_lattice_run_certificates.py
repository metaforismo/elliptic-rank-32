from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATES = (
    ROOT / "certificates" / "rank17_exact_neighbor_chain_run_32673843229.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_32674002260.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_32677113429.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_34161893470.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_34164507132.json",
    ROOT / "certificates" / "e8_a2_target_neighbor_bridge_run_34174072896.json",
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
        for path, moves, rejections, isometry_checks in (
            (CERTIFICATES[2], 7868, 4, 0),
            (CERTIFICATES[3], 7872, 0, 32),
        ):
            record = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(record["search"]["counters"]["successful_moves"], moves)
            self.assertEqual(record["search"]["total_classified_p2_boundary_failures"], rejections)
            self.assertEqual(record["search"]["isometry_comparisons"]["qfisom_checked"], isometry_checks)
            self.assertTrue(record["search"]["source_checkout"]["git_head_verified"])
            self.assertFalse(record["claim_boundary"]["rank32_curve_found"])
            self.assertFalse(record["claim_boundary"]["exact_lattice_bridge_found"])

    def test_positive_bridge_is_lattice_only_and_sampling_is_replayed(self) -> None:
        record = json.loads(CERTIFICATES[4].read_text())
        self.assertEqual(record["search"]["status"], "bridge_found")
        self.assertEqual(record["search"]["counters"]["successful_moves"], 4595)
        self.assertEqual(record["search"]["isometry_comparisons"]["qfisom_checked"], 5439)
        self.assertTrue(record["claim_boundary"]["exact_lattice_bridge_found"])
        self.assertFalse(record["claim_boundary"]["rank32_curve_found"])
        self.assertEqual(record["search"]["bridge"]["neighbor_step_count"], 7)
        self.assertTrue(record["search"]["bridge"]["complete_bridge_composite_verified"])
        self.assertTrue(record["search"]["sampling_audit"]["all_recorded_draw_streams_replayed"])
        self.assertEqual(record["search"]["sampling_audit"]["windows_replayed"], 144)
        self.assertEqual(record["search"]["sampling_audit"]["sampling_costs"]["raw_draws"], 15976)

    def test_frontier_points_to_all_audits_and_remains_unsolved(self) -> None:
        frontier = json.loads((ROOT / "STATUS.frontier.json").read_text(encoding="utf-8"))
        supporting = set(frontier["newest_result"]["supporting_certificates"])
        for path in CERTIFICATES:
            self.assertIn(str(path.relative_to(ROOT)), supporting)
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

    def test_all_p2_bridge_is_transparent_not_rootless_and_sampling_is_replayed(self) -> None:
        record = json.loads(CERTIFICATES[5].read_text())
        self.assertEqual(record["search"]["status"], "bridge_found")
        self.assertEqual(record["search"]["counters"]["successful_moves"], 4986)
        self.assertEqual(record["search"]["counters"]["failed_moves"], 0)
        self.assertEqual(record["search"]["isometry_comparisons"]["qfisom_checked"], 6509)
        self.assertEqual(record["search"]["bridge"]["endpoint"], "transparent")
        self.assertEqual(record["search"]["bridge"]["neighbor_step_count"], 15)
        self.assertTrue(record["search"]["bridge"]["complete_bridge_composite_verified"])
        sampling = record["search"]["sampling_audit"]
        self.assertTrue(sampling["all_recorded_draw_streams_replayed"])
        self.assertEqual(sampling["windows_replayed"], 156)
        self.assertEqual(sampling["sampling_costs"]["raw_draws"], 10175)
        self.assertEqual(sampling["selected_lines_not_attempted_after_terminal_stop"], 6)
        bridge = json.loads((ROOT / "certificates/e8_a2_transparent_p2_bridge_34174072896/exact-bridge.json").read_text())
        moves = bridge["origin_forward_moves"] + bridge["endpoint_forward_moves_to_invert"]
        self.assertEqual(len(moves), 15)
        self.assertTrue(all(move["prime"] == 2 for move in moves))
        self.assertFalse(record["claim_boundary"]["rank32_curve_found"])
        self.assertFalse(record["claim_boundary"]["explicit_K3_fibration_switch_completed"])


if __name__ == "__main__":
    unittest.main()
