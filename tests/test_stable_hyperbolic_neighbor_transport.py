import copy
import json
import unittest
from fractions import Fraction
from pathlib import Path

from sympy import Matrix

from research import stable_hyperbolic_neighbor_transport as transport


ROOT = Path(__file__).resolve().parents[1]


def independent_stabilization(gram):
    result = Matrix.diag(Matrix([[0, 1], [1, 0]]), -Matrix(gram))
    assert result.shape == (19, 19)
    return result


class StableHyperbolicTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = transport.build_certificate()
        cls.bridge = json.loads((ROOT / transport.BRIDGE_PATH).read_text())

    def test_committed_certificate_recomputes(self):
        committed = json.loads((ROOT / "certificates/stable_hyperbolic_rootless_bridge_34164507132.json").read_text())
        self.assertEqual(committed, self.certificate)

    def test_all_seven_edge_matrices_have_independent_sympy_checks(self):
        for name, chain in (
            ("origin", self.bridge["origin_chain"]),
            ("rootless", self.bridge["endpoint_chain_endpoint_to_meeting"]),
        ):
            for index, record in enumerate(self.certificate["per_edge_transports"][name]):
                matrix = Matrix(record["child_stabilized_basis_in_parent_coordinates"])
                old = independent_stabilization(chain[index]["gram"])
                new = independent_stabilization(chain[index + 1]["gram"])
                self.assertEqual(matrix.T * old * matrix, new)
                self.assertEqual(int(matrix.det()), record["determinant"])
                self.assertEqual(abs(record["determinant"]), 1)
                fiber = Matrix(record["abstract_new_fiber_class"])
                section = Matrix(record["abstract_new_zero_section_class"])
                self.assertEqual((fiber.T * old * fiber)[0], 0)
                self.assertEqual((section.T * old * fiber)[0], 1)
                self.assertEqual((section.T * old * section)[0], -2)
                self.assertEqual(fiber[1], record["prime"])

    def test_end_to_end_matrix_has_independent_sympy_identity_and_inverse(self):
        matrix = Matrix(self.certificate["rootless_initial_to_origin_initial_stable_matrix"])
        old = independent_stabilization(self.bridge["origin_chain"][0]["gram"])
        new = independent_stabilization(self.bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"])
        self.assertEqual(matrix.T * old * matrix, new)
        self.assertEqual(matrix.det(), 1)
        self.assertEqual(self.certificate["stable_matrix_determinant"], 1)
        self.assertTrue(all(value.q == 1 for value in matrix.inv()))

    def test_tampered_or_fractional_stable_matrix_is_rejected(self):
        parent = self.bridge["origin_chain"][0]["gram"]
        child = self.bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
        for increment in (1, Fraction(1, 2)):
            matrix = copy.deepcopy(self.certificate["rootless_initial_to_origin_initial_stable_matrix"])
            matrix[0][0] += increment
            with self.assertRaises(transport.audit.VerificationError):
                transport.verify_stable_matrix(matrix, parent, child)

    def test_additional_p2_and_p5_examples_replay_despite_inconclusive_source_run(self):
        data = json.loads((ROOT / "baseline/stable_hyperbolic_neighbor_examples_20260907.json").read_text())
        self.assertFalse(data["provenance"]["source_run_is_bounded_negative_certificate"])
        for example in data["examples"]:
            record = transport.extend_move(example["move"], example["parent_gram"])
            matrix = Matrix(record["child_stabilized_basis_in_parent_coordinates"])
            self.assertEqual(matrix.T * independent_stabilization(example["parent_gram"]) * matrix,
                             independent_stabilization(example["move"]["child_gram"]))
            self.assertEqual(record["determinant"], -1)

    def test_radical_p2_line_is_not_silently_promoted(self):
        gram = self.bridge["origin_chain"][0]["gram"]
        with self.assertRaisesRegex(transport.audit.VerificationError, "radical line"):
            transport.lift_isotropic_vector(gram, [0] * 16 + [1], 2)

    def test_bezout_witness_is_exact_and_nonprimitive_input_rejected(self):
        values = [0, -6, 10, 15]
        witness = transport.bezout_one(values)
        self.assertEqual(sum(x * y for x, y in zip(values, witness)), 1)
        with self.assertRaisesRegex(transport.audit.VerificationError, "divisibility"):
            transport.bezout_one([0, 6, 10])

    def test_abstract_marking_is_not_geometric_or_rank_evidence(self):
        claims = self.certificate["claim_boundary"]
        self.assertTrue(claims["abstract_U_stabilized_lattice_isometry_constructed"])
        self.assertTrue(claims["all_seven_neighbor_edges_extended_integrally"])
        for key in (
            "effective_or_nef_divisor_classes_certified", "period_or_ample_cone_compatibility_certified",
            "marked_geometric_K3_transport_completed", "riemann_roch_pencils_or_birational_maps_constructed",
            "rational_P3_found", "rank32_curve_found",
        ):
            self.assertFalse(claims[key])


if __name__ == "__main__":
    unittest.main()
