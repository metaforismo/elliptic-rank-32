import copy
import itertools
import json
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

from sympy import Matrix

from research import reduce_stable_hyperbolic_marking as short
from research import reduce_hyperbolic_fiber_degree as eichler
from research import reduce_hyperbolic_root_reflections as roots
from research import stable_degree_two_bridge as degree_two

ROOT = Path(__file__).resolve().parents[1]
stable = short.stable


def independent_gram(gram):
    return Matrix.diag(Matrix([[0, 1], [1, 0]]), -Matrix(gram))


class MarkingReductionCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bridge = json.loads((ROOT / stable.BRIDGE_PATH).read_text())
        cls.origin = cls.bridge["origin_chain"][0]["gram"]
        cls.endpoint = cls.bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
        cls.short = short.build_certificate()
        cls.eichler = eichler.build_certificate()
        cls.roots = roots.build_certificate()

    def test_all_three_committed_certificates_recompute_exactly(self):
        for filename, fresh in (
            ("short_hyperbolic_rootless_marking_20260908.json", self.short),
            ("reduced_hyperbolic_fiber_degree_20260908.json", self.eichler),
            ("root_reflection_hyperbolic_marking_20260908.json", self.roots),
        ):
            with self.subTest(filename=filename):
                self.assertEqual(json.loads((ROOT / "certificates" / filename).read_text()), fresh)

    def test_seven_short_edge_markings_have_independent_sympy_gram_and_determinant_checks(self):
        count = 0
        for side, chain in (
            ("origin", self.bridge["origin_chain"]),
            ("rootless", self.bridge["endpoint_chain_endpoint_to_meeting"]),
        ):
            for index, record in enumerate(self.short["per_edge_transports"][side]):
                matrix = Matrix(record["matrix"])
                sg = independent_gram(chain[index]["gram"])
                self.assertEqual(matrix.T * sg * matrix, independent_gram(chain[index + 1]["gram"]))
                self.assertEqual(int(matrix.det()), record["determinant"])
                self.assertEqual(abs(record["determinant"]), 1)
                fiber, section = Matrix(record["abstract_fiber"]), Matrix(record["abstract_zero_section"])
                self.assertEqual((fiber.T * sg * fiber)[0], 0)
                self.assertEqual((section.T * sg * fiber)[0], 1)
                self.assertEqual((section.T * sg * section)[0], -2)
                self.assertEqual(fiber[1], record["prime"])
                count += 1
        self.assertEqual(count, 7)

    def test_all_three_end_to_end_matrices_have_independent_sympy_checks(self):
        for record, key in (
            (self.short, "rootless_to_origin_stable_matrix"),
            (self.eichler, "reduced_rootless_to_origin_stable_matrix"),
            (self.roots, "rootless_to_origin_stable_matrix"),
        ):
            matrix = Matrix(record[key])
            self.assertEqual(matrix.T * independent_gram(self.origin) * matrix, independent_gram(self.endpoint))
            self.assertEqual(matrix.det(), 1)
            self.assertEqual(record["determinant"], 1)
            self.assertTrue(all(value.q == 1 for value in matrix.inv()))

    def test_eichler_trace_replays_with_independent_matrices(self):
        matrix = Matrix(self.short["rootless_to_origin_stable_matrix"])
        for event in self.eichler["trace"]:
            gram = Matrix(self.origin if event["side"] == "left" else self.endpoint)
            shift = Matrix(event["shift"])
            action = Matrix.eye(19)
            action[0, 1] = (shift.T * gram * shift)[0] / 2
            for i in range(17):
                action[0, i + 2] = (gram * shift)[i]
                action[i + 2, 1] = shift[i]
            self.assertEqual(action.T * independent_gram(gram.tolist()) * action, independent_gram(gram.tolist()))
            self.assertEqual(action.det(), 1)
            matrix = action * matrix if event["side"] == "left" else matrix * action
            if event["swap_performed"]:
                if event["side"] == "left":
                    matrix.row_swap(0, 1)
                else:
                    matrix.col_swap(0, 1)
            self.assertEqual(matrix[1, 0], event["fiber_intersection_after"])
        self.assertEqual(matrix, Matrix(self.eichler["reduced_rootless_to_origin_stable_matrix"]))

    def test_every_reflection_is_an_independent_integral_involution(self):
        matrix = Matrix(self.eichler["reduced_rootless_to_origin_stable_matrix"])
        for event in self.roots["reflection_trace"]:
            sg = independent_gram(self.origin if event["side"] == "left" else self.endpoint)
            root = Matrix(event["root"])
            self.assertEqual((root.T * sg * root)[0], -2)
            action = Matrix.eye(19) + root * (root.T * sg)
            self.assertEqual(action.T * sg * action, sg)
            self.assertEqual(action * action, Matrix.eye(19))
            self.assertEqual(action.det(), -1)
            self.assertEqual(matrix[1, 0], event["degree_before"])
            fiber = matrix[:, 0] if event["side"] == "left" else matrix.inv()[:, 0]
            self.assertEqual((fiber.T * sg * root)[0], event["fiber_root_pairing"])
            matrix = action * matrix if event["side"] == "left" else matrix * action
            self.assertEqual(matrix[1, 0], event["degree_after"])
            self.assertLess(event["degree_after"], event["degree_before"])
        actual = eichler.replay_trace(
            [[int(x) for x in row] for row in matrix.tolist()],
            self.roots["postprocessing_eichler_trace"], self.origin, self.endpoint,
        )
        self.assertEqual(actual, self.roots["rootless_to_origin_stable_matrix"])
        self.assertEqual(len(self.roots["reflection_trace"]), 22)

    def test_measured_reduction_is_not_global_optimality_or_geometry(self):
        self.assertEqual(self.short["comparison"]["new_fiber_intersection"], 8701436)
        self.assertEqual(self.eichler["comparison"]["reduced_fiber_intersection"], 6050128)
        self.assertEqual(self.roots["comparison"]["reduced_fiber_intersection"], 710298)
        self.assertEqual(self.roots["comparison"]["reduced_maximum_absolute_entry"], 1320133)
        self.assertEqual(self.roots["configuration"]["maximum_root_height"], 64)
        self.assertEqual(self.roots["search_nodes_visited"], 125793)
        self.assertFalse(any(event["node_budget_exhausted"] for event in self.roots["root_searches"]))
        for record in (self.short, self.eichler, self.roots):
            for claim in ("effective_or_nef_divisor_classes_certified", "marked_geometric_K3_transport_completed",
                          "rational_P3_found", "rank32_curve_found"):
                self.assertFalse(record["claim_boundary"][claim])
        self.assertFalse(self.short["claim_boundary"]["global_minimum_marking_or_intersection_proved"])
        self.assertFalse(self.eichler["claim_boundary"]["global_minimum_fiber_intersection_proved"])
        self.assertFalse(self.roots["claim_boundary"]["all_negative_pairing_roots_excluded"])

    def test_tampered_input_hash_or_trace_is_rejected(self):
        for module, record in ((short, self.short), (eichler, self.eichler), (roots, self.roots)):
            with patch.object(module.audit, "load_json", return_value=record):
                with self.assertRaisesRegex(stable.audit.VerificationError, "drift"):
                    module.build_certificate()
        trace = copy.deepcopy(self.eichler["trace"])
        trace[0]["fiber_intersection_after"] += 1
        with self.assertRaises(stable.audit.VerificationError):
            eichler.replay_trace(self.short["rootless_to_origin_stable_matrix"], trace, self.origin, self.endpoint)
        trace = copy.deepcopy(self.roots["reflection_trace"])
        trace[0]["fiber_root_pairing"] += 1
        with self.assertRaises(stable.audit.VerificationError):
            roots.replay_reflections(self.eichler["reduced_rootless_to_origin_stable_matrix"], trace, self.origin, self.endpoint)

    def test_zero_reduction_budgets_are_inconclusive_and_preserve_the_matrix(self):
        initial = self.short["rootless_to_origin_stable_matrix"]
        result, trace, _, termination = eichler.reduce_matrix(initial, self.origin, self.endpoint, max_iterations=0)
        self.assertEqual(result, initial)
        self.assertEqual(trace, [])
        self.assertEqual(termination, "iteration_budget_exhausted")
        result, trace, _, termination = roots.reduce_roots(initial, self.origin, self.endpoint, maximum_reflections=0)
        self.assertEqual(result, initial)
        self.assertEqual(trace, [])
        self.assertEqual(termination, "reflection_budget_exhausted")


class ExactReductionPrimitiveTests(unittest.TestCase):
    def test_gram_schmidt_reconstructs_the_gram_independently(self):
        gram = [[4, 2, 0], [2, 6, 2], [0, 2, 4]]
        mu, diagonal = eichler.gram_schmidt(gram)
        lower = Matrix(mu) + Matrix.eye(3)
        self.assertEqual(lower * Matrix.diag(*diagonal) * lower.T, Matrix(gram))

    def test_invalid_grams_and_search_inputs_fail_closed(self):
        for gram in ([], [[2, 0]], [[2, 1], [0, 2]], [[1]], [[0]], [[-2]], [[True]], [[Fraction(2)]]):
            with self.assertRaises(stable.audit.VerificationError):
                eichler.gram_schmidt(gram)
        for divisor, nodes in ((0, 5), (-1, 5), (True, 5), (2, -1)):
            with self.assertRaises(stable.audit.VerificationError):
                eichler.bounded_shift([[4]], [1], divisor, max_nodes=nodes)
        for prime in (0, 1, 4, 257):
            with self.assertRaises(stable.audit.VerificationError):
                short.short_lifts([[4]], [1], prime)
        with self.assertRaises(stable.audit.VerificationError):
            short.short_lifts([[4]], [0], 2)
        with self.assertRaises(stable.audit.VerificationError):
            short.short_lifts([[2]], [1], 2)
        for root in ([0, 0, 0], [Fraction(1, 2), 1, 0]):
            with self.assertRaises(stable.audit.VerificationError):
                roots.reflection_matrix([[4]], root)

    def test_exact_cvp_small_exhaustive_case_matches_brute_force(self):
        gram, coordinates, divisor = [[10, 2], [2, 12]], [1, 1], 2
        shift, report = eichler.bounded_shift(gram, coordinates, divisor)
        target = Matrix([Fraction(x, divisor) for x in coordinates])
        best = min(((target + Matrix(x)).T * Matrix(gram) * (target + Matrix(x)))[0]
                   for x in itertools.product(range(-3, 4), repeat=2))
        selected = ((target + Matrix(shift)).T * Matrix(gram) * (target + Matrix(shift)))[0]
        self.assertEqual(selected, best)
        self.assertEqual(Fraction(report["selected_squared_norm"]), best)
        self.assertTrue(report["enumeration_finished_without_early_target_or_budget"])
        self.assertFalse(report["strict_norm_below_two_found"])

    def test_cvp_early_target_and_budget_reporting(self):
        _, early = eichler.bounded_shift([[4, 2], [2, 6]], [7, -4], 5)
        self.assertTrue(early["strict_norm_below_two_found"])
        self.assertLess(Fraction(early["selected_squared_norm"]), 2)
        _, bounded = eichler.bounded_shift([[10, 2], [2, 12]], [1, 1], 2, max_nodes=0)
        self.assertEqual(bounded["nodes_visited"], 0)
        self.assertTrue(bounded["node_budget_exhausted"])

    def test_root_search_agrees_with_toy_brute_force_including_strict_boundary(self):
        gram = [[4, 0], [0, 4]]
        for fiber in ([1, 4, 1, 1], [2, 2, 1, 1]):
            for height in (1, 2, 3):
                witnesses = []
                for u in itertools.product(range(-3, 5), repeat=2):
                    q = 2 * sum(x * x for x in u)
                    if q % height != 1 % height:
                        continue
                    root = [(q - 1) // height, height, *u]
                    pairing = (Matrix(fiber).T * independent_gram(gram) * Matrix(root))[0]
                    if pairing < 0 and fiber[1] + pairing * height > 0:
                        witnesses.append(root)
                found, report = roots.decreasing_root(gram, fiber, height, 10000)
                self.assertEqual(found is not None, bool(witnesses))
                self.assertFalse(report["node_budget_exhausted"])
                if found is not None:
                    self.assertIn(found, witnesses)
        found, _ = roots.decreasing_root(gram, [2, 2, 1, 1], 1, 10000)
        self.assertIsNone(found)  # The shortest errors have norm exactly two.

    def test_root_search_budget_and_zero_cusp_boundary_are_not_promoted(self):
        found, report = roots.decreasing_root([[4]], [1, 2, 1], 1, 0)
        self.assertIsNone(found)
        self.assertTrue(report["node_budget_exhausted"])
        with self.assertRaisesRegex(stable.audit.VerificationError, "zero-height cusp"):
            roots.decreasing_root([[4]], [0, 2, 0], 1, 100)


class StableDegreeTwoBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bridge = json.loads((ROOT / degree_two.BRIDGE_PATH).read_text())
        cls.certificate = degree_two.build_certificate()

    def test_committed_certificate_recomputes(self):
        committed = json.loads((ROOT / "certificates/stable_degree_two_bridge_34174072896.json").read_text())
        self.assertEqual(committed, self.certificate)

    def test_all_fifteen_edges_have_independent_sympy_checks(self):
        count = 0
        for side, chain in (
            ("origin", self.bridge["origin_chain"]),
            ("transparent", self.bridge["endpoint_chain_endpoint_to_meeting"]),
        ):
            for index, record in enumerate(self.certificate["per_edge_transports"][side]):
                matrix = Matrix(record["matrix"])
                sg = independent_gram(chain[index]["gram"])
                self.assertEqual(matrix.T * sg * matrix, independent_gram(chain[index + 1]["gram"]))
                self.assertEqual(matrix.det(), record["determinant"])
                self.assertEqual(abs(record["determinant"]), 1)
                fiber, section = Matrix(record["abstract_fiber"]), Matrix(record["abstract_zero_section"])
                self.assertEqual((fiber.T * sg * fiber)[0], 0)
                self.assertEqual((section.T * sg * fiber)[0], 1)
                self.assertEqual((section.T * sg * section)[0], -2)
                self.assertEqual(record["prime"], 2)
                self.assertEqual(fiber[1], 2)
                count += 1
        self.assertEqual(count, 15)

    def test_end_to_end_integral_isometry_is_not_a_geometric_degree_two_map(self):
        matrix = Matrix(self.certificate["transparent_to_origin_stable_matrix"])
        self.assertEqual(matrix.T * independent_gram(self.bridge["origin_chain"][0]["gram"]) * matrix,
                         independent_gram(self.bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]))
        self.assertEqual(abs(matrix.det()), 1)
        self.assertEqual(matrix.det(), self.certificate["determinant"])
        self.assertTrue(all(x.q == 1 for x in matrix.inv()))
        self.assertEqual(matrix[1, 0], self.certificate["end_to_end_fiber_intersection"])
        self.assertNotEqual(matrix[1, 0], 2)
        claims = self.certificate["claim_boundary"]
        self.assertTrue(claims["all_fifteen_edges_extended_integrally"])
        for claim in ("endpoint_is_the_rootless_lattice", "end_to_end_fiber_intersection_is_two",
                      "effective_or_nef_divisor_classes_certified", "period_or_ample_cone_compatibility_certified",
                      "marked_geometric_K3_transport_completed", "riemann_roch_pencils_or_birational_maps_constructed",
                      "rational_P3_found", "rank32_curve_found"):
            self.assertFalse(claims[claim])


if __name__ == "__main__":
    unittest.main()
