from __future__ import annotations

import itertools
import unittest

from research import sample_projective_lines as sampler
from research import search_e8_a2_target_neighbor_bridge as search


class ProjectiveLineSamplerTests(unittest.TestCase):
    def test_toy_isotropic_spaces_match_exhaustive_enumeration(self):
        gram = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        for prime in (2, 5):
            expected = {
                sampler.normalize_projective(list(vector), prime)
                for vector in itertools.product(range(prime), repeat=3)
                if any(vector) and sampler.q_value(gram, list(vector)) % prime == 0
            }
            window = sampler.sample_window(gram, prime, "toy", 0, len(expected), 5000)
            selected = {tuple(line["vector"]) for line in window["selected_lines"]}
            self.assertEqual(selected, expected)
            self.assertEqual(len(selected), len(window["selected_lines"]))

    def test_offset_is_a_slice_of_the_same_accepted_stream(self):
        gram = search.FROZEN_ROOTLESS_GRAM
        full = sampler.sample_window(gram, 5, "slice", 0, 64, 4096)
        suffix = sampler.sample_window(gram, 5, "slice", 32, 32, 4096)
        self.assertEqual(full["selected_lines"][32:], suffix["selected_lines"])
        self.assertEqual(full["counters"], suffix["counters"])

    def test_determinism_and_seed_separation(self):
        args = (search.target_essential_lattice(), 5, "rank32-20260907-v1", 0, 32, 4096)
        first = sampler.sample_window(*args)
        self.assertEqual(first, sampler.sample_window(*args))
        other = sampler.sample_window(args[0], 5, "different-seed", 0, 32, 4096)
        self.assertNotEqual(first["selected_lines"], other["selected_lines"])
        core = {key: value for key, value in first.items() if key != "window_sha256"}
        self.assertEqual(first["window_sha256"], sampler.payload_sha256(core))

    def test_all_coordinates_are_exercised_on_each_actual_initial_lattice(self):
        for gram in (search.target_essential_lattice(), search.transparent_seed_lattice(), search.FROZEN_ROOTLESS_GRAM):
            for prime in (2, 5):
                record = sampler.sample_window(gram, prime, "rank32-20260907-v1", 0, 32, 4096)
                vectors = [line["vector"] for line in record["selected_lines"]]
                self.assertTrue(all(any(vector[column] for vector in vectors) for column in range(17)))
                for vector in vectors:
                    self.assertEqual(sampler.q_value(gram, vector) % prime, 0)
                    self.assertEqual(next(x for x in reversed(vector) if x), 1)
                costs = record["counters"]
                self.assertEqual(costs["raw_draws"], sum(costs[key] for key in (
                    "zero_vectors", "nonisotropic_vectors", "duplicate_isotropic_lines",
                    "accepted_unique_isotropic_lines",
                )))

    def test_normalization_identifies_all_nonzero_scalar_multiples(self):
        vector = [3, 2, 0, 4, 0]
        for multiple in range(1, 5):
            self.assertEqual(
                sampler.normalize_projective([multiple * x for x in vector], 5),
                sampler.normalize_projective(vector, 5),
            )
        self.assertIsNone(sampler.normalize_projective([0, 0], 5))

    def test_quadratic_form_handles_two_without_modular_division(self):
        gram = [[2, 3], [3, 4]]
        self.assertEqual(sampler.q_value(gram, [1, 1]), 6)
        self.assertEqual(sampler.q_value(gram, [1, 1]) % 2, 0)

    def test_no_isotropic_line_exhausts_budget_instead_of_returning_negative(self):
        with self.assertRaisesRegex(sampler.SamplingBudgetError, "no mathematical conclusion"):
            sampler.sample_window([[2]], 5, "anisotropic", 0, 1, 10)

    def test_duplicate_lines_are_charged_and_do_not_satisfy_the_limit(self):
        with self.assertRaisesRegex(sampler.SamplingBudgetError, "duplicate_isotropic_lines"):
            sampler.sample_window([[0]], 2, "one-line", 0, 2, 30)

    def test_invalid_parameters_fail_closed(self):
        for prime, seed, offset, limit, draws in (
            (4, "a", 0, 1, 10), (257, "a", 0, 1, 10),
            (True, "a", 0, 1, 10), (2, "", 0, 1, 10),
            (2, "bad\nseed", 0, 1, 10), (2, "a", -1, 1, 10),
            (2, "a", 0, 0, 10), (2, "a", 5, 10, 10),
            (2, "a", 0, 1, 1000001),
        ):
            with self.assertRaises(ValueError):
                sampler.sample_window([[2]], prime, seed, offset, limit, draws)

    def test_invalid_grams_fail_closed(self):
        for gram in ([], [[1]], [[2, 1], [0, 2]], [[2, 0]], [[True]], [[2.0]]):
            with self.assertRaises(ValueError):
                sampler.sample_window(gram, 2, "gram", 0, 1, 10)

    @unittest.skipUnless(search.SAGE_AVAILABLE, "requires pinned Sage neighbor implementation")
    def test_hash_sampling_real_sage_neighbor_api(self):
        period = search.Matrix(search.ZZ, search.PERIOD_LATTICE)
        genus = search.IntegralLattice(period).discriminant_group().genus((17, 0))
        for side, gram in (
            ("origin", search.target_essential_lattice()),
            ("transparent", search.transparent_seed_lattice()),
            ("rootless", search.FROZEN_ROOTLESS_GRAM),
        ):
            parent = search.make_initial_state(side, "sampler API test", gram, genus)
            for prime in (2, 5):
                window = sampler.sample_window(gram, prime, "rank32-20260907-v1", 0, 2, 4096)
                for line in window["selected_lines"]:
                    vector = search.sage_vector(search.ZZ, line["vector"])
                    form = search.QuadraticForm(search.ZZ, parent["gram_object"])
                    self.assertEqual(int(form(vector)), sampler.q_value(gram, line["vector"]))
                    built, classification = search.call_neighbor_builder_fail_closed(
                        prime, lambda: search.build_neighbor(
                            parent, prime, line["ordinal_one_based"], vector, genus,
                        ),
                    )
                    if classification is None:
                        self.assertIsNotNone(built)
                        self.assertTrue(built[1]["p_neighbor_relation"]["gram_identity_verified"])
                    else:
                        self.assertEqual(prime, 2)


if __name__ == "__main__":
    unittest.main()
