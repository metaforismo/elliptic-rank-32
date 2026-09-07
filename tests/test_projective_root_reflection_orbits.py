import itertools
import json
import unittest
from pathlib import Path

from research import projective_root_reflection_orbits as audit


class RootReflectionOrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = audit.build_certificate()

    def test_committed_certificate_recomputes(self):
        path = Path(__file__).resolve().parents[1] / "certificates/initial_f2_root_reflection_orbits.json"
        self.assertEqual(json.loads(path.read_text()), self.certificate)

    def test_exact_initial_counts_and_no_rank_promotion(self):
        records = self.certificate["initial_lattices"]
        self.assertEqual({name: item["nonzero_isotropic_lines"] for name, item in records.items()}, dict.fromkeys(records, 65279))
        self.assertEqual({name: item["orbit_count"] for name, item in records.items()}, {"origin": 95, "transparent": 34, "rootless": 65279})
        self.assertFalse(self.certificate["claim_boundary"]["neighbors_constructed"])
        self.assertFalse(self.certificate["claim_boundary"]["rank32_curve_found"])

    def test_toy_mask_evaluator_matches_integer_quadratic_form(self):
        gram = [[2, -1, 0], [-1, 4, 1], [0, 1, 6]]
        for vector in itertools.product((0, 1), repeat=3):
            mask = sum(value << i for i, value in enumerate(vector))
            self.assertEqual(audit.q_mod2(mask, audit.quadratic_masks(gram)), audit.sampler.q_value(gram, list(vector)) % 2)

    def test_reflection_action_matches_exact_integer_matrix(self):
        gram = [[2, -1, 0], [-1, 2, -1], [0, -1, 4]]
        record, labels = audit.enumerate_orbits(gram)
        for index, matrix in zip(record["coordinate_root_indices_zero_based"], record["integral_reflection_matrices"]):
            for vector in itertools.product((0, 1), repeat=3):
                image = [sum(matrix[i][j] * vector[j] for j in range(3)) % 2 for i in range(3)]
                mask = sum(value << i for i, value in enumerate(vector))
                image_mask = sum(value << i for i, value in enumerate(image))
                if mask in labels:
                    self.assertEqual(labels[mask], labels[image_mask])

    def test_rootless_subgroup_is_not_mislabeled_full_automorphism_group(self):
        rootless = self.certificate["initial_lattices"]["rootless"]
        self.assertEqual(rootless["integral_reflection_matrices"], [])
        self.assertIsNone(rootless["explicit_orbits"])
        self.assertFalse(self.certificate["claim_boundary"]["full_integral_automorphism_group_computed"])
        self.assertEqual(rootless["orbit_size_histogram"], {"1": 65279})

    def test_hash_sampler_does_not_cover_all_root_reflection_orbits(self):
        records = self.certificate["initial_lattices"]
        self.assertEqual({name: item["initial_hash_window"]["distinct_reflection_subgroup_orbits_touched"] for name, item in records.items()}, {"origin": 23, "transparent": 10, "rootless": 32})


if __name__ == "__main__":
    unittest.main()
