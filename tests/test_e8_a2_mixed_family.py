from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "derive_e8_a2_mixed_family.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_mixed_family.json"
SEMISTABLE_MODULE_PATH = ROOT / "research" / "derive_e8_a2_semistable_family.py"
SEMISTABLE_CERTIFICATE_PATH = (
    ROOT / "certificates" / "e8_a2_semistable_family.json"
)
SMALL_PRIME_CERTIFICATE_PATH = (
    ROOT / "certificates" / "e8_a2_mixed_target_small_primes.json"
)
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location("derive_e8_a2_mixed_family", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

SEMISTABLE_SPEC = importlib.util.spec_from_file_location(
    "derive_e8_a2_semistable_family", SEMISTABLE_MODULE_PATH
)
assert SEMISTABLE_SPEC is not None and SEMISTABLE_SPEC.loader is not None
SEMISTABLE_MODULE = importlib.util.module_from_spec(SEMISTABLE_SPEC)
sys.modules[SEMISTABLE_SPEC.name] = SEMISTABLE_MODULE
SEMISTABLE_SPEC.loader.exec_module(SEMISTABLE_MODULE)


class E8A2MixedFamilyTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact polynomial algebra unavailable")
    def test_certificate_recomputes_exactly(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(
            computed["family"]["kodaira_data"]["t=infinity"]["type"], "II*"
        )
        self.assertEqual(
            computed["target_gram"]["neron_severi_absolute_discriminant"], 948
        )
        self.assertFalse(computed["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact polynomial algebra unavailable")
    def test_modular_examples_cover_only_the_three_diagonal_profiles(self) -> None:
        examples = MODULE.compute_certificate()["modular_engine_checks"]["examples"]
        self.assertEqual(
            {example["shioda_height"] for example in examples},
            {"8/3", "10/3", "4"},
        )
        self.assertEqual(len({(item["c"], item["lambda"]) for item in examples}), 3)
        for example in examples:
            self.assertEqual(example["identity_residual"], [0] * 13)

    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact polynomial algebra unavailable")
    def test_full_semistable_chart_recomputes_and_contains_mixed_chart(self) -> None:
        computed = SEMISTABLE_MODULE.compute_certificate()
        committed = json.loads(
            SEMISTABLE_CERTIFICATE_PATH.read_text(encoding="utf-8")
        )
        self.assertEqual(computed, committed)
        self.assertEqual(computed["parameter_count"]["moduli_dimension"], 4)
        self.assertEqual(computed["mixed_subfamily"]["dimension"], 2)
        self.assertEqual(
            computed["kodaira_open_conditions"]["configuration"],
            "II* + I3 + I3 + I3 + 5 I1",
        )
        self.assertFalse(computed["claim_boundary"]["rank31_curve_found"])

    def test_small_prime_target_search_totals_are_internally_consistent(self) -> None:
        payload = json.loads(
            SMALL_PRIME_CERTIFICATE_PATH.read_text(encoding="utf-8")
        )
        runs = payload["runs"]
        self.assertEqual(
            payload["totals"]["section_tests"],
            sum(sum(run["tests"].values()) for run in runs),
        )
        self.assertEqual(
            payload["totals"]["retained_sections"],
            sum(sum(run["sections"].values()) for run in runs),
        )
        self.assertEqual(
            [run["complete_gram_triples"] for run in runs], [0, 0, 0]
        )


if __name__ == "__main__":
    unittest.main()
