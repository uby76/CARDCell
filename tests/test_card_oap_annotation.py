#!/usr/bin/env python3
import csv
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location("calculate_results", ROOT / "lib/calculate_results.py")
CALC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CALC)


class CardOapAnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db, cls.size = CALC.card_oap_index(ROOT / "resources/CARD_OAP_full_6059_annotation.csv")

    def test_database_size_and_unique_aro_keys(self):
        self.assertEqual(self.size, 6059)
        self.assertEqual(len(self.db), 6059)

    def test_known_class_assignments(self):
        expected = {
            "ARO:3000498": "MLS",
            "ARO:3000412": "Sulfonamide",
            "ARO:3000501": "Rifamycin",
            "ARO:3000378": "Multidrug",
        }
        for aro, arg_class in expected.items():
            self.assertEqual(self.db[aro]["Class"], arg_class)

    def test_example_results_use_database_classes(self):
        result = ROOT / "examples/SRR6468562/expected_results/06_ARG_copies_per_cell.tsv"
        with result.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertGreater(len(rows), 0)
        for row in rows:
            aro = CALC.normalize_aro(row["ARO_Accession"])
            self.assertIn(aro, self.db)
            self.assertEqual(row["Drug_Class"], self.db[aro]["Class"])

    def test_example_class_total_is_conserved(self):
        result = ROOT / "examples/SRR6468562/expected_results/06_ARG_copies_per_cell.tsv"
        classes = ROOT / "examples/SRR6468562/expected_results/08_Drug_Class_abundance.tsv"
        with result.open(newline="") as handle:
            per_reference = sum(float(r["ARG_copies_per_cell"]) for r in csv.DictReader(handle, delimiter="\t"))
        with classes.open(newline="") as handle:
            per_class = sum(float(r["ARG_copies_per_cell"]) for r in csv.DictReader(handle, delimiter="\t"))
        self.assertAlmostEqual(per_reference, per_class, places=11)


if __name__ == "__main__":
    unittest.main()
