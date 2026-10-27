import importlib.util
import math
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np

path = Path(__file__).resolve().parents[1] / "scripts/calculate_pocket_ef.py"
spec = importlib.util.spec_from_file_location("ef", path)
ef = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ef)


class EFTests(unittest.TestCase):
    def test_exact_counts_other_is_in_background(self):
        records = [dict(parent_id=str(i), source_entry=i, docking_score=float(i),
                        pocket="A" if i < 2 else "Other") for i in range(10)]
        result = ef.calculate_ef(records, ["A", "B"], [20])
        self.assertEqual(result[0]["EF"], 5)
        self.assertEqual((result[0]["N"], result[0]["K"], result[0]["n_p"], result[0]["k_p"]), (10, 2, 2, 2))
        self.assertIsNone(result[1]["EF"])

    def test_boundary_ties(self):
        records = [dict(parent_id=str(i), source_entry=i, docking_score=s, pocket="A")
                   for i, s in enumerate([-8, -7, -7, -6])]
        result = ef.calculate_ef(records, ["A"], [50])[0]
        self.assertEqual(result["K"], 3)
        self.assertEqual(result["EF"], 1)

    def test_duplicate_best_pose(self):
        best = {}
        for i, score in enumerate([-5, -7, -6]):
            ef.retain_best(best, dict(parent_id="x", source_entry=i, docking_score=score))
        self.assertEqual(best["x"]["source_entry"], 1)

    def test_outside_and_ambiguous(self):
        defs = [dict(name="A", center=[0, 0, 0], radius=8), dict(name="B", center=[4, 0, 0], radius=8)]
        self.assertEqual(ef.assign_region([2, 0, 0], defs)[0], "Ambiguous")
        self.assertEqual(ef.assign_region([0, 0, 0], defs)[0], "A")
        self.assertEqual(ef.assign_region([20, 0, 0], defs)[0], "Other")

    def test_unique_identity_and_percentages(self):
        row = dict(parent_id="a", source_entry=1, docking_score=-6, pocket="A")
        with self.assertRaises(ValueError):
            ef.calculate_ef([row, row], ["A"], [5])
        with self.assertRaises(ValueError):
            ef.calculate_ef([row], ["A"], [0])

    def test_anchor_missing_or_multiple(self):
        with self.assertRaises(ValueError):
            ef.default_definitions(SimpleNamespace(atom=[]))

    @staticmethod
    def receptor(chain="A"):
        atoms = []
        for group_index, (_, residues) in enumerate(ef.POCKET_ANCHORS):
            for index, (residue, number) in enumerate(residues):
                atoms.append(SimpleNamespace(resnum=number, pdbres=residue, pdbname=" CA ",
                                             chain=chain, x=group_index*20+index*2.,
                                             y=index*4., z=index*6.))
        return SimpleNamespace(atom=atoms)

    def test_correct_residue_sets_and_group_centroids(self):
        self.assertEqual(ef.POCKET_ANCHORS, [
            ("ARM2_3", [("MET", 202), ("THR", 205), ("LYS", 242)]),
            ("N_terminal_bottom", [("MET", 189), ("GLU", 226), ("LYS", 233)]),
            ("ARM4_5", [("LEU", 286), ("ASN", 326)]),
            ("C_terminal", [("LEU", 644), ("SER", 681)])])
        receptor = self.receptor()
        # Non-CA atoms do not affect the group center.
        receptor.atom.append(SimpleNamespace(resnum=202, pdbres="MET", pdbname="CB",
                                            chain="A", x=999., y=999., z=999.))
        definitions = ef.default_definitions(receptor)
        ef.validate_definitions(definitions)
        np.testing.assert_allclose(definitions[0]["center"], [2., 4., 6.])
        np.testing.assert_allclose(definitions[2]["center"], [41., 2., 3.])
        self.assertEqual([len(d["anchor_residues"]) for d in definitions], [3, 3, 2, 2])
        self.assertTrue(all(d["anchor_chain"] == "A" for d in definitions))

    def test_chain_selection_and_missing_anchor(self):
        receptor = self.receptor()
        receptor.atom += self.receptor("B").atom
        with self.assertRaises(ValueError):
            ef.default_definitions(receptor)
        definitions = ef.default_definitions(receptor, chain="B", radius=10.)
        self.assertTrue(all(d["anchor_chain"] == "B" and d["radius"] == 10 for d in definitions))
        receptor.atom = [a for a in receptor.atom if not (a.chain == "B" and a.resnum == 681)]
        with self.assertRaisesRegex(ValueError, "SER681"):
            ef.default_definitions(receptor, chain="B")

    def test_chain_mixing_is_rejected(self):
        receptor = self.receptor()
        receptor.atom[-1].chain = "B"
        with self.assertRaisesRegex(ValueError, "multiple chains"):
            ef.default_definitions(receptor)

    def test_pocket_definitions(self):
        defs = [dict(name=str(i), center=[i, 0, 0], radius=8) for i in range(4)]
        ef.validate_definitions(defs)
        defs[0]["radius"] = -1
        with self.assertRaises(ValueError):
            ef.validate_definitions(defs)


if __name__ == "__main__":
    unittest.main()
