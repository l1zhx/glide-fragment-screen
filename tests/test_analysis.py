import importlib
import json
import sys
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_top_hits import ranked_indices
from fixed_clusters import fixed_clusters, tanimoto_condensed
from sequential_selection import annotate_review_flags, sequential_pool, select_candidates
import export_manual_top30 as manual
import cluster_glide_top1w as cluster
import analyze_global_glide_1w as density


class AnalysisTests(unittest.TestCase):
    @staticmethod
    def record(entry, score, hac, cid, size):
        return dict(entry_number=entry, docking_score=score, hac=hac,
                    ligand_efficiency=-score/hac, cluster_id=cid, cluster_size=size,
                    review_flags="", title=f"T{entry}")

    def test_exact_100_nonempty_clusters(self):
        from rdkit import DataStructs
        rng = np.random.default_rng(73)
        fps = []
        for _ in range(125):
            fp = DataStructs.ExplicitBitVect(512)
            for bit in rng.choice(512, 32, replace=False):
                fp.SetBit(int(bit))
            fps.append(fp)
        records = [{} for _ in fps]
        groups, tree = fixed_clusters(records, fps, 100)
        self.assertEqual(len(groups), 100)
        self.assertEqual(sum(map(len, groups)), 125)
        self.assertTrue(all(groups))
        self.assertEqual(len({r["cluster_id"] for r in records}), 100)
        self.assertEqual(sum(r["is_cluster_medoid"] for r in records), 100)

    def test_scipy_condensed_order(self):
        from rdkit import Chem, DataStructs
        from scipy.spatial.distance import squareform
        records = [{"mol": Chem.MolFromSmiles(s)} for s in ["CCO", "CCCC", "c1ccccc1", "CCC"]]
        fps = cluster.make_fingerprints(records, 2, 2048)
        matrix = squareform(tanimoto_condensed(fps))
        for i in range(4):
            for j in range(i+1, 4):
                self.assertAlmostEqual(matrix[i, j], 1-DataStructs.TanimotoSimilarity(fps[i], fps[j]))

    def test_top_hits_sort_all_not_prefix(self):
        self.assertEqual(ranked_indices([(-2, 1), (-8, 2), (-8, 4), (-7, 3)], 2), [2, 4])

    def test_gates_cannot_be_rescued_by_cluster_size(self):
        records = [self.record(1, -8, 40, 1, 10), self.record(2, -7, 10, 2, 20),
                   self.record(3, -6, 20, 3, 10000)]
        eligible, report = sequential_pool(records, 2, 1)
        self.assertEqual([r["entry_number"] for r in eligible], [2])
        self.assertEqual(report["score_pool"], 2)
        self.assertEqual(report["le_pool"], 1)
        self.assertIsNone(report["weights"])

    def test_gate_ties_are_all_retained(self):
        records = [self.record(i, -6, 20, i, 1) for i in range(1, 5)]
        eligible, report = sequential_pool(records, 1, 1)
        self.assertEqual(len(eligible), 4)

    def test_cluster_size_order_and_strict_cap(self):
        records = [self.record(1, -9, 20, 1, 5), self.record(2, -8, 20, 2, 50),
                   self.record(3, -7, 20, 2, 50)]
        eligible, _ = sequential_pool(records, 3, 3)
        self.assertEqual([r["entry_number"] for r in select_candidates(eligible, 60)], [2, 1])
        self.assertEqual([r["entry_number"] for r in select_candidates(eligible, 60, 2)], [2, 1, 3])

    def test_alerts_annotate_unless_exclusion_requested(self):
        records = [self.record(1, -8, 10, 1, 1)]
        records[0]["pains"] = "test alert"
        annotate_review_flags(records)
        self.assertEqual(len(sequential_pool(records, 1, 1)[0]), 1)
        self.assertEqual(len(sequential_pool(records, 1, 1, True)[0]), 0)

    def test_invalid_cluster_number(self):
        with self.assertRaises(ValueError):
            fixed_clusters([{}], [None], 100)

    def test_manual_mapping_and_missing_rank(self):
        ranks = manual.parse_rows(manual.DEFAULT_ROWS)
        labels = manual.parse_pockets(manual.DEFAULT_POCKETS)
        self.assertEqual(len(ranks), 30)
        self.assertEqual(Counter(labels), {"Ser246": 16, "Asn326": 11, "Thr150": 2, "Lys233": 1})
        archived = json.loads((Path(__file__).resolve().parents[1] / "config/manual_selection.json").read_text())
        self.assertEqual(archived["ranks"], ranks)
        self.assertEqual(archived["pockets"], labels)
        with self.assertRaises(ValueError):
            manual.choose_records([], ranks)
        with self.assertRaises(ValueError):
            manual.assign_pockets([{}], labels)

    def test_manual_le_consistency_and_labels(self):
        row = dict(selection_rank=3, title="DEMO", docking_score=-6,
                   ligand_efficiency=.3, hac=20, cluster_id=4, cluster_size=2)
        selected = manual.choose_records([row], [3])
        manual.assign_pockets(selected, ["Asn326"])
        selected[0]["display_title"] = "DEMO"
        text = manual.make_legends(selected)[0]
        self.assertIn("Asn326", text)
        self.assertNotIn("Top60", text)
        self.assertNotIn("Selected", text)
        row["ligand_efficiency"] = .4
        with self.assertRaises(ValueError):
            manual.choose_records([row], [3])

    def test_identical_graphs_are_together_when_possible(self):
        from rdkit import Chem
        records = [{"mol": Chem.MolFromSmiles(s)} for s in ["Cc1ccccc1", "Cc1ccccc1", "CCO"]]
        fps = cluster.make_fingerprints(records, 2, 2048)
        groups, _ = cluster.cluster_records(records, fps, 2)
        self.assertEqual(sorted(map(len, groups)), [1, 2])
        self.assertEqual(records[0]["cluster_id"], records[1]["cluster_id"])

    def test_density_contains_ligands(self):
        records = [{"coords": np.array([[10., 0., 0.]])}]
        origin, counts = density.grid_definition(np.array([[0., 0., 0.]]), records, 1, 1)
        self.assertGreaterEqual((origin + counts - 1)[0], 11)

    def test_numpy_smoothing_keeps_small_grid_shape(self):
        with patch.dict(sys.modules, {"scipy": None, "scipy.ndimage": None}):
            result = density.gaussian_smooth(np.ones((2, 3, 4)), 2.5)
        self.assertEqual(result.shape, (2, 3, 4))
        self.assertTrue(np.all(np.isfinite(result)))


if __name__ == "__main__":
    unittest.main()
