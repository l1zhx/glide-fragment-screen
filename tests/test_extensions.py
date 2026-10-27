import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import review_cluster_pockets as review
import plot_volume_occupancy as volume
import plot_le_analysis as le


class ExtensionTests(unittest.TestCase):
    def test_review_hash_and_missing_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            src = root/'source.csv'
            rows = [dict(entry_number=i, cluster_id=1, parent_smiles=f'C{i}') for i in [1, 2]]
            review.write_csv(src, rows, list(rows[0]))
            sheet = review.make_sheet(src)
            sheet[0]['pocket_label'] = 'ARM4_5'
            target = root/'sheet.csv'
            review.write_csv(target, sheet, list(sheet[0]))
            summary = review.summarize(src, target)[0]
            self.assertEqual(summary['review_coverage'], .5)
            self.assertEqual(summary['majority_fraction_reviewed'], 1.)
            src.write_text(src.read_text(encoding='utf-8-sig')+'\n', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'DIFFERENT'):
                review.summarize(src, target)

    def test_review_rejects_duplicate_entry(self):
        with tempfile.TemporaryDirectory() as temp:
            src = Path(temp)/'s.csv'
            review.write_csv(src, [{'entry_number':1}, {'entry_number':1}], ['entry_number'])
            with self.assertRaises(ValueError):
                review.make_sheet(src)

    def test_time_join_rejects_bad_times_and_origin(self):
        for states in [[1, 1, 2], [1, float('nan')], [2, 3]]:
            with self.assertRaises(ValueError):
                volume.match_times([1, 2], states)
        np.testing.assert_array_equal(volume.match_times([1, 11], [1, 2, 11]), [0, 2])
        with self.assertRaises(ValueError):
            volume.parse_bool('')

    def test_annotation_preserves_zeros_and_unknown_other_states(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            volume.write_csv(root/'v.csv', [dict(time_ns=i, volume_A3=v, nb_AS=n)
                for i, v, n in [(1, 100, 5), (11, 0, 0), (21, 200, 10), (31, 300, 12)]])
            volume.write_csv(root/'o.csv', [{'time_ns':i, 'MED_ARM4-5':b} for i, b in [(1,1),(11,1),(21,0),(31,0)]])
            args = SimpleNamespace(volume=root/'v.csv', occupancy=root/'o.csv', out=root/'out',
                state_column='MED_ARM4-5', ligand='MED', time_tolerance_ns=.001, descriptors=None)
            with patch.object(volume, 'plots'):
                volume.annotate(args)
            _, joined = volume.read_csv(root/'out/volume_occupancy_joined.csv')
            self.assertEqual(float(joined[1]['volume_A3']), 0.)
            self.assertEqual(joined[0]['other_fragments_contact_defined'], 'unknown')
            _, stats = volume.read_csv(root/'out/conditional_volume_summary.csv')
            self.assertEqual(float(stats[0]['detection_fraction']), .5)
            self.assertEqual(float(stats[0]['median_detected_positive_A3']), 100.)

    def test_plot_omits_reference_and_keeps_zero_marker(self):
        from matplotlib.figure import Figure
        t=np.array([1.,11.,21.,31.]); v=np.array([200.,0.,250.,300.]); bound=np.array([1,1,0,0],bool)
        detected=v>0
        stats=[volume.conditional_summary(v,detected,s,label) for label,s in [('present',bound),('absent',~bound)]]
        with patch.object(Figure, 'savefig', autospec=True) as save:
            volume.plots(Path('unused'),t,v,detected,t,bound,None,bound,stats,'MED',reference_volume=0.)
        fig=save.call_args_list[0].args[0]
        labels=[line.get_label() for ax in fig.axes for line in ax.lines]
        self.assertFalse(any('Static reference' in label for label in labels))
        self.assertTrue(any('Detected cavity' in label for label in labels))

    def test_le_recomputed_from_best_parent_score(self):
        with tempfile.TemporaryDirectory() as temp:
            src=Path(temp)/'s.csv'
            review.write_csv(src, [dict(parent_id='x',docking_score=-5,hac=10),
                                  dict(parent_id='x',docking_score=-7,hac=10)],
                             ['parent_id','docking_score','hac'])
            records,audit=le.load_records(src)
            self.assertEqual(records[0]['LE_score'], .7)
            self.assertEqual(audit['removed_duplicates'], 1)


if __name__ == '__main__':
    unittest.main()
