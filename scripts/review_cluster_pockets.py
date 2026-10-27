#!/usr/bin/env python3
"""Create a Maestro review sheet; summarize reviewed chemical-cluster/pocket labels.

Labels are supplied by the researcher, never inferred from cluster identity.
The source CSV hash and identities must still match the original review sheet.
"""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from prepare_top_hits import file_sha256

POCKETS = ['ARM2_3', 'N_terminal_bottom', 'ARM4_5', 'C_terminal', 'Other', 'Ambiguous']


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with Path(path).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def make_sheet(source):
    digest = file_sha256(source)
    rows = read_csv(source)
    seen = set()
    for row in rows:
        key = int(row['entry_number'])
        if key in seen:
            raise ValueError('Duplicate source entry_number')
        seen.add(key)
        row.update(source_csv_sha256=digest, pocket_label='', keep='', review_notes='')
    return rows


def summarize(source, sheet):
    original = {int(r['entry_number']): r for r in read_csv(source)}
    digest = file_sha256(source)
    labels, seen = {}, set()
    for row in read_csv(sheet):
        key = int(row['entry_number'])
        if key in seen or key not in original:
            raise ValueError('Duplicate or unknown review entry_number')
        seen.add(key)
        if row['source_csv_sha256'] != digest:
            raise ValueError('Review sheet belongs to a DIFFERENT source CSV/run')
        ref = original[key]
        if row['parent_smiles'] != ref['parent_smiles'] or int(row['cluster_id']) != int(ref['cluster_id']):
            raise ValueError('Review identity or cluster was edited')
        label = row.get('pocket_label', '').strip()
        if label and label not in POCKETS:
            raise ValueError(f'Unknown label {label!r}; allowed: {POCKETS}')
        labels[key] = label
    grouped = defaultdict(list)
    for key, ref in original.items():
        grouped[int(ref['cluster_id'])].append(labels.get(key, ''))
    result = []
    for cluster_id, group in sorted(grouped.items()):
        counts = Counter(x for x in group if x)
        n = sum(counts.values())
        # Other/Ambiguous remain explicit categories in the descriptive summary.
        entropy = -sum((v/n)*math.log2(v/n) for v in counts.values()) if n else None
        row = dict(cluster_id=cluster_id, n_in_source=len(group), n_reviewed=n,
                   review_coverage=n/len(group),
                   majority_label=counts.most_common(1)[0][0] if n else '',
                   majority_fraction_reviewed=max(counts.values())/n if n else None,
                   label_entropy_bits=entropy)
        row.update({p: counts[p] for p in POCKETS})
        result.append(row)
    return result


def plot(rows, path):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matrix = np.array([[r[p]/r['n_reviewed'] if r['n_reviewed'] else np.nan for p in POCKETS] for r in rows])
    fig, ax = plt.subplots(figsize=(9, max(4, len(rows)*.17)), constrained_layout=True)
    im = ax.imshow(matrix, aspect='auto', vmin=0, vmax=1, cmap='Blues')
    ax.set(xticks=range(len(POCKETS)), xticklabels=POCKETS,
           yticks=range(len(rows)), yticklabels=[f"C{r['cluster_id']} ({r['n_reviewed']}/{r['n_in_source']})" for r in rows],
           title='Reviewed pocket-label fractions by chemical cluster', ylabel='Cluster (reviewed / in source)')
    ax.tick_params(axis='x', rotation=30)
    fig.colorbar(im, ax=ax, label='Fraction among reviewed members')
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    init = sub.add_parser('init')
    init.add_argument('--csv', type=Path, required=True)
    init.add_argument('--out', type=Path, required=True)
    summary = sub.add_parser('summarize')
    summary.add_argument('--csv', type=Path, required=True)
    summary.add_argument('--review', type=Path, required=True)
    summary.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'init':
        rows = make_sheet(args.csv)
        if not rows:
            raise ValueError('No source molecules')
        if args.out.exists():
            raise FileExistsError(args.out)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        write_csv(args.out, rows, list(rows[0]))
    else:
        rows = summarize(args.csv, args.review)
        if not rows:
            raise ValueError('No source molecules')
        args.out.mkdir(parents=True, exist_ok=False)
        write_csv(args.out/'cluster_pocket_summary.csv', rows, list(rows[0]))
        plot(rows, args.out/'cluster_pocket_heatmap.png')
        (args.out/'review_audit.json').write_text(json.dumps(dict(
            source_sha256=file_sha256(args.csv), review_sha256=file_sha256(args.review),
            source_population='Exactly the supplied CSV; Top60 is NOT a whole-cluster census',
            interpretation='Descriptive annotation only; incomplete review can bias label fractions'), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
