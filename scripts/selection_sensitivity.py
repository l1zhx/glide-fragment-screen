#!/usr/bin/env python3
"""Reuse saved 100-group assignments; compare sequential gates without reclustering."""
import argparse
import csv
import json
from pathlib import Path
from sequential_selection import sequential_pool, select_candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--score-tops", default="1000,2000,3000")
    parser.add_argument("--le-tops", default="300,500,1000")
    parser.add_argument("--number", type=int, default=30)
    parser.add_argument("--max-per-cluster", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    with args.csv.open(encoding="utf-8-sig", newline="") as f:
        records = list(csv.DictReader(f))
    for r in records:
        for key in ("entry_number", "hac", "cluster_id", "cluster_size"):
            r[key] = int(r[key])
        for key in ("docking_score", "ligand_efficiency"):
            r[key] = float(r[key])
    summaries, candidate_sets = [], {}
    for score_top in map(int, args.score_tops.split(",")):
        for le_top in map(int, args.le_tops.split(",")):
            eligible, summary = sequential_pool(records, score_top, le_top)
            chosen = select_candidates(eligible, args.number, args.max_per_cluster)
            name = f"score{score_top}_le{le_top}"
            candidate_sets[name] = {r["entry_number"] for r in chosen}
            summary.update(setting=name, selected=len(chosen),
                           candidate_entries=[r["entry_number"] for r in chosen],
                           candidate_names=[r["title"] for r in chosen])
            summaries.append(summary)
    (args.out / "sensitivity.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    with (args.out / "overlap_jaccard.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["setting"] + list(candidate_sets))
        for name, a in candidate_sets.items():
            writer.writerow([name] + [len(a & b) / len(a | b) if a | b else ""
                                      for b in candidate_sets.values()])
    print(f"Compared {len(summaries)} settings; empty-set overlap is left blank.")


if __name__ == "__main__":
    main()
