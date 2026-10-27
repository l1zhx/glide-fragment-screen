#!/usr/bin/env python3
"""Audit a complete pose file and export score-sorted Top N without changing it.

Run with licensed Schrodinger Python. Sorting is by scored POSE, not unique parent.
"""
import argparse
import csv
import hashlib
import heapq
import json
import math
from collections import Counter
from pathlib import Path


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ranked_indices(items, number):
    """Stable score ordering, source entry breaks exact ties."""
    if number < 1:
        raise ValueError("number must be positive")
    return [index for score, index in heapq.nsmallest(number, items)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--number", type=int, default=10000)
    parser.add_argument("--score-property", default="r_i_docking_score")
    parser.add_argument("--receptor-entry", type=int,
                        help="Explicit 1-based unscored receptor entry; otherwise detect protein residues")
    args = parser.parse_args()
    if args.number < 1 or not args.input.is_file():
        parser.error("Use a valid input file and positive --number")
    from schrodinger import structure
    args.out.mkdir(parents=True, exist_ok=False)
    before = file_sha256(args.input)
    items, names, properties = [], {}, Counter()
    receptor = None
    receptor_index = None
    invalid = unscored = total = 0
    for total, st in enumerate(structure.StructureReader(str(args.input)), 1):
        properties.update(st.property.keys())
        if args.score_property not in st.property:
            unscored += 1
            # Protein heuristic is explicitly reported; unusual constructs can override it.
            residues = {str(r.pdbres).strip().upper() for r in st.residue}
            amino = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                     "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                     "THR", "TRP", "TYR", "VAL"}
            chosen = total == args.receptor_entry if args.receptor_entry else (
                receptor is None and st.atom_total > 100 and len(residues & amino) >= 5)
            if chosen:
                receptor, receptor_index = st.copy(), total
            continue
        try:
            score = float(st.property[args.score_property])
        except (ValueError, TypeError):
            invalid += 1
            continue
        if not math.isfinite(score):
            invalid += 1
            continue
        items.append((score, total))
        names[total] = str(st.title or "")
    if not items:
        raise RuntimeError(f"No finite {args.score_property} entries; detected: {sorted(properties)}")
    if args.receptor_entry is not None and receptor is None:
        raise ValueError("Explicit receptor entry not found or it has a score property")
    indices = ranked_indices(items, args.number)
    wanted = set(indices)
    retained = {}
    for index, st in enumerate(structure.StructureReader(str(args.input)), 1):
        if index in wanted:
            retained[index] = st.copy()
    if before != file_sha256(args.input):
        raise RuntimeError("Input changed while reading")
    if set(retained) != wanted:
        raise RuntimeError("Missing entries on second pass")
    for filename, order in [("top1w_SP_pv.maegz", indices), ("top100_SP_pv.maegz", indices[:100])]:
        with structure.StructureWriter(str(args.out / filename)) as writer:
            if receptor is not None:
                writer.append(receptor)
            for index in order:
                st = retained[index]
                st.property["i_user_original_entry"] = index
                writer.append(st)
    score_by_index = {index: score for score, index in items}
    with (args.out / "pose_index.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["score_rank", "original_entry", "title", "docking_score"])
        for rank, index in enumerate(indices, 1):
            writer.writerow([rank, index, names[index], score_by_index[index]])
    report = {"input_sha256": before, "score_property": args.score_property,
              "total_entries": total, "finite_scored_poses": len(items),
              "unscored_entries": unscored, "invalid_scores": invalid,
              "input_sorted_by_score": all(a[0] <= b[0] for a, b in zip(items, items[1:])),
              "selected_poses": len(indices), "requested_poses": args.number,
              "receptor_entry": receptor_index, "deduplication": "none: pose-level",
              "scope": "Top N within the supplied file, not proof of full-library coverage"}
    (args.out / "input_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
