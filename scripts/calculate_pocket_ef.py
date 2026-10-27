#!/usr/bin/env python3
"""Read Glide MAEGZ directly and calculate unadjusted regional EF.

Use licensed Schrodinger Python with RDKit and NumPy. Input files are unchanged.
Default geometry: ligand heavy-atom centroid within 8 A of the geometric mean
of a specified residue group's C-alpha coordinates, with a 2 A nearest/second-
nearest distance gap. These are EXPLORATORY residue-group neighborhoods,
not established anatomical pocket boundaries.
Use --definitions for reviewed centers/radii, or --label-property for existing
pose annotations. Missing property labels become Unassigned, never inferred.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import platform
from collections import Counter
from pathlib import Path
import numpy as np

SCRIPT_VERSION = "2.0-residue-groups-2026-10-22"
POCKET_ANCHORS = [
    ("ARM2_3", [("MET", 202), ("THR", 205), ("LYS", 242)]),
    ("N_terminal_bottom", [("MET", 189), ("GLU", 226), ("LYS", 233)]),
    ("ARM4_5", [("LEU", 286), ("ASN", 326)]),
    ("C_terminal", [("LEU", 644), ("SER", 681)]),
]


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def coordinates(st):
    return np.array([(a.x, a.y, a.z) for a in st.atom
                     if str(a.element).strip().upper() not in {"H", "D"}], dtype=float).reshape(-1, 3)


def default_definitions(st, chain=None, radius=8.):
    definitions = []
    ca_atoms = [a for a in st.atom if str(a.pdbname).strip().upper() == "CA"]
    for name, residues in POCKET_ANCHORS:
        anchors = []
        for residue, number in residues:
            candidates = [a for a in ca_atoms if int(a.resnum) == number
                          and str(a.pdbres).strip().upper() == residue
                          and (chain is None or str(a.chain).strip() == chain)]
            if len(candidates) != 1:
                raise ValueError(f"{name}: {residue}{number} expected one C-alpha; "
                                 f"found {len(candidates)}. Check construct/residue numbering, "
                                 "use --chain or reviewed --definitions. No missing anchor is omitted.")
            atom = candidates[0]
            anchors.append(dict(residue=residue, number=number, chain=str(atom.chain).strip(),
                                ca=[float(atom.x), float(atom.y), float(atom.z)]))
        chains = {a["chain"] for a in anchors}
        if len(chains) != 1:
            raise ValueError(f"{name}: anchor residues span multiple chains; select the correct --chain.")
        xyz = np.asarray([a["ca"] for a in anchors], dtype=float)
        if not np.all(np.isfinite(xyz)):
            raise ValueError(f"{name}: nonfinite anchor coordinates")
        center = xyz.mean(axis=0)
        definitions.append(dict(name=name, center=center.tolist(), radius=radius,
                                anchor_chain=anchors[0]["chain"], anchor_residues=anchors,
                                max_anchor_ca_distance=float(np.linalg.norm(xyz-center, axis=1).max()),
                                definition="EXPLORATORY sphere at residue-group C-alpha geometric centroid"))
    if len({d["anchor_chain"] for d in definitions}) != 1:
        raise ValueError("The four regions resolve to different protein chains; select one --chain.")
    return definitions


def validate_definitions(definitions):
    if len(definitions) != 4 or len({d["name"] for d in definitions}) != 4:
        raise ValueError("Define four distinctly named regions")
    reserved = {"Other", "Ambiguous", "Unassigned"}
    for d in definitions:
        center = np.asarray(d["center"], dtype=float)
        if (d["name"] in reserved or center.shape != (3,) or not np.all(np.isfinite(center))
                or not math.isfinite(float(d["radius"])) or float(d["radius"]) <= 0):
            raise ValueError("Each region needs a non-reserved name, finite xyz center and positive radius")


def assign_region(center, definitions, min_gap=2.):
    """Return Other outside all spheres; Ambiguous if multiple close candidates."""
    distances = {d["name"]: float(np.linalg.norm(np.asarray(center) - d["center"]))
                 for d in definitions}
    inside = sorted((distances[d["name"]], d["name"]) for d in definitions
                    if distances[d["name"]] <= float(d["radius"]))
    if not inside:
        return "Other", distances
    if len(inside) > 1 and inside[1][0] - inside[0][0] <= min_gap:
        return "Ambiguous", distances
    return inside[0][1], distances


def retain_best(best, record):
    key = record["parent_id"]
    previous = best.get(key)
    if previous is None or (record["docking_score"], record["source_entry"]) < (
            previous["docking_score"], previous["source_entry"]):
        best[key] = record


def calculate_ef(records, pocket_names, percentages):
    if not records:
        raise ValueError("No scored parent molecules")
    if not all(math.isfinite(r["docking_score"]) for r in records):
        raise ValueError("Nonfinite DockingScore")
    if len({r["parent_id"] for r in records}) != len(records):
        raise ValueError("EF requires one globally best pose per parent identity")
    ordered = sorted(records, key=lambda r: (r["docking_score"], r["source_entry"]))
    n = len(ordered)
    counts = Counter(r["pocket"] for r in ordered)
    results = []
    for percentage in percentages:
        if not math.isfinite(percentage) or not 0 < percentage <= 100:
            raise ValueError("Percentages must be in (0,100]")
        requested = max(1, math.ceil(n * percentage / 100))
        cutoff = ordered[requested - 1]["docking_score"]
        hits = [r for r in ordered if r["docking_score"] <= cutoff]
        k = len(hits)
        hit_counts = Counter(r["pocket"] for r in hits)
        for pocket in pocket_names:
            npocket, kpocket = counts[pocket], hit_counts[pocket]
            results.append(dict(pocket=pocket, top_percent_requested=percentage,
                                N=n, K_requested=requested, K=k, n_p=npocket, k_p=kpocket,
                                docking_score_cutoff=cutoff, background_fraction=npocket/n,
                                top_fraction=kpocket/k, expected_top_count=k*npocket/n,
                                EF=kpocket*n/(k*npocket) if npocket else None))
    return results


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def plot_ef(path, results, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    for name in dict.fromkeys(r["pocket"] for r in results):
        rows = sorted((r for r in results if r["pocket"] == name), key=lambda r: r["top_percent_requested"])
        ax.plot([r["top_percent_requested"] for r in rows],
                [np.nan if r["EF"] is None else r["EF"] for r in rows], marker="o", label=name)
    ax.axhline(1, linestyle="--", color="#777777", label="Background EF = 1")
    ax.set(xlabel="Requested global top fraction (%)", ylabel="EF", title=title)
    ax.legend()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def read_maegz(args):
    from schrodinger import structure
    from schrodinger.structutils.smiles import SmilesGenerator
    from rdkit import Chem, rdBase
    from rdkit.Chem.MolStandardize import rdMolStandardize
    smiles_generator = SmilesGenerator(unique=True, safe=True, wantAllH=False)
    uncharger = rdMolStandardize.Uncharger()
    enumerator = rdMolStandardize.TautomerEnumerator()
    best, issues, stats = {}, [], Counter()
    receptor = None
    reader = structure.StructureReader(str(args.input))
    with reader:
        for index, st in enumerate(reader, 1):
            stats["all_entries"] += 1
            if args.score_property not in st.property:
                stats["unscored_entries"] += 1
                residues = {str(a.pdbres).strip().upper() for a in st.atom}
                if receptor is None and st.atom_total > 100 and len(residues &
                        {"ALA", "ARG", "ASN", "GLU", "LEU", "LYS", "SER", "THR", "VAL"}) >= 5:
                    receptor = st.copy()
                    stats["detected_receptor_entry"] = index
                continue
            try:
                score = float(st.property[args.score_property])
                if not math.isfinite(score):
                    raise ValueError("Nonfinite DockingScore")
                xyz = coordinates(st)
                if len(xyz) == 0:
                    raise ValueError("No heavy atoms")
                if args.id_property:
                    parent = str(st.property[args.id_property]).strip()
                    if not parent:
                        raise ValueError("Empty requested parent ID property")
                else:
                    mol = Chem.MolFromSmiles(smiles_generator.getSmiles(st))
                    if mol is None:
                        raise ValueError("Could not parse generated SMILES")
                    mol = rdMolStandardize.FragmentParent(Chem.RemoveHs(mol))
                    mol = enumerator.Canonicalize(uncharger.uncharge(mol))
                    parent = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                    if not parent or not mol.GetNumHeavyAtoms():
                        raise ValueError("Empty standardized parent")
                raw_label = str(st.property.get(args.label_property, "") or "").strip() if args.label_property else ""
                retain_best(best, dict(parent_id=parent, source_entry=index,
                                      title=str(st.title or ""), docking_score=score,
                                      hac=len(xyz),
                                      center=xyz.mean(axis=0).tolist(), raw_label=raw_label))
                stats["usable_scored_poses"] += 1
            except Exception as exc:
                stats["rejected_scored_entries"] += 1
                issues.append(dict(source_entry=index, title=str(st.title or ""), error=str(exc)))
    if args.receptor:
        with structure.StructureReader(str(args.receptor)) as reader:
            receptor = next(iter(reader)).copy()
    stats["unique_parent_molecules"] = len(best)
    stats["python"] = platform.python_version()
    stats["rdkit"] = rdBase.rdkitVersion
    return list(best.values()), receptor, issues, dict(stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receptor", type=Path, help="Prepared receptor in ORIGINAL docking coordinates if absent from input")
    parser.add_argument("--definitions", type=Path, help="Reviewed JSON list of four {name,center:[x,y,z],radius} regions")
    parser.add_argument("--chain", help="Protein chain for all residue-group anchors; omit only if unambiguous")
    parser.add_argument("--radius", type=float, default=8.)
    parser.add_argument("--min-gap", type=float, default=2.)
    parser.add_argument("--top-percent", default="1,2,5,10")
    parser.add_argument("--score-property", default="r_i_docking_score")
    parser.add_argument("--id-property", help="Optional reliable original-parent ID; default uses standardized parent SMILES")
    parser.add_argument("--label-property", help="Use existing pose labels instead of geometry; unlabelled poses remain Unassigned")
    parser.add_argument("--pocket-names", default=",".join(name for name, _ in POCKET_ANCHORS))
    parser.add_argument("--cohort-label", default="input-file cohort; full-campaign coverage unverified")
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    percentages = list(dict.fromkeys(float(v) for v in args.top_percent.split(",")))
    names = [s.strip() for s in args.pocket_names.split(",")]
    if (not percentages or any(not math.isfinite(q) or not 0 < q <= 100 for q in percentages)
            or not math.isfinite(args.radius) or args.radius <= 0
            or not math.isfinite(args.min_gap) or args.min_gap < 0):
        parser.error("Invalid percentage, radius or distance gap")
    if not args.input.is_file():
        parser.error("Input file does not exist")
    if args.label_property and (args.definitions or args.receptor):
        parser.error("Use label-property or geometry options, not both")
    if args.label_property and (len(names) != 4 or len(set(names)) != 4 or any(not v for v in names)
                                or set(names) & {"Other", "Ambiguous", "Unassigned"}):
        parser.error("Label mode requires four distinct non-reserved pocket names")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    before = sha256(args.input)
    print("Reading complete input and retaining globally best pose per parent...", flush=True)
    records, receptor, issues, stats = read_maegz(args)
    write_csv(args.output_dir / "rejected_entries.csv", issues, ["source_entry", "title", "error"])
    if not records:
        raise ValueError("No usable scored parent molecules; check score field and rejected_entries.csv")
    if args.label_property:
        definitions = None
        for r in records:
            label = r["raw_label"]
            if label and label not in names + ["Other", "Ambiguous", "Unassigned"]:
                raise ValueError(f"Unknown pocket label {label!r}; check --pocket-names")
            r["pocket"] = label or "Unassigned"
        region_note = "Existing pose labels; completeness must be verified"
    else:
        if args.definitions:
            definitions = json.loads(args.definitions.read_text(encoding="utf-8-sig"))
            region_note = "User-supplied center/radius regions; verify geometric boundaries"
        else:
            if receptor is None:
                raise ValueError("No receptor found. Supply original prepared --receptor or reviewed --definitions.")
            definitions = default_definitions(receptor, args.chain, args.radius)
            region_note = "EXPLORATORY residue-group C-alpha centroid neighborhoods; verify pocket boundaries"
        validate_definitions(definitions)
        names = [d["name"] for d in definitions]
        for r in records:
            r["pocket"], distances = assign_region(r["center"], definitions, args.min_gap)
            r.update({"distance_to_"+name: value for name, value in distances.items()})
        (args.output_dir / "regions_used.json").write_text(json.dumps(definitions, indent=2), encoding="utf-8")
    results = calculate_ef(records, names, percentages)
    if before != sha256(args.input):
        raise RuntimeError("Input changed during analysis")
    for r in records:
        r.update(center_x=r["center"][0], center_y=r["center"][1], center_z=r["center"][2])
    fields = ["parent_id", "source_entry", "title", "docking_score", "hac", "pocket", "center_x", "center_y", "center_z"]
    if definitions:
        fields += ["distance_to_"+name for name in names]
    write_csv(args.output_dir / "molecule_pocket_assignments.csv", records, fields)
    write_csv(args.output_dir / "pocket_EF.csv", results, list(results[0]))
    audit = dict(script_version=SCRIPT_VERSION, stats=stats, input_sha256=before, cohort=args.cohort_label,
                 identity=args.id_property or "RDKit standardized canonical isomeric parent SMILES",
                 score_property=args.score_property, region_definition=region_note,
                 pocket_counts=dict(Counter(r["pocket"] for r in records)),
                 top_percent=percentages, include_boundary_ties=True, min_gap=args.min_gap,
                 region_centers_and_anchors=definitions,
                 background="ALL retained parents, including Other/Ambiguous/Unassigned",
                 note="Winning-pose enrichment, not per-pocket binding hit rate or calibrated druggability")
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(region_note, flush=True)
    print(json.dumps(audit["pocket_counts"], indent=2), flush=True)
    if not args.skip_plot:
        plot_ef(args.output_dir / "pocket_EF.png", results,
                "Residue-group neighborhood EF (exploratory)" if not args.definitions and not args.label_property else "Regional docking enrichment")
    print("Finished: " + str(args.output_dir), flush=True)


if __name__ == "__main__":
    main()
