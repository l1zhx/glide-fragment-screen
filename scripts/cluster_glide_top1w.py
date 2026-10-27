#!/usr/bin/env python3
"""Cluster the top Glide hits by 2D chemical similarity.

Designed for Schrödinger 2024-3 and later. Run with Schrödinger Python:

    $SCHRODINGER/run python3 cluster_glide_top1w.py

The script:
  1. reads scored ligand entries from a Maestro/pose-viewer file;
  2. keeps the best-scoring pose for each standardized parent molecule;
  3. cuts average-linkage Tanimoto clustering into exactly 100 groups;
  4. writes CSV tables, plots, SVG structure grids, and annotated MAEGZ files;
  5. applies sequential score/LE gates, then cluster-size review order.

Plot axes omit units. LE_score is -DockingScore/HAC, not experimental LE.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import json
from pathlib import Path
from fixed_clusters import fixed_clusters, make_hierarchy_plots, plot_selection_stages
from sequential_selection import annotate_review_flags, sequential_pool, select_candidates

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

try:
    from rdkit.Chem import rdFingerprintGenerator
except ImportError:
    rdFingerprintGenerator = None

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    rdMolStandardize = None

try:
    from schrodinger import structure
except ImportError:
    structure = None  # Pure RDKit demo/tests do not read Maestro files.

try:
    from schrodinger.rdkit import rdkit_adapter
except ImportError:
    rdkit_adapter = None


DEFAULT_INPUT = (
    "data/private/truncated_hits..maegz"
)
DEFAULT_OUTPUT = "results/cluster_top1w"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fixed-count Tanimoto clustering and sequential candidate selection."
    )
    parser.add_argument("-i", "--input", default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output-dir", default=DEFAULT_OUTPUT)
    parser.add_argument("--input-csv", type=Path,
                        help="Recluster saved parent_smiles/docking_score CSV without Schrodinger; no MAEGZ export")
    parser.add_argument(
        "--max-molecules",
        type=int,
        default=10000,
        help="Maximum number of scored ligand entries to read (default: 10000).",
    )
    parser.add_argument("--n-clusters", type=int, default=100)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--fp-size", type=int, default=2048)
    parser.add_argument("--score-top", type=int, default=2000,
                        help="DockingScore gate: retain best N, including boundary ties")
    parser.add_argument("--le-top", type=int, default=500,
                        help="Within score gate: retain best N LE values including ties")
    parser.add_argument("--exclude-review-flags", action="store_true",
                        help="Explicitly exclude alerts after both gates; default annotates only")
    parser.add_argument(
        "--top30", type=int, default=30, help="Size of automatic diverse set."
    )
    parser.add_argument(
        "--review-set",
        type=int,
        default=60,
        help="Size of larger diversity-first manual review set.",
    )
    parser.add_argument(
        "--max-per-cluster",
        type=int,
        default=1,
        help="Maximum molecules per cluster in automatic selections.",
    )
    parser.add_argument(
        "--skip-pains",
        action="store_true",
        help="Skip PAINS annotation to save a little time.",
    )
    return parser.parse_args()


def log(message):
    print(f"[cluster_glide_top1w] {message}", flush=True)


def score_from_structure(st):
    preferred = (
        "r_i_docking_score",
        "r_i_glide_gscore",
        "r_glide_gscore",
        "r_glide_docking_score",
    )
    for key in preferred:
        try:
            value = float(st.property[key])
            if math.isfinite(value):
                return value, key
        except (KeyError, TypeError, ValueError):
            pass

    for key in st.property.keys():
        low = key.lower()
        if "docking_score" in low or low.endswith("glide_gscore"):
            try:
                value = float(st.property[key])
                if math.isfinite(value):
                    return value, key
            except (TypeError, ValueError):
                pass
    return None, None


def structure_to_rdkit(st):
    if rdkit_adapter is not None:
        mol = rdkit_adapter.to_rdkit(
            st,
            implicitH=True,
            include_properties=False,
            sanitize=True,
            include_stereo=True,
        )
    else:
        from schrodinger import adapter

        mol = adapter.to_rdkit(st)
        mol = Chem.RemoveHs(mol, sanitize=True)
    if mol is None:
        raise ValueError("Schrodinger-to-RDKit conversion returned None")
    return mol


def standardized_parent(mol):
    mol = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
    if rdMolStandardize is None:
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fragments:
            mol = max(fragments, key=lambda item: item.GetNumHeavyAtoms())
        return mol

    try:
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    except Exception:
        # Standardization is for clustering only; retain a usable molecule if an
        # unusual valence or metal prevents one of the optional standardizers.
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fragments:
            mol = max(fragments, key=lambda item: item.GetNumHeavyAtoms())
    Chem.SanitizeMol(mol)
    return mol


def make_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def pains_annotation(mol, catalog):
    if catalog is None:
        return ""
    match = catalog.GetFirstMatch(mol)
    return "" if match is None else match.GetDescription()


def molecular_properties(mol):
    return {
        "hac": int(mol.GetNumHeavyAtoms()),
        "mw": float(Descriptors.MolWt(mol)),
        "clogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotb": int(Lipinski.NumRotatableBonds(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
    }


def read_and_deduplicate(input_path, max_molecules, pains_catalog):
    by_parent = {}
    scored_seen = 0
    conversion_failures = 0
    score_property_names = set()

    with structure.StructureReader(str(input_path)) as reader:
        for entry_number, st in enumerate(reader, start=1):
            score, score_key = score_from_structure(st)
            if score is None:
                continue
            scored_seen += 1
            score_property_names.add(score_key)
            if scored_seen > max_molecules:
                break

            try:
                original_mol = structure_to_rdkit(st)
                parent = standardized_parent(original_mol)
                parent_smiles = Chem.MolToSmiles(
                    parent, canonical=True, isomericSmiles=True
                )
                if not parent_smiles:
                    raise ValueError("empty canonical SMILES")
            except Exception as exc:
                conversion_failures += 1
                if conversion_failures <= 10:
                    log(f"Skipping entry {entry_number}: {exc}")
                continue

            props = molecular_properties(parent)
            title = str(getattr(st, "title", "") or "").strip()
            if not title:
                title = f"entry_{entry_number}"
            ligand_efficiency = -score / props["hac"] if props["hac"] else math.nan

            record = {
                "entry_number": entry_number,
                "title": title,
                "docking_score": float(score),
                "ligand_efficiency": float(ligand_efficiency),
                "parent_smiles": parent_smiles,
                "mol": parent,
                "pains": pains_annotation(parent, pains_catalog),
                "state_count": 1,
                **props,
            }

            previous = by_parent.get(parent_smiles)
            if previous is None:
                by_parent[parent_smiles] = record
            else:
                state_count = previous["state_count"] + 1
                if score < previous["docking_score"]:
                    record["state_count"] = state_count
                    by_parent[parent_smiles] = record
                else:
                    previous["state_count"] = state_count

    records = list(by_parent.values())
    records.sort(key=lambda item: item["docking_score"])
    for rank, record in enumerate(records, start=1):
        record["score_rank"] = rank

    log(f"Scored entries read: {min(scored_seen, max_molecules)}")
    log(f"Unique standardized parent molecules: {len(records)}")
    log(f"Conversion failures: {conversion_failures}")
    log(f"DockingScore property detected: {sorted(score_property_names)}")
    return records


def make_fingerprints(records, radius, fp_size):
    if rdFingerprintGenerator is not None:
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=fp_size
        )
        return [generator.GetFingerprint(record["mol"]) for record in records]

    from rdkit.Chem import AllChem

    return [
        AllChem.GetMorganFingerprintAsBitVect(
            record["mol"], radius, nBits=fp_size
        )
        for record in records
    ]


def read_parent_csv(path, pains_catalog=None):
    """Reuse parent chemistry and scores; recompute descriptors/LE and deduplicate."""
    by_parent = {}
    seen_entries = set()
    catalog = pains_catalog
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            entry = int(row["entry_number"])
            if entry in seen_entries:
                raise ValueError(f"Duplicate CSV entry_number {entry}")
            seen_entries.add(entry)
            mol = Chem.MolFromSmiles(row["parent_smiles"])
            if mol is None:
                raise ValueError(f"Invalid parent_smiles at entry {entry}")
            mol = standardized_parent(mol)
            smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            score = float(row["docking_score"])
            if not math.isfinite(score) or not mol.GetNumHeavyAtoms():
                raise ValueError(f"Invalid score or HAC at entry {entry}")
            props = molecular_properties(mol)
            record = dict(entry_number=entry, title=row.get("title", "") or f"entry_{entry}",
                          docking_score=score, ligand_efficiency=-score / props["hac"],
                          parent_smiles=smiles, mol=mol, pains=pains_annotation(mol, catalog),
                          state_count=int(row.get("state_count") or 1), **props)
            previous = by_parent.get(smiles)
            if previous is None or (score, entry) < (previous["docking_score"], previous["entry_number"]):
                by_parent[smiles] = record
    records = sorted(by_parent.values(), key=lambda r: (r["docking_score"], r["entry_number"]))
    for rank, r in enumerate(records, 1):
        r["score_rank"] = rank
    return records


def cluster_records(records, fingerprints, n_clusters=100):
    return fixed_clusters(records, fingerprints, n_clusters)


CSV_FIELDS = [
    "score_rank",
    "entry_number",
    "title",
    "docking_score",
    "hac",
    "ligand_efficiency",
    "cluster_id",
    "cluster_size",
    "leaf_index",
    "is_cluster_medoid",
    "similarity_to_medoid",
    "cluster_mean_pair_similarity",
    "pass_score_gate",
    "pass_le_gate",
    "pass_review_gate",
    "le_rank_within_score_pool",
    "eligible_cluster_size",
    "state_count",
    "mw",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotb",
    "formal_charge",
    "pains",
    "review_flags",
    "parent_smiles",
]


def write_records_csv(path, records, include_selection_rank=False):
    fields = (["selection_rank"] if include_selection_rank else []) + CSV_FIELDS
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            row = record.copy()
            for key in (
                "docking_score",
                "ligand_efficiency",
                "similarity_to_medoid",
                "cluster_mean_pair_similarity",
                "mw",
                "clogp",
                "tpsa",
            ):
                row[key] = format(float(row[key]), ".17g")
            writer.writerow(row)


def cluster_summary(records):
    grouped = {}
    for record in records:
        grouped.setdefault(record["cluster_id"], []).append(record)

    summaries = []
    for cluster_id in sorted(grouped):
        members = grouped[cluster_id]
        scores = np.asarray([item["docking_score"] for item in members])
        les = np.asarray([item["ligand_efficiency"] for item in members])
        centroid = next(item for item in members if item["is_cluster_medoid"])
        best = min(members, key=lambda item: item["docking_score"])
        summaries.append(
            {
                "cluster_id": cluster_id,
                "cluster_size": len(members),
                "best_title": best["title"],
                "best_docking_score": float(scores.min()),
                "median_docking_score": float(np.median(scores)),
                "mean_docking_score": float(scores.mean()),
                "best_ligand_efficiency": float(les.max()),
                "median_ligand_efficiency": float(np.median(les)),
                "medoid_title": centroid["title"],
                "mean_pair_tanimoto": float(members[0]["cluster_mean_pair_similarity"]),
                "eligible_members": sum(item.get("pass_review_gate", 0) for item in members),
                "medoid_smiles": centroid["parent_smiles"],
            }
        )
    return summaries


def write_summary_csv(path, summaries):
    fields = [
        "cluster_id",
        "cluster_size",
        "best_title",
        "best_docking_score",
        "median_docking_score",
        "mean_docking_score",
        "best_ligand_efficiency",
        "median_ligand_efficiency",
        "medoid_title",
        "medoid_smiles",
        "mean_pair_tanimoto",
        "eligible_members",
    ]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            row = summary.copy()
            for key in fields:
                if isinstance(row[key], float):
                    row[key] = f"{row[key]:.4f}"
            writer.writerow(row)


def make_overview_plot(path, records, summaries):
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)

    top_summaries = summaries[: min(30, len(summaries))]
    cluster_labels = [str(item["cluster_id"]) for item in top_summaries]
    axes[0, 0].bar(
        cluster_labels,
        [item["cluster_size"] for item in top_summaries],
        color="#4C78A8",
    )
    axes[0, 0].set_title("Largest clusters")
    axes[0, 0].set_xlabel("Cluster ID")
    axes[0, 0].set_ylabel("Number of molecules")
    axes[0, 0].tick_params(axis="x", labelrotation=90)

    sizes = [item["cluster_size"] for item in summaries]
    max_size = max(sizes) if sizes else 1
    bins = np.arange(0.5, min(max_size, 50) + 1.5, 1.0)
    clipped_sizes = np.minimum(sizes, 50)
    axes[0, 1].hist(clipped_sizes, bins=bins, color="#72B7B2", edgecolor="none")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("Cluster-size distribution (50 includes 50+)")
    axes[0, 1].set_xlabel("Cluster size")
    axes[0, 1].set_ylabel("Number of clusters (log scale)")

    top20_ids = [item["cluster_id"] for item in summaries[:20]]
    grouped_scores = [
        [
            record["docking_score"]
            for record in records
            if record["cluster_id"] == cluster_id
        ]
        for cluster_id in top20_ids
    ]
    axes[1, 0].boxplot(grouped_scores, showfliers=False)
    axes[1, 0].set_title("DockingScore in the 20 largest clusters")
    axes[1, 0].set_xlabel("Cluster ID")
    axes[1, 0].set_ylabel("DockingScore")
    axes[1, 0].set_xticks(range(1, len(top20_ids) + 1), top20_ids, rotation=90)

    top10_ids = {item["cluster_id"] for item in summaries[:10]}
    colors = [
        record["cluster_id"] if record["cluster_id"] in top10_ids else 0
        for record in records
    ]
    scatter = axes[1, 1].scatter(
        [record["docking_score"] for record in records],
        [record["ligand_efficiency"] for record in records],
        c=colors,
        s=12,
        alpha=0.55,
        cmap="tab10",
        linewidths=0,
    )
    axes[1, 1].set_title("DockingScore vs ligand efficiency")
    axes[1, 1].set_xlabel("DockingScore")
    axes[1, 1].set_ylabel("Ligand efficiency")

    fig.suptitle("Glide top-hit chemical-similarity clustering", fontsize=16)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def ascii_text(value):
    text = str(value).encode("ascii", "replace").decode("ascii")
    return re.sub(r"\s+", " ", text).strip()


def write_structure_grid_svg(path, records, legends, mols_per_row=5):
    if not records:
        return
    mols = []
    for record in records:
        mol = Chem.Mol(record["mol"])
        try:
            rdMolDescriptors.CalcNumRings(mol)
            from rdkit.Chem import rdDepictor

            rdDepictor.Compute2DCoords(mol)
        except Exception:
            pass
        mols.append(mol)
    svg = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=(300, 250),
        legends=[ascii_text(item) for item in legends],
        useSVG=True,
    )
    path.write_text(str(svg), encoding="utf-8")


def write_annotated_maegz(input_path, output_path, records):
    by_entry = {record["entry_number"]: record for record in records}
    receptor_written = False
    with structure.StructureWriter(str(output_path), overwrite=True) as writer:
        with structure.StructureReader(str(input_path)) as reader:
            for entry_number, st in enumerate(reader, start=1):
                record = by_entry.get(entry_number)
                if record is not None:
                    st.property["i_user_cluster_id"] = int(record["cluster_id"])
                    st.property["i_user_cluster_size"] = int(record["cluster_size"])
                    st.property["i_user_score_rank"] = int(record["score_rank"])
                    st.property["r_user_ligand_efficiency"] = float(
                        record["ligand_efficiency"]
                    )
                    st.property["r_user_similarity_to_medoid"] = float(record["similarity_to_medoid"])
                    st.property["i_user_pass_score_gate"] = int(record["pass_score_gate"])
                    st.property["i_user_pass_le_gate"] = int(record["pass_le_gate"])
                    st.property["i_user_pass_review_gate"] = int(record["pass_review_gate"])
                    st.property["i_user_eligible_cluster_size"] = int(record["eligible_cluster_size"])
                    for old_key in ("r_user_selection_priority", "r_user_similarity_to_centroid"):
                        if old_key in st.property:
                            del st.property[old_key]
                    st.property["s_user_review_flags"] = record["review_flags"]
                    writer.append(st)
                    continue

                if not receptor_written:
                    score, _ = score_from_structure(st)
                    atom_total = int(getattr(st, "atom_total", 0) or 0)
                    if score is None and atom_total > 100:
                        writer.append(st)
                        receptor_written = True


def write_selected_maegz(input_path, output_path, selected_records):
    wanted = {record["entry_number"]: record for record in selected_records}
    found = {}
    receptor = None
    with structure.StructureReader(str(input_path)) as reader:
        for entry_number, st in enumerate(reader, start=1):
            if entry_number in wanted:
                found[entry_number] = st.copy()
            elif receptor is None:
                score, _ = score_from_structure(st)
                atom_total = int(getattr(st, "atom_total", 0) or 0)
                if score is None and atom_total > 100:
                    receptor = st.copy()
            if len(found) == len(wanted) and receptor is not None:
                break

    with structure.StructureWriter(str(output_path), overwrite=True) as writer:
        if receptor is not None:
            writer.append(receptor)
        for record in selected_records:
            st = found.get(record["entry_number"])
            if st is None:
                continue
            st.property["i_user_selection_rank"] = int(record["selection_rank"])
            st.property["i_user_cluster_id"] = int(record["cluster_id"])
            st.property["i_user_cluster_size"] = int(record["cluster_size"])
            st.property["r_user_ligand_efficiency"] = float(
                record["ligand_efficiency"]
            )
            st.property["i_user_eligible_cluster_size"] = int(record["eligible_cluster_size"])
            for old_key in ("r_user_selection_priority", "r_user_similarity_to_centroid"):
                if old_key in st.property:
                    del st.property[old_key]
            st.property["s_user_review_flags"] = record["review_flags"]
            writer.append(st)


def main():
    args = parse_args()
    if structure is None and not args.input_csv:
        raise RuntimeError("MAEGZ input requires licensed Schrodinger Python")
    if min(args.max_molecules, args.top30, args.review_set, args.max_per_cluster,
           args.fp_size, args.n_clusters, args.score_top, args.le_top) < 1 or args.radius < 0:
        raise ValueError("Counts must be positive and radius non-negative")
    input_path = Path(args.input_csv or args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    log(f"Input: {input_path}")
    log(f"Output directory: {output_dir}")
    log(
        f"Fingerprint: Morgan radius={args.radius}, bits={args.fp_size}; "
        f"average linkage; target clusters={args.n_clusters}"
    )

    pains_catalog = None if args.skip_pains else make_pains_catalog()
    records = (read_parent_csv(input_path, pains_catalog) if args.input_csv else
               read_and_deduplicate(input_path, args.max_molecules, pains_catalog))
    if len(records) < 2:
        raise RuntimeError("Fewer than two valid scored ligand structures were found")

    fingerprints = make_fingerprints(records, args.radius, args.fp_size)
    clusters, tree = cluster_records(records, fingerprints, args.n_clusters)
    annotate_review_flags(records)
    eligible, selection_report = sequential_pool(records, args.score_top, args.le_top,
                                                 args.exclude_review_flags)
    top30 = select_candidates(eligible, args.top30, args.max_per_cluster)
    review = select_candidates(eligible, args.review_set, args.max_per_cluster)
    selection_report.update(top30_requested=args.top30, top30_returned=len(top30),
                            review_requested=args.review_set, review_returned=len(review),
                            max_per_cluster=args.max_per_cluster)
    (output_dir / "selection_report.json").write_text(
        json.dumps(selection_report, indent=2), encoding="utf-8")
    np.save(output_dir / "linkage.npy", tree)
    make_hierarchy_plots(output_dir, records, fingerprints, clusters, tree)
    plot_selection_stages(output_dir, records, selection_report, top30)
    log(f"Score pool={selection_report['score_pool']}; LE pool={selection_report['le_pool']}; "
        f"eligible clusters={selection_report['eligible_clusters']}")
    if len(review) < args.review_set or len(top30) < args.top30:
        log("SHORTFALL: insufficient eligible clusters/candidates; gates and cap were NOT relaxed.")
    from rdkit import rdBase
    import scipy
    from prepare_top_hits import file_sha256
    metadata = {"arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, "python": sys.version,
                "numpy": np.__version__, "rdkit": rdBase.rdkitVersion,
                "scipy": scipy.__version__,
                "unique_parents": len(records), "clusters": len(clusters),
                "method": "average-linkage Tanimoto; exact cut_tree",
                "input_sha256": file_sha256(input_path),
                "selection": "sequential gates; no weighted priority"}
    (output_dir / "clustering_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8")

    summaries = cluster_summary(records)
    write_records_csv(output_dir / "cluster_assignments.csv", records)
    write_records_csv(output_dir / "eligible_pool.csv", eligible)
    representatives = select_candidates(eligible, len(records), 1)
    write_records_csv(output_dir / "eligible_cluster_representatives.csv", representatives, True)
    write_records_csv(output_dir / "stage1_score_pool.csv", [r for r in records if r["pass_score_gate"]])
    write_records_csv(output_dir / "stage2_le_pool.csv", [r for r in records if r["pass_le_gate"]])
    write_summary_csv(output_dir / "cluster_summary.csv", summaries)
    make_overview_plot(output_dir / "cluster_overview.png", records, summaries)

    centroids = [
        next(records[i] for i in group if records[i]["is_cluster_medoid"])
        for group in clusters
    ]
    centroid_legends = [
        (
            f"C{record['cluster_id']} n={record['cluster_size']} "
            f"medoid {record['title']}"
        )
        for record in centroids
    ]
    write_structure_grid_svg(
        output_dir / "cluster_medoids.svg",
        centroids,
        centroid_legends,
    )

    write_records_csv(
        output_dir / "diverse_top30.csv", top30, include_selection_rank=True
    )
    write_records_csv(
        output_dir / "diverse_top60_for_review.csv",
        review,
        include_selection_rank=True,
    )
    write_structure_grid_svg(
        output_dir / "diverse_top30.svg",
        top30,
        [
            (
                f"#{record['selection_rank']} C{record['cluster_id']} "
                f"{record['title']} DS={record['docking_score']:.2f} "
                f"LE={record['ligand_efficiency']:.3f}"
            )
            for record in top30
        ],
    )

    if not args.input_csv:
        log("Writing annotated Maestro files")
        write_annotated_maegz(input_path, output_dir / "clustered_top1w.maegz", records)
        write_selected_maegz(input_path, output_dir / "diverse_top30.maegz", top30)
        write_selected_maegz(input_path, output_dir / "diverse_top60_for_review.maegz", review)
    else:
        log("CSV mode: original entry numbers retained; no MAEGZ files generated.")

    log("Finished. Main outputs:")
    for filename in (
        "cluster_assignments.csv",
        "cluster_summary.csv",
        "cluster_overview.png",
        "cluster_medoids.svg",
        "clustered_top1w.maegz",
        "diverse_top30.csv",
        "diverse_top30.svg",
        "diverse_top30.maegz",
        "diverse_top60_for_review.csv",
        "diverse_top60_for_review.maegz",
    ):
        if (output_dir / filename).exists():
            log(f"  {output_dir / filename}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
