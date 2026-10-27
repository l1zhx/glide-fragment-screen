#!/usr/bin/env python3
"""Export the user's manually selected candidates from a saved Top60 result.

Run on the server with:
    $SCHRODINGER/run python3 export_manual_top30.py

The numbers refer to CSV selection_rank / Maestro i_user_selection_rank,
NOT cluster_id, raw file entry number, or Maestro's global row number.
Existing ranks and scores are read from disk; clustering is not repeated.
Figures show structure, name, DockingScore, LE, cluster and pocket label.
Selection ranks are retained in data files for traceability, not in figures.
"""

from __future__ import annotations

import argparse
import csv
import math
import textwrap
import json
from pathlib import Path
from prepare_top_hits import file_sha256


DEFAULT_ROWS = (
    "1,3,4,5,6,7,11,12,14,17,18,19,20,24,25,28,29,30,31,34,"
    "36,38,39,41,42,43,44,45,46,47"
)
DEFAULT_POCKETS = (
    "Ser246,Ser246,Asn326,Asn326,Asn326,Asn326,Ser246,Ser246,Ser246,Ser246,"
    "Ser246,Asn326,Ser246,Thr150,Asn326,Ser246,Thr150,Ser246,Asn326,Ser246,"
    "Asn326,Ser246,Asn326,Ser246,Lys233,Asn326,Ser246,Ser246,Ser246,Asn326"
)
DEFAULT_CSV = (
    "results/cluster_top1w/"
    "diverse_top60_for_review.csv"
)


def log(message):
    print(f"[export_manual_top30] {message}", flush=True)


def parse_rows(text):
    rows = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not rows or any(item < 1 for item in rows):
        raise ValueError("Selection ranks must be positive integers.")
    if len(rows) != len(set(rows)):
        raise ValueError("Duplicate selection ranks were provided.")
    return rows


def parse_pockets(text):
    labels = [item.strip() for item in text.replace("，", ",").split(",")]
    if not labels or any(not label for label in labels):
        raise ValueError("Pocket labels must be non-empty comma-separated values.")
    return labels


def assign_pockets(selected, pocket_labels):
    if len(selected) != len(pocket_labels):
        raise ValueError(
            f"Selected {len(selected)} molecules, but received "
            f"{len(pocket_labels)} pocket labels. Counts must match."
        )
    for record, pocket in zip(selected, pocket_labels):
        record["pocket_label"] = pocket


def choose_records(all_records, requested_rows):
    by_rank = {}
    for source_record in all_records:
        rank = int(source_record["selection_rank"])
        if rank in by_rank:
            raise ValueError(f"Duplicate selection_rank {rank} in input CSV.")
        by_rank[rank] = dict(source_record)
    missing = [rank for rank in requested_rows if rank not in by_rank]
    if missing:
        raise ValueError(f"Selection ranks missing from input CSV: {missing}")

    selected = []
    for manual_rank, top60_rank in enumerate(requested_rows, start=1):
        record = by_rank[top60_rank]
        record["manual_rank"] = manual_rank
        record["original_top60_rank"] = top60_rank
        score = float(record["docking_score"])
        le = float(record["ligand_efficiency"])
        hac = int(record["hac"])
        if hac <= 0 or not math.isfinite(score) or not math.isfinite(le):
            raise ValueError(f"Invalid score/LE/HAC at Top60 #{top60_rank}.")
        if not math.isclose(le, -score / hac, rel_tol=0, abs_tol=0.0006):
            raise ValueError(
                f"LE is inconsistent with -DockingScore/HAC at Top60 #{top60_rank}."
            )
        record["cluster_id"] = int(record["cluster_id"])
        record["cluster_size"] = int(record["cluster_size"])
        selected.append(record)
    return selected


def load_selection(csv_path, rows):
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        required = {
            "selection_rank", "title", "docking_score", "ligand_efficiency",
            "hac", "cluster_id", "cluster_size", "parent_smiles",
        }
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"CSV is missing required fields: {sorted(missing)}")
        selected = choose_records(list(reader), rows)
    return selected, fields


def docking_score(st):
    keys = ["r_i_docking_score", "r_glide_gscore", "r_glide_docking_score"]
    keys.extend(
        key for key in st.property.keys()
        if "docking_score" in key.lower() or key.lower().endswith("glide_gscore")
    )
    for key in keys:
        try:
            value = float(st.property[key])
            if math.isfinite(value):
                return value
        except (KeyError, ValueError, TypeError):
            pass
    return None


def molecules_from_maegz(path, selected):
    try:
        from schrodinger import structure
    except ImportError as exc:
        raise RuntimeError(
            "Use $SCHRODINGER/run python3 for MAEGZ input, or --csv-only "
            "to draw standardized parent SMILES without a Maestro file."
        ) from exc
    from rdkit import Chem

    try:
        from schrodinger.rdkit import rdkit_adapter
    except ImportError:
        rdkit_adapter = None

    wanted = {int(item["original_top60_rank"]): item for item in selected}
    found = {}
    receptor = None
    annotated_seen = set()
    with structure.StructureReader(str(path)) as reader:
        for st in reader:
            try:
                rank = int(st.property["i_user_selection_rank"])
            except (KeyError, TypeError, ValueError):
                if receptor is None and docking_score(st) is None:
                    if int(getattr(st, "atom_total", 0) or 0) > 100:
                        receptor = st.copy()
                continue
            if rank in annotated_seen:
                raise ValueError(f"Duplicate i_user_selection_rank {rank} in MAEGZ.")
            annotated_seen.add(rank)
            if rank not in wanted:
                continue

            record = wanted[rank]
            actual_score = docking_score(st)
            if actual_score is None or not math.isclose(
                actual_score, float(record["docking_score"]),
                rel_tol=0, abs_tol=0.0006,
            ):
                raise ValueError(
                    f"CSV/MAEGZ score mismatch at Top60 #{rank}. "
                    "Use the CSV and MAEGZ from the SAME clustering run."
                )
            try:
                cluster_id = int(st.property["i_user_cluster_id"])
            except (KeyError, ValueError, TypeError) as exc:
                raise ValueError(f"MAEGZ has no cluster ID at Top60 #{rank}.") from exc
            if cluster_id != int(record["cluster_id"]):
                raise ValueError(f"CSV/MAEGZ cluster mismatch at Top60 #{rank}.")

            if rdkit_adapter is not None:
                mol = rdkit_adapter.to_rdkit(
                    st, implicitH=True, include_properties=False,
                    sanitize=True, include_stereo=True,
                )
            else:
                from schrodinger import adapter

                mol = adapter.to_rdkit(st)
            if mol is None:
                raise ValueError(f"Cannot convert molecule at Top60 #{rank}.")
            mol = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
            record["depiction_source"] = "original_docked_structure"
            record["depicted_smiles"] = Chem.MolToSmiles(mol, isomericSmiles=True)
            actual_title = str(getattr(st, "title", "") or "").strip()
            record["display_title"] = actual_title or record["title"]
            if actual_title and actual_title != record["title"]:
                log(f"Top60 #{rank}: retaining modified MAEGZ title {actual_title!r}.")
            if mol.GetNumHeavyAtoms() != int(record["hac"]):
                log(
                    f"Top60 #{rank}: depicted species has {mol.GetNumHeavyAtoms()} "
                    f"heavy atoms; saved LE uses parent HAC={record['hac']}."
                )
            found[rank] = (mol, st.copy())

    missing = [rank for rank in wanted if rank not in found]
    if missing:
        raise ValueError(
            f"Selected ranks absent from annotated Top60 MAEGZ: {missing}. "
            "Input must be diverse_top60_for_review.maegz, not the raw Top1w file."
        )
    mols = [found[int(item["original_top60_rank"])][0] for item in selected]
    structures = [found[int(item["original_top60_rank"])][1] for item in selected]
    return mols, structures, receptor


def molecules_from_csv(selected):
    from rdkit import Chem

    mols = []
    for item in selected:
        mol = Chem.MolFromSmiles(item["parent_smiles"])
        if mol is None:
            raise ValueError(
                f"Invalid parent SMILES at Top60 #{item['original_top60_rank']}."
            )
        item["depiction_source"] = "standardized_parent_smiles"
        item["depicted_smiles"] = item["parent_smiles"]
        item["display_title"] = item["title"]
        mols.append(mol)
    return mols


def make_legends(selected):
    legends = []
    for item in selected:
        # ASCII avoids the Unicode font issue encountered in the old renderer.
        title = str(item["display_title"]).encode("ascii", "replace").decode("ascii")
        name_lines = textwrap.wrap(" ".join(title.split()), width=34) or ["Unnamed"]
        if len(name_lines) > 3:
            name_lines = name_lines[:3]
            name_lines[-1] = name_lines[-1][:31] + "..."
        legends.append(
            "\n".join(name_lines)
            + f"\nDockingScore {float(item['docking_score']):.3f}"
            + f" | LE {float(item['ligand_efficiency']):.3f}"
            + f"\nCluster {item['cluster_id']} | Pocket "
            + str(item["pocket_label"]).encode("ascii", "replace").decode("ascii")
        )
    return legends


def configure_drawer(drawer, scale):
    options = drawer.drawOptions()
    options.legendFontSize = int(17 * scale)
    if hasattr(options, "legendFraction"):
        options.legendFraction = 0.32
    if hasattr(options, "bondLineWidth"):
        options.bondLineWidth = 1.6 * scale
    if hasattr(options, "maxFontSize"):
        options.maxFontSize = int(40 * scale)
    if hasattr(options, "minFontSize"):
        options.minFontSize = int(6 * scale)
    options.padding = 0.08


def draw_grid(mols, legends, prefix, columns, tile_width, tile_height, png_scale):
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor

    copies = []
    for mol in mols:
        copy = Chem.Mol(mol)
        copy.RemoveAllConformers()
        rdDepictor.Compute2DCoords(copy)
        copies.append(copy)
    rows = math.ceil(len(copies) / columns)
    width, height = columns * tile_width, rows * tile_height
    svg_path = prefix.with_suffix(".svg")
    png_path = prefix.with_suffix(".png")

    try:
        from rdkit.Chem.Draw import rdMolDraw2D

        svg_drawer = rdMolDraw2D.MolDraw2DSVG(width, height, tile_width, tile_height)
        configure_drawer(svg_drawer, 1)
        svg_drawer.DrawMolecules(copies, legends=legends)
        svg_drawer.FinishDrawing()
        svg = svg_drawer.GetDrawingText()
    except (ImportError, AttributeError):
        svg = Draw.MolsToGridImage(
            copies, molsPerRow=columns, subImgSize=(tile_width, tile_height),
            legends=legends, useSVG=True,
        )
    svg_path.write_text(str(svg), encoding="utf-8")
    log(f"SVG: {svg_path}")

    png_bytes = None
    try:
        from rdkit.Chem.Draw import rdMolDraw2D

        drawer = rdMolDraw2D.MolDraw2DCairo(
            width * png_scale, height * png_scale,
            tile_width * png_scale, tile_height * png_scale,
        )
        configure_drawer(drawer, png_scale)
        drawer.DrawMolecules(copies, legends=legends)
        drawer.FinishDrawing()
        png_bytes = drawer.GetDrawingText()
    except Exception as exc:
        log(f"Native PNG renderer unavailable ({exc}); trying SVG conversion.")

    if png_bytes is None:
        try:
            import cairosvg

            png_bytes = cairosvg.svg2png(
                bytestring=str(svg).encode("utf-8"),
                output_width=width * png_scale,
                output_height=height * png_scale,
            )
        except Exception:
            # Compatibility fallback: all legend strings above are ASCII.
            try:
                png_image = Draw.MolsToGridImage(
                    copies, molsPerRow=columns,
                    subImgSize=(tile_width * png_scale, tile_height * png_scale),
                    legends=legends, useSVG=False,
                )
                if hasattr(png_image, "save"):
                    png_image.save(str(png_path))
                    log(f"PNG: {png_path}")
                    return
                if isinstance(png_image, bytes):
                    png_bytes = png_image
                else:
                    raise RuntimeError("Unsupported image object returned by RDKit.")
            except Exception as exc:
                log(
                    f"PNG could not be rendered: {exc}. "
                    f"The completed SVG is available at {svg_path}."
                )
                return
    png_path.write_bytes(png_bytes)
    log(f"PNG: {png_path}")


def save_csv(path, selected, source_fields):
    fields = ["manual_rank", "original_top60_rank"] + source_fields
    fields += ["pocket_label", "display_title", "depiction_source", "depicted_smiles"]
    fields = list(dict.fromkeys(fields))
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    log(f"CSV: {path}")


def save_maegz(path, selected, structures, receptor):
    from schrodinger import structure

    with structure.StructureWriter(str(path), overwrite=True) as writer:
        if receptor is not None:
            writer.append(receptor)
        for item, st in zip(selected, structures):
            # Preserve old selection_rank so the saved Top60 numbering is traceable.
            st.property["i_user_manual_rank"] = int(item["manual_rank"])
            st.property["i_user_original_top60_rank"] = int(item["original_top60_rank"])
            st.property["s_user_pocket_label"] = item["pocket_label"]
            writer.append(st)
    log(f"MAEGZ: {path}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument(
        "--maegz", default=None,
        help="Annotated Top60 MAEGZ. Default: diverse_top60_for_review.maegz beside CSV.",
    )
    parser.add_argument(
        "--rows", default=None,
        help="Comma-separated original Top60 selection_rank values, in display order.",
    )
    parser.add_argument(
        "--pockets", default=None,
        help=(
            "Comma-separated pocket labels in the same order as --rows. "
            "The built-in 30 labels apply only to the original default selection."
        ),
    )
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--use-archived-selection", action="store_true",
                        help="Apply historical 30 ranks/pockets ONLY to their original archived Top60")
    parser.add_argument("--columns", type=int, default=5)
    parser.add_argument("--tile-width", type=int, default=480)
    parser.add_argument("--tile-height", type=int, default=430)
    parser.add_argument("--png-scale", type=int, default=2)
    parser.add_argument(
        "--csv-only", action="store_true",
        help="Draw standardized parent SMILES; do not read or export MAEGZ.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.use_archived_selection:
        if args.rows is not None or args.pockets is not None:
            raise ValueError("Use either archived selection or explicit --rows/--pockets")
        args.rows, args.pockets = DEFAULT_ROWS, DEFAULT_POCKETS
    elif args.rows is None or args.pockets is None:
        raise ValueError("New Top60 needs explicitly reviewed --rows and --pockets. "
                         "For the ORIGINAL historical Top60 only, use --use-archived-selection.")
    if min(args.columns, args.tile_width, args.tile_height, args.png_scale) < 1:
        raise ValueError("Columns, tile dimensions and PNG scale must be positive.")
    csv_path = Path(args.csv).expanduser().resolve()
    if not csv_path.is_file():
        raise FileNotFoundError(f"Saved Top60 CSV not found: {csv_path}")
    requested_rows = parse_rows(args.rows)
    if args.pockets is not None:
        pocket_labels = parse_pockets(args.pockets)
    elif requested_rows == parse_rows(DEFAULT_ROWS):
        pocket_labels = parse_pockets(DEFAULT_POCKETS)
    else:
        raise ValueError(
            "When changing --rows, also provide --pockets for the new selection."
        )
    selected, source_fields = load_selection(csv_path, requested_rows)
    assign_pockets(selected, pocket_labels)
    log(f"Selected {len(selected)} molecules by saved Top60 selection_rank.")
    log("Original Top60 ranks: " + ",".join(str(item["original_top60_rank"]) for item in selected))
    for item in selected:
        log(
            f"Selected {item['manual_rank']:02d} <- Top60 #{item['original_top60_rank']}: "
            f"{item['title']} | DockingScore {float(item['docking_score']):.3f} "
            f"| LE {float(item['ligand_efficiency']):.3f} "
            f"| Cluster {item['cluster_id']} | Pocket {item['pocket_label']}"
        )

    maegz_path = (
        Path(args.maegz).expanduser().resolve() if args.maegz
        else csv_path.with_name("diverse_top60_for_review.maegz")
    )
    structures = receptor = None
    if args.csv_only:
        mols = molecules_from_csv(selected)
        log("Drawing standardized parents; these may differ from docked protonation/tautomer states.")
    else:
        if not maegz_path.is_file():
            raise FileNotFoundError(
                f"Annotated Top60 MAEGZ not found: {maegz_path}. "
                "Supply --maegz or use --csv-only for a parent-structure image."
            )
        mols, structures, receptor = molecules_from_maegz(maegz_path, selected)
        log("Drawing the original docked chemical species from MAEGZ.")

    output_dir = (
        Path(args.output_dir).expanduser().resolve() if args.output_dir
        else csv_path.parent / "manual_selected30"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"manual_selected{len(selected)}"
    provenance = {"csv_sha256": file_sha256(csv_path),
                  "maegz_sha256": None if args.csv_only else file_sha256(maegz_path),
                  "source_ranks": requested_rows, "pocket_labels": pocket_labels,
                  "depiction": "parent_smiles" if args.csv_only else "docked_species"}
    (output_dir / "selection_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")
    save_csv(prefix.with_suffix(".csv"), selected, source_fields)
    if structures is not None:
        save_maegz(prefix.with_suffix(".maegz"), selected, structures, receptor)
    draw_grid(
        mols, make_legends(selected), prefix,
        args.columns, args.tile_width, args.tile_height, args.png_scale,
    )
    log("Finished. Scores and cluster assignments were not recalculated.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
