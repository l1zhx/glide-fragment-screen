#!/usr/bin/env python3
"""Draw the DockingScore top-100 Glide poses as one large 2D figure.

Run on a machine with Schrodinger installed:

    $SCHRODINGER/run python3 draw_top100_2d.py

Default input:
    data/private/top100_SP_pv.maegz

The receptor is skipped automatically because it has no r_i_docking_score.
Scored entries are sorted from the lowest to the highest DockingScore, and the
first 100 are drawn. DockingScore is displayed without a unit.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from schrodinger import structure
from schrodinger.structutils import smiles as schrodinger_smiles

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
except ImportError as exc:
    raise SystemExit(
        "RDKit could not be imported. Run this script with "
        "'$SCHRODINGER/run python3', or use a Python environment containing RDKit."
    ) from exc


DEFAULT_INPUT = Path(
    "data/private/top100_SP_pv.maegz"
)
DEFAULT_OUTPUT = Path(
    "results/top100_2D.png"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw the DockingScore top-100 Glide poses in one image."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--number", type=int, default=100)
    parser.add_argument("--columns", type=int, default=10)
    parser.add_argument("--cell-width", type=int, default=380)
    parser.add_argument("--cell-height", type=int, default=300)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--label-property",
        default=None,
        help=(
            "Optional Maestro property to use as the compound label. "
            "The structure title is used by default."
        ),
    )
    parser.add_argument(
        "--score-decimals",
        type=int,
        default=2,
        help="Number of decimal places shown for DockingScore.",
    )
    parser.add_argument(
        "--deduplicate-smiles",
        action="store_true",
        help=(
            "Keep only the best-scoring occurrence of each canonical SMILES. "
            "By default, the first 100 scored entries are drawn without deduplication."
        ),
    )
    parser.add_argument(
        "--no-svg",
        action="store_true",
        help="Do not generate an additional vector SVG figure.",
    )
    return parser.parse_args()


def float_property(st: Any, name: str) -> Optional[float]:
    value = st.property.get(name)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compound_label(st: Any, index: int, label_property: Optional[str]) -> str:
    if label_property:
        value = st.property.get(label_property)
        if value not in (None, ""):
            return str(value)
    if st.title:
        return str(st.title)
    return f"pose_{index}"


def shorten(text: str, maximum: int = 30) -> str:
    text = " ".join(text.split())
    if len(text) <= maximum:
        return text
    return text[: maximum - 3] + "..."


def safe_legend_text(text: str) -> str:
    """Return text accepted by the legacy RDKit/Pillow Latin-1 renderer.

    Some Schrodinger 2024 installations route RDKit grid legends through the
    legacy PIL canvas. It cannot render characters such as the Unicode ellipsis
    or arbitrary CJK text. ASCII replacement prevents a late drawing failure.
    """
    return text.encode("ascii", errors="replace").decode("ascii")


def make_rdkit_molecule(st: Any, generator: Any) -> Tuple[Any, str]:
    """Convert a Schrodinger Structure to a sanitized RDKit molecule."""
    smi = generator.getSmiles(st)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError("RDKit could not parse the generated SMILES")
    mol = Chem.RemoveHs(mol)

    # CoordGen generally gives clean medicinal-chemistry depictions. Fall back
    # to the standard RDKit depiction engine when CoordGen is unavailable.
    try:
        from rdkit.Chem import rdCoordGen

        rdCoordGen.AddCoords(mol)
    except (ImportError, AttributeError, RuntimeError):
        rdDepictor.Compute2DCoords(mol)

    return mol, smi


def read_scored_entries(
    input_file: Path, label_property: Optional[str]
) -> Tuple[List[Dict[str, Any]], int]:
    generator = schrodinger_smiles.SmilesGenerator(
        unique=True,
        safe=True,
        wantAllH=False,
    )
    entries: List[Dict[str, Any]] = []
    skipped = 0

    for index, st in enumerate(structure.StructureReader(str(input_file)), 1):
        score = float_property(st, "r_i_docking_score")
        if score is None:
            # Normally the receptor entry in a PoseViewer file.
            continue

        try:
            mol, smi = make_rdkit_molecule(st, generator)
        except Exception as exc:
            skipped += 1
            print(
                f"Warning: skipped entry {index} ({st.title!r}): {exc}",
                file=sys.stderr,
            )
            continue

        entries.append(
            {
                "source_index": index,
                "label": compound_label(st, index, label_property),
                "score": score,
                "smiles": smi,
                "mol": mol,
            }
        )

    entries.sort(key=lambda item: item["score"])
    return entries, skipped


def deduplicate_by_smiles(entries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    unique: List[Dict[str, Any]] = []
    seen = set()
    for entry in entries:
        canonical = Chem.MolToSmiles(entry["mol"], canonical=True, isomericSmiles=True)
        if canonical in seen:
            continue
        seen.add(canonical)
        unique.append(entry)
    return unique


def drawing_options() -> Any:
    options = Draw.MolDrawOptions()
    options.addStereoAnnotation = True
    options.legendFontSize = 18
    options.minFontSize = 12
    options.maxFontSize = 28
    options.bondLineWidth = 2.0
    options.padding = 0.08
    return options


def save_png(image: Any, output_file: Path, dpi: int) -> None:
    """Save either a PIL image or PNG bytes returned by RDKit."""
    if isinstance(image, bytes):
        from PIL import Image

        image = Image.open(io.BytesIO(image))
        image.load()
    image.save(str(output_file), dpi=(dpi, dpi))


def write_csv(entries: List[Dict[str, Any]], output_file: Path) -> None:
    with output_file.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "label", "docking_score", "smiles", "source_index"],
        )
        writer.writeheader()
        for rank, entry in enumerate(entries, 1):
            writer.writerow(
                {
                    "rank": rank,
                    "label": entry["label"],
                    "docking_score": entry["score"],
                    "smiles": entry["smiles"],
                    "source_index": entry["source_index"],
                }
            )


def main() -> None:
    args = parse_args()

    print("[draw_top100_2d] Script started", flush=True)
    print(f"[draw_top100_2d] Input:  {args.input}", flush=True)
    print(f"[draw_top100_2d] Output: {args.output}", flush=True)

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if args.number <= 0 or args.columns <= 0:
        raise ValueError("--number and --columns must be positive integers")
    if args.cell_width < 100 or args.cell_height < 100:
        raise ValueError("Cell dimensions are too small")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    entries, skipped = read_scored_entries(args.input, args.label_property)
    if args.deduplicate_smiles:
        entries = deduplicate_by_smiles(entries)

    selected = entries[: args.number]
    if not selected:
        raise RuntimeError("No drawable structures with r_i_docking_score were found")

    legends = []
    for rank, entry in enumerate(selected, 1):
        legend = (
            f"{rank}. {shorten(entry['label'])}\n"
            f"DockingScore = {entry['score']:.{args.score_decimals}f}"
        )
        legends.append(safe_legend_text(legend))

    mols = [entry["mol"] for entry in selected]
    options = drawing_options()

    png_image = Draw.MolsToGridImage(
        mols,
        molsPerRow=args.columns,
        subImgSize=(args.cell_width, args.cell_height),
        legends=legends,
        useSVG=False,
        returnPNG=False,
        drawOptions=options,
    )
    save_png(png_image, args.output, args.dpi)

    if not args.no_svg:
        svg_text = Draw.MolsToGridImage(
            mols,
            molsPerRow=args.columns,
            subImgSize=(args.cell_width, args.cell_height),
            legends=legends,
            useSVG=True,
            drawOptions=options,
        )
        svg_file = args.output.with_suffix(".svg")
        if isinstance(svg_text, bytes):
            svg_text = svg_text.decode("utf-8")
        svg_file.write_text(str(svg_text), encoding="utf-8")

    csv_file = args.output.with_name(args.output.stem + "_index.csv")
    write_csv(selected, csv_file)

    rows = (len(selected) + args.columns - 1) // args.columns
    print(f"Scored/drawable structures found: {len(entries)}")
    print(f"Structures drawn:                 {len(selected)}")
    print(f"Conversion failures skipped:      {skipped}")
    print(f"Layout:                           {args.columns} columns x {rows} rows")
    print(f"PNG:                              {args.output}")
    if not args.no_svg:
        print(f"SVG:                              {args.output.with_suffix('.svg')}")
    print(f"Index CSV:                        {csv_file}")


if __name__ == "__main__":
    main()
