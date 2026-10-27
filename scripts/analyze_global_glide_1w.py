#!/usr/bin/env python3
"""Analyze the top-10,000 Glide SP poses and build global pose-density maps.

Run on a host with Schrodinger installed:

    $SCHRODINGER/run python3 analyze_global_glide_1w.py

The default input is:
    data/private/truncated_hits..maegz

Outputs include a CSV table, three requested plots, one additional LE-HAC plot,
and normalized OpenDX maps for ligand-center and ligand-heavy-atom density.
DockingScore axes intentionally do not display a unit.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from schrodinger import structure
except ImportError:
    structure = None


DEFAULT_INPUT = Path(
    "data/private/truncated_hits..maegz"
)
DEFAULT_OUTPUT = Path(
    "results/global_analysis_1w"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze top-10,000 Glide SP poses and generate DX density maps."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--id-property",
        default=None,
        help=(
            "Maestro property containing the original parent-compound ID. "
            "If omitted, structure title is used."
        ),
    )
    parser.add_argument(
        "--keep-all-poses",
        action="store_true",
        help="Do not collapse records to the best DockingScore per parent ID.",
    )
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--sigma", type=float, default=2.5)
    parser.add_argument("--padding", type=float, default=5.0)
    parser.add_argument(
        "--skip-density",
        action="store_true",
        help="Generate CSV and plots only; skip the OpenDX maps.",
    )
    return parser.parse_args()


def get_float_property(st: Any, name: str) -> Optional[float]:
    value = st.property.get(name)
    if value is None:
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def heavy_atom_coordinates(st: Any) -> np.ndarray:
    coords = []
    for atom in st.atom:
        element = str(atom.element).strip().upper()
        if element in {"H", "D"}:
            continue
        coords.append((atom.x, atom.y, atom.z))
    if not coords:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(coords, dtype=np.float64)


def parent_id(st: Any, pose_index: int, id_property: Optional[str]) -> str:
    if id_property:
        value = st.property.get(id_property)
        if value not in (None, ""):
            return str(value)
    if st.title:
        return str(st.title)
    return f"pose_{pose_index}"


def read_pose_file(
    path: Path, id_property: Optional[str]
) -> Tuple[Optional[np.ndarray], List[Dict[str, Any]]]:
    receptor_coords = None
    records: List[Dict[str, Any]] = []
    first_pose_properties: Optional[List[str]] = None

    for index, st in enumerate(structure.StructureReader(str(path)), 1):
        docking_score = get_float_property(st, "r_i_docking_score")

        # A PoseViewer receptor normally has no DockingScore.
        if docking_score is None:
            coords = heavy_atom_coordinates(st)
            if receptor_coords is None and len(coords) > 100:
                receptor_coords = coords
            continue

        coords = heavy_atom_coordinates(st)
        if len(coords) == 0:
            continue

        if first_pose_properties is None:
            first_pose_properties = sorted(st.property.keys())

        hac = int(len(coords))
        records.append(
            {
                "source_index": index,
                "parent_id": parent_id(st, index, id_property),
                "title": str(st.title or ""),
                "docking_score": docking_score,
                "glide_gscore": get_float_property(st, "r_i_glide_gscore"),
                "emodel": get_float_property(st, "r_i_glide_emodel"),
                "evdw": get_float_property(st, "r_i_glide_evdw"),
                "ecoul": get_float_property(st, "r_i_glide_ecoul"),
                "hbond": get_float_property(st, "r_i_glide_hbond"),
                "lipo": get_float_property(st, "r_i_glide_lipo"),
                "state_penalty": get_float_property(
                    st, "r_epik_State_Penalty"
                ),
                "hac": hac,
                "le": -docking_score / hac,
                "coords": coords,
                # Geometric center of heavy atoms, not a mass-weighted COM.
                "center": coords.mean(axis=0),
            }
        )

    if first_pose_properties:
        print("\nProperties present on the first scored pose:")
        for prop in first_pose_properties:
            print(f"  {prop}")

    return receptor_coords, records


def choose_best_per_parent(
    records: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Keep the lowest DockingScore per parent; break ties with Emodel."""
    best: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        old = best.get(rec["parent_id"])
        new_emodel = rec["emodel"] if rec["emodel"] is not None else math.inf
        new_key = (rec["docking_score"], new_emodel)
        if old is None:
            best[rec["parent_id"]] = rec
            continue
        old_emodel = old["emodel"] if old["emodel"] is not None else math.inf
        old_key = (old["docking_score"], old_emodel)
        if new_key < old_key:
            best[rec["parent_id"]] = rec
    return list(best.values())


def write_records_csv(records: List[Dict[str, Any]], output_file: Path) -> None:
    fields = [
        "rank",
        "parent_id",
        "title",
        "source_index",
        "docking_score",
        "glide_gscore",
        "emodel",
        "evdw",
        "ecoul",
        "hbond",
        "lipo",
        "state_penalty",
        "hac",
        "le",
        "center_x",
        "center_y",
        "center_z",
    ]
    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, rec in enumerate(records, 1):
            writer.writerow(
                {
                    "rank": rank,
                    "parent_id": rec["parent_id"],
                    "title": rec["title"],
                    "source_index": rec["source_index"],
                    "docking_score": rec["docking_score"],
                    "glide_gscore": rec["glide_gscore"],
                    "emodel": rec["emodel"],
                    "evdw": rec["evdw"],
                    "ecoul": rec["ecoul"],
                    "hbond": rec["hbond"],
                    "lipo": rec["lipo"],
                    "state_penalty": rec["state_penalty"],
                    "hac": rec["hac"],
                    "le": rec["le"],
                    "center_x": rec["center"][0],
                    "center_y": rec["center"][1],
                    "center_z": rec["center"][2],
                }
            )


def percentile_summary(records: List[Dict[str, Any]], output_file: Path) -> None:
    scores = np.asarray([r["docking_score"] for r in records], dtype=float)
    hac = np.asarray([r["hac"] for r in records], dtype=float)
    le = np.asarray([r["le"] for r in records], dtype=float)
    with output_file.open("w") as handle:
        handle.write(f"Scored records: {len(records)}\n")
        handle.write(f"DockingScore min: {scores.min():.6f}\n")
        handle.write(f"DockingScore mean: {scores.mean():.6f}\n")
        handle.write(f"DockingScore median: {np.median(scores):.6f}\n")
        handle.write(f"DockingScore max: {scores.max():.6f}\n")
        handle.write(f"HAC mean: {hac.mean():.3f}\n")
        handle.write(f"HAC median: {np.median(hac):.3f}\n")
        handle.write(f"LE mean: {le.mean():.6f}\n")
        handle.write(f"LE median: {np.median(le):.6f}\n")
        handle.write("\nDockingScore percentiles:\n")
        for p in (0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100):
            handle.write(f"  p{p}: {np.percentile(scores, p):.6f}\n")
        handle.write("\nLE percentiles:\n")
        for p in (0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100):
            handle.write(f"  p{p}: {np.percentile(le, p):.6f}\n")


def plot_histogram(records: List[Dict[str, Any]], output_file: Path) -> None:
    scores = np.asarray([r["docking_score"] for r in records])
    plt.figure(figsize=(7.2, 5.4))
    plt.hist(scores, bins=60, color="#4C78A8", alpha=0.86, edgecolor="white")
    plt.axvline(
        np.median(scores),
        color="#D62728",
        linestyle="--",
        linewidth=1.7,
        label=f"Median = {np.median(scores):.2f}",
    )
    plt.xlabel("DockingScore")
    plt.ylabel("Number of analyzed records")
    plt.title("Distribution of Glide SP DockingScore")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_score_vs_hac(records: List[Dict[str, Any]], output_file: Path) -> None:
    hac = np.asarray([r["hac"] for r in records])
    scores = np.asarray([r["docking_score"] for r in records])
    plt.figure(figsize=(7.2, 5.4))
    plt.scatter(hac, scores, s=10, color="#4C78A8", alpha=0.28, linewidths=0)
    if len(records) > 2 and len(np.unique(hac)) > 1:
        slope, intercept = np.polyfit(hac, scores, 1)
        xline = np.linspace(hac.min(), hac.max(), 100)
        plt.plot(
            xline,
            slope * xline + intercept,
            color="#D62728",
            linewidth=1.7,
            label=f"Linear fit; slope = {slope:.3f}",
        )
        plt.legend(frameon=False)
    plt.xlabel("Heavy atom count (HAC)")
    plt.ylabel("DockingScore")
    plt.title("DockingScore vs Heavy Atom Count")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_score_vs_le(records: List[Dict[str, Any]], output_file: Path) -> None:
    scores = np.asarray([r["docking_score"] for r in records])
    le = np.asarray([r["le"] for r in records])
    plt.figure(figsize=(7.2, 5.4))
    plt.scatter(scores, le, s=10, color="#59A14F", alpha=0.30, linewidths=0)
    plt.xlabel("DockingScore")
    plt.ylabel("Ligand efficiency (-DockingScore / HAC)")
    plt.title("DockingScore vs Ligand Efficiency")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_le_vs_hac(records: List[Dict[str, Any]], output_file: Path) -> None:
    hac = np.asarray([r["hac"] for r in records])
    le = np.asarray([r["le"] for r in records])
    plt.figure(figsize=(7.2, 5.4))
    plt.scatter(hac, le, s=10, color="#F28E2B", alpha=0.30, linewidths=0)
    plt.xlabel("Heavy atom count (HAC)")
    plt.ylabel("Ligand efficiency (-DockingScore / HAC)")
    plt.title("Ligand Efficiency vs Heavy Atom Count")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def grid_definition(
    receptor_coords: Optional[np.ndarray],
    records: List[Dict[str, Any]],
    spacing: float,
    padding: float,
) -> Tuple[np.ndarray, np.ndarray]:
    # Include ligand extents: otherwise points outside the receptor bounding box
    # are silently clipped to an edge voxel.
    parts = [r["coords"] for r in records]
    if receptor_coords is not None:
        parts.append(receptor_coords)
    reference = np.vstack(parts)
    origin = reference.min(axis=0) - padding
    upper = reference.max(axis=0) + padding
    counts = np.ceil((upper - origin) / spacing).astype(int) + 1
    return origin, counts


def voxel_indices(
    points: np.ndarray,
    origin: np.ndarray,
    counts: np.ndarray,
    spacing: float,
) -> np.ndarray:
    indices = np.rint((points - origin) / spacing).astype(int)
    return np.clip(indices, 0, counts - 1)


def gaussian_smooth(density: np.ndarray, sigma_voxels: float) -> np.ndarray:
    """Use SciPy when available, with a NumPy separable fallback."""
    try:
        from scipy.ndimage import gaussian_filter

        print("Density smoothing backend: scipy.ndimage.gaussian_filter")
        return gaussian_filter(density, sigma=sigma_voxels, mode="constant")
    except ImportError:
        print("SciPy unavailable; using NumPy separable Gaussian convolution")
        radius = max(1, int(math.ceil(3.0 * sigma_voxels)))
        x = np.arange(-radius, radius + 1, dtype=float)
        kernel = np.exp(-(x * x) / (2.0 * sigma_voxels * sigma_voxels))
        kernel /= kernel.sum()
        result = density.astype(np.float64, copy=True)
        for axis in range(3):
            result = np.apply_along_axis(
                lambda values: np.convolve(
                    np.pad(values, (radius, radius)), kernel, mode="valid"),
                axis,
                result,
            )
        return result


def normalize_density(density: np.ndarray) -> np.ndarray:
    maximum = float(density.max())
    if maximum > 0:
        density = density / maximum
    return density.astype(np.float32)


def make_centroid_density(
    records: List[Dict[str, Any]],
    origin: np.ndarray,
    counts: np.ndarray,
    spacing: float,
    sigma: float,
) -> np.ndarray:
    density = np.zeros(tuple(counts), dtype=np.float32)
    centers = np.vstack([r["center"] for r in records])
    indices = voxel_indices(centers, origin, counts, spacing)
    np.add.at(density, tuple(indices.T), 1.0)
    return normalize_density(gaussian_smooth(density, sigma / spacing))


def make_heavy_atom_density(
    records: List[Dict[str, Any]],
    origin: np.ndarray,
    counts: np.ndarray,
    spacing: float,
    sigma: float,
) -> np.ndarray:
    density = np.zeros(tuple(counts), dtype=np.float32)
    for rec in records:
        indices = voxel_indices(rec["coords"], origin, counts, spacing)
        # Each ligand contributes a total weight of one, independent of HAC.
        weights = np.full(len(indices), 1.0 / len(indices), dtype=np.float32)
        np.add.at(density, tuple(indices.T), weights)
    return normalize_density(gaussian_smooth(density, sigma / spacing))


def write_dx(
    output_file: Path,
    density: np.ndarray,
    origin: np.ndarray,
    spacing: float,
) -> None:
    nx, ny, nz = density.shape
    total = nx * ny * nz
    with output_file.open("w") as handle:
        handle.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        handle.write(
            f"origin {origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}\n"
        )
        handle.write(f"delta {spacing:.6f} 0 0\n")
        handle.write(f"delta 0 {spacing:.6f} 0\n")
        handle.write(f"delta 0 0 {spacing:.6f}\n")
        handle.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        handle.write(
            f"object 3 class array type double rank 0 items {total} data follows\n"
        )
        flat = density.ravel(order="C")
        for start in range(0, len(flat), 3):
            handle.write(" ".join(f"{v:.7e}" for v in flat[start : start + 3]))
            handle.write("\n")
        handle.write('attribute "dep" string "positions"\n')
        handle.write('object "density" class field\n')
        handle.write('component "positions" value 1\n')
        handle.write('component "connections" value 2\n')
        handle.write('component "data" value 3\n')


def main() -> None:
    args = parse_args()
    if structure is None:
        raise RuntimeError("MAEGZ input requires licensed Schrodinger Python")
    args.output.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if args.spacing <= 0 or args.sigma <= 0 or args.padding < 0:
        raise ValueError("spacing and sigma must be > 0; padding must be >= 0")

    print(f"Reading: {args.input}")
    receptor_coords, raw_records = read_pose_file(args.input, args.id_property)
    if not raw_records:
        raise RuntimeError("No structures with r_i_docking_score were found.")

    records = (
        raw_records if args.keep_all_poses else choose_best_per_parent(raw_records)
    )
    records.sort(key=lambda r: r["docking_score"])

    print("\nCounts")
    print(f"  Raw scored poses: {len(raw_records)}")
    print(f"  Records used:     {len(records)}")
    print(f"  Receptor found:   {receptor_coords is not None}")
    print(f"  ID field:         {args.id_property or 'structure title'}")

    write_records_csv(records, args.output / "top1w_analyzed_poses.csv")
    percentile_summary(records, args.output / "top1w_summary.txt")
    plot_histogram(records, args.output / "docking_score_histogram.png")
    plot_score_vs_hac(records, args.output / "docking_score_vs_hac.png")
    plot_score_vs_le(records, args.output / "docking_score_vs_le.png")
    plot_le_vs_hac(records, args.output / "le_vs_hac.png")

    if not args.skip_density:
        origin, counts = grid_definition(
            receptor_coords, records, args.spacing, args.padding
        )
        grid_points = int(np.prod(counts, dtype=np.int64))
        print("\nDensity grid")
        print(f"  Origin:  {origin}")
        print(f"  Counts:  {counts}")
        print(f"  Points:  {grid_points:,}")
        print(f"  Spacing: {args.spacing}")
        print(f"  Sigma:   {args.sigma}")

        centroid_density = make_centroid_density(
            records, origin, counts, args.spacing, args.sigma
        )
        write_dx(
            args.output / "top1w_centroid_density.dx",
            centroid_density,
            origin,
            args.spacing,
        )
        del centroid_density

        heavy_atom_density = make_heavy_atom_density(
            records, origin, counts, args.spacing, args.sigma
        )
        write_dx(
            args.output / "top1w_heavyatom_density.dx",
            heavy_atom_density,
            origin,
            args.spacing,
        )
        del heavy_atom_density

    print(f"\nFinished. Outputs: {args.output}")


if __name__ == "__main__":
    main()
