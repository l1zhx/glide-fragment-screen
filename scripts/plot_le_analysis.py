#!/usr/bin/env python3
"""Scientific LE figures. Experimental theory never converts docking scores to Kd.

--theory: formula-derived figures only (298.15 K, C0=1 M), no molecular results.
--csv: molecular score/HAC/pocket table, e.g. updated EF assignments or cluster CSV.
The plot labels LE_score and DockingScore deliberately omit physical units.
"""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
import numpy as np

R_KCAL = 0.00198720425864083
T_K = 298.15
ALPHA = R_KCAL * T_K * math.log(10)


def experimental_kd(le, hac, temperature=T_K):
    if hac <= 0 or temperature <= 0 or not all(math.isfinite(v) for v in (le, hac, temperature)):
        raise ValueError("Finite LE and positive HAC/temperature required")
    return 10 ** (-le * hac / (R_KCAL * temperature * math.log(10)))


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    return plt


def save(fig, path, plt):
    fig.savefig(path.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def theory(out):
    plt = style()
    x = np.linspace(.1, .6, 101)
    fig, ax = plt.subplots(figsize=(9, 5.4), layout="constrained")
    for hac in [10, 15, 20, 25, 30]:
        ax.plot(x, x*hac/ALPHA, label=f"HAC = {hac}")
    ax.set(xlabel="Experimental LE (kcal mol$^{-1}$ heavy-atom$^{-1}$)",
           ylabel="pKd = -log10(Kd / 1 M)",
           title="Thermodynamic relation at 298.15 K, standard state 1 M")
    ax.legend(ncol=2)
    ax.text(.02, .98, "Formula-derived curves; experimental LE only", transform=ax.transAxes,
            va="top", fontsize=10)
    save(fig, out / "experimental_LE_to_pKd", plt)
    has = [10, 15, 20, 25, 30]
    les = [.2, .3, .4, .5]
    rows = [{"HAC": h, "experimental_LE": le, "temperature_K": T_K,
             "Kd_M": experimental_kd(le, h), "pKd": le*h/ALPHA} for h in has for le in les]
    with (out / "experimental_conversion.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Experimental conversion constant: {ALPHA:.9f}")


def first(row, keys):
    for k in keys:
        if row.get(k) not in (None, ""):
            return row[k]
    return None


def load_records(path):
    best, raw, duplicate = {}, 0, 0
    with path.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            raw += 1
            score_text = first(row, ["docking_score", "r_i_docking_score", "DockingScore"])
            hac_text = first(row, ["hac", "HAC", "heavy_atom_count"])
            identity = first(row, ["parent_smiles", "parent_id", "molecule_id", "compound_id"])
            if identity is None or score_text is None:
                raise ValueError(f"Row {raw}: need a parent identity and DockingScore")
            score = float(score_text)
            if hac_text is None:
                smiles = row.get("parent_smiles")
                if not smiles:
                    raise ValueError("Missing HAC: rerun the updated EF script or supply parent_smiles")
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Row {raw}: invalid parent_smiles")
                hac = mol.GetNumHeavyAtoms()
            else:
                hac_float = float(hac_text)
                if not math.isfinite(hac_float) or not hac_float.is_integer():
                    raise ValueError(f"Row {raw}: HAC must be a positive integer")
                hac = int(hac_float)
            if hac <= 0 or not math.isfinite(score):
                raise ValueError(f"Row {raw}: invalid score or HAC")
            item = dict(parent_id=identity, title=row.get("title", identity),
                        docking_score=score, hac=hac, LE_score=-score/hac,
                        pocket=first(row, ["pocket", "pocket_label"]) or "Unassigned")
            if identity in best:
                duplicate += 1
            if identity not in best or score < best[identity]["docking_score"]:
                best[identity] = item
    if not best:
        raise ValueError("No molecular records")
    return list(best.values()), dict(input_rows=raw, retained_parents=len(best), removed_duplicates=duplicate)


def pareto_indices(records):
    """Minimize BOTH HAC and DockingScore. This is a ranking aid, not affinity."""
    result = []
    for i, r in enumerate(records):
        if not any((o["hac"] <= r["hac"] and o["docking_score"] <= r["docking_score"]
                    and (o["hac"] < r["hac"] or o["docking_score"] < r["docking_score"]))
                   for o in records):
            result.append(i)
    return result


def docking_figures(records, audit, out):
    plt = style()
    pockets = list(dict.fromkeys(r["pocket"] for r in records))
    palette = plt.get_cmap("tab10")
    counts = Counter(r["pocket"] for r in records)
    for metric, filename, y_label in [("docking_score", "score_vs_HAC", "DockingScore"),
                                     ("LE_score", "LE_score_vs_HAC", "LE_score = -DockingScore / HAC")]:
        fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
        for index, pocket in enumerate(pockets):
            subset = [r for r in records if r["pocket"] == pocket]
            ax.scatter([r["hac"] for r in subset], [r[metric] for r in subset], s=12,
                       alpha=.4, linewidths=0, color=palette(index % 10), label=f"{pocket} (n={len(subset)})")
        if metric == "docking_score":
            x = np.array([min(r["hac"] for r in records), max(r["hac"] for r in records)])
            for value in [.2, .3, .4]:
                ax.plot(x, -value*x, linestyle="--", linewidth=1, alpha=.6,
                        label=f"LE_score = {value:.1f}")
        ax.set(xlabel="Heavy atom count (HAC)", ylabel=y_label,
               title="Docking data: size-normalized score analysis")
        ax.legend(fontsize=8, ncol=2)
        save(fig, out / filename, plt)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    for ax, metric in zip(axes, ["docking_score", "LE_score"]):
        groups = [[r[metric] for r in records if r["pocket"] == p] for p in pockets]
        ax.boxplot(groups, showfliers=False)
        ax.set_xticks(range(1, len(pockets)+1), [f"{p}\nn={counts[p]}" for p in pockets], rotation=25, ha="right")
        ax.set(ylabel=metric, title="Pocket comparison (unadjusted for HAC)")
    save(fig, out / "pocket_score_and_LE", plt)
    bands = [(1, 10), (11, 15), (16, 20), (21, 25), (26, math.inf)]
    fig, axes = plt.subplots(1, 5, figsize=(19, 5), layout="constrained", sharey=True)
    for ax, (low, high) in zip(axes, bands):
        subset = [r for r in records if low <= r["hac"] <= high]
        present = [p for p in pockets if any(r["pocket"] == p for r in subset)]
        if present:
            ax.boxplot([[r["LE_score"] for r in subset if r["pocket"] == p] for p in present], showfliers=False)
            ax.set_xticks(range(1, len(present)+1), present, rotation=65, ha="right", fontsize=8)
        else:
            ax.text(.5, .5, "No data", transform=ax.transAxes, ha="center")
            ax.set_xticks([])
        ax.set_title(f"HAC {low}-{int(high)}" if math.isfinite(high) else "HAC >= 26")
    axes[0].set_ylabel("LE_score")
    fig.suptitle("HAC-stratified pocket comparison; descriptive, no affinity conversion")
    save(fig, out / "pocket_LE_HAC_strata", plt)
    # Avoid quadratic all-pairs Pareto computation by traversing HAC groups.
    best_so_far = math.inf
    front = []
    for hac in sorted({r["hac"] for r in records}):
        group = [r for r in records if r["hac"] == hac]
        score = min(r["docking_score"] for r in group)
        if score < best_so_far:
            front.extend(r for r in group if r["docking_score"] == score)
            best_so_far = score
    with (out / "pareto_candidates.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(front)
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    ax.scatter([r["hac"] for r in records], [r["docking_score"] for r in records], s=10, alpha=.25, color="#999999")
    ax.scatter([r["hac"] for r in front], [r["docking_score"] for r in front], s=65, marker="*", color="#BC493E", label="Pareto candidates")
    ax.set(xlabel="HAC (lower favored)", ylabel="DockingScore (lower favored)", title="Score-size trade-off; affinity uncalibrated")
    ax.legend()
    save(fig, out / "score_HAC_Pareto", plt)
    audit.update(pocket_counts=dict(counts), LE_definition="-DockingScore / HAC", affinity_conversion="NONE")
    (out / "analysis_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--theory", action="store_true")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if not args.csv and not args.theory:
        parser.error("Supply --csv or --theory")
    args.out.mkdir(parents=True, exist_ok=False)
    if args.theory:
        theory(args.out)
    if args.csv:
        records, audit = load_records(args.csv)
        docking_figures(records, audit, args.out)
    print(f"Finished: {args.out}. Docking scores were NOT converted to affinity.")


if __name__ == "__main__":
    main()
