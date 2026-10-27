#!/usr/bin/env python3
"""Plot the archived researcher-supplied labels, not automated binding predictions."""
import argparse
import json
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/manual_selection.json"))
    parser.add_argument("--output", type=Path, default=Path("results/manual_pocket_counts.png"))
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if len(config["ranks"]) != len(config["pockets"]) or len(set(config["ranks"])) != len(config["ranks"]):
        raise ValueError("Invalid rank-to-pocket mapping")
    counts = Counter(config["pockets"])
    ordered = counts.most_common()
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    bars = ax.bar([p for p, _ in ordered], [n for _, n in ordered],
                  color=["#377E9D", "#57A98D", "#E0A44C", "#9C80AE"])
    ax.bar_label(bars, padding=3)
    ax.set(xlabel="Manually assigned pocket anchor", ylabel="Selected candidates",
           ylim=(0, max(counts.values()) * 1.25),
           title=f"Researcher-reported manual shortlist (n={sum(counts.values())})")
    ax.text(.5, .95, "Residue labels are annotations; binding is not established",
            transform=ax.transAxes, ha="center", fontsize=9, color="#555555")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(dict(ordered))


if __name__ == "__main__":
    main()
