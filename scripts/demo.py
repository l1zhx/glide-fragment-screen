#!/usr/bin/env python3
"""Public synthetic example: actual RDKit clustering, invented scores, no docking."""
import argparse
import json
from pathlib import Path
from rdkit import Chem, rdBase
import cluster_glide_top1w as pipeline
from export_manual_top30 import choose_records, assign_pockets, draw_grid, make_legends

TOY_SMILES = [
    "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "Oc1ccccc1",
    "Nc1ccccc1", "COc1ccccc1", "c1ccncc1", "Cc1ccncc1",
    "C1CCCCC1", "OC1CCCCC1", "C1CCNCC1", "CN1CCCCC1",
    "CC(=O)O", "CCC(=O)O", "CC(N)C(=O)O", "CCO",
    "CNC(=O)c1ccc(OC)cc1", "CCNC(=O)c1ccc(OC)cc1",
    "CCCNC(=O)c1ccc(OC)cc1", "CCCCNC(=O)c1ccc(OC)cc1",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("demo_output"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    records = []
    for index, smiles in enumerate(TOY_SMILES, 1):
        mol = pipeline.standardized_parent(Chem.MolFromSmiles(smiles))
        props = pipeline.molecular_properties(mol)
        score = -3.0 - ((index * 7) % 19) / 5  # Invented; no affinity meaning.
        records.append(dict(entry_number=index, title=f"SYNTHETIC_{index:02d}",
                            docking_score=score, ligand_efficiency=-score / props["hac"],
                            parent_smiles=Chem.MolToSmiles(mol), mol=mol, pains="",
                            state_count=1, **props))
    records.sort(key=lambda r: r["docking_score"])
    for rank, record in enumerate(records, 1):
        record["score_rank"] = rank
    fps = pipeline.make_fingerprints(records, 2, 2048)
    clusters, tree = pipeline.cluster_records(records, fps, 5)
    pipeline.annotate_review_flags(records)
    eligible, stage_report = pipeline.sequential_pool(records, 16, 12)
    pipeline.make_hierarchy_plots(args.out, records, fps, clusters, tree)
    (args.out / "selection_report.json").write_text(json.dumps(stage_report, indent=2), encoding="utf-8")
    summaries = pipeline.cluster_summary(records)
    pipeline.write_records_csv(args.out / "cluster_assignments.csv", records)
    pipeline.write_summary_csv(args.out / "cluster_summary.csv", summaries)
    pipeline.make_overview_plot(args.out / "synthetic_cluster_overview.png", records, summaries)
    review = pipeline.select_candidates(eligible, 8, 2)
    pipeline.plot_selection_stages(args.out, records, stage_report, review)
    pipeline.write_records_csv(args.out / "synthetic_review.csv", review, True)
    chosen = choose_records(review, [1, 3, 4, 6])
    assign_pockets(chosen, ["DEMO_A", "DEMO_B", "DEMO_A", "DEMO_C"])
    for row in chosen:
        row["display_title"] = row["title"]
    draw_grid([r["mol"] for r in chosen], make_legends(chosen),
              args.out / "synthetic_selected4", 2, 480, 430, 1)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    scores = [r["docking_score"] for r in records]
    hac = [r["hac"] for r in records]
    le = [r["ligand_efficiency"] for r in records]
    axes[0].hist(scores, bins=8, color="#447799")
    axes[0].set(xlabel="DockingScore", ylabel="Molecules")
    axes[1].scatter(hac, scores, color="#447799")
    axes[1].set(xlabel="HAC", ylabel="DockingScore")
    axes[2].scatter(scores, le, color="#228866")
    axes[2].set(xlabel="DockingScore", ylabel="LE = -DockingScore / HAC")
    fig.suptitle("SYNTHETIC DEMO: invented scores, no docking or binding evidence")
    fig.savefig(args.out / "synthetic_score_plots.png", dpi=160)
    plt.close(fig)
    # Invented region labels exercise EF arithmetic and review joins only.
    import calculate_pocket_ef as ef
    import review_cluster_pockets as pockets
    import plot_le_analysis as le_plots
    labels = [name for name, _ in ef.POCKET_ANCHORS] + ['Other']
    ef_records = [dict(parent_id=r['parent_smiles'], source_entry=r['entry_number'],
        docking_score=r['docking_score'], pocket=labels[i % len(labels)]) for i, r in enumerate(records)]
    ef_rows = ef.calculate_ef(ef_records, labels[:4], [10, 25, 50])
    ef.write_csv(args.out/'synthetic_regional_EF.csv', ef_rows, list(ef_rows[0]))
    ef.plot_ef(args.out/'synthetic_regional_EF.png', ef_rows,
               'SYNTHETIC DEMO: invented scores and regional labels')
    sheet = pockets.make_sheet(args.out/'cluster_assignments.csv')
    for i, row in enumerate(sheet):
        row['pocket_label'] = labels[i % len(labels)]
    pockets.write_csv(args.out/'synthetic_pocket_review.csv', sheet, list(sheet[0]))
    summary = pockets.summarize(args.out/'cluster_assignments.csv', args.out/'synthetic_pocket_review.csv')
    pockets.write_csv(args.out/'synthetic_cluster_pockets.csv', summary, list(summary[0]))
    le_out = args.out/'le'
    le_out.mkdir()
    le_records, le_audit = le_plots.load_records(args.out/'cluster_assignments.csv')
    le_plots.docking_figures(le_records, le_audit, le_out)
    report = {"synthetic": True, "molecules": len(records), "clusters": len(clusters),
              "review_size": len(review), "manual_size": len(chosen),
              "selection": "score top16 -> LE top12 -> cluster size; cap2 for small demo", "rdkit": rdBase.rdkitVersion,
              "regional_EF_and_cluster_pocket_labels": "SYNTHETIC; arithmetic and join demo only",
              "warning": "Toy molecules, invented scores and pocket labels; not research results"}
    (args.out / "demo_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
