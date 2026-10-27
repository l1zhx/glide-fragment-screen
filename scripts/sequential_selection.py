"""Explicit score gate -> LE gate -> cluster-size review order; no weighted score."""
from collections import Counter
import math


def annotate_review_flags(records):
    for r in records:
        flags = []
        for condition, name in [(bool(r.get("pains")), "PAINS"),
                                (abs(r.get("formal_charge", 0)) > 1, "ABS_CHARGE_GT_1"),
                                (r["hac"] < 6, "HAC_LT_6"),
                                (r.get("mw", 0) > 350, "MW_GT_350"),
                                (r.get("clogp", 0) > 4, "CLOGP_GT_4")]:
            if condition:
                flags.append(name)
        r["review_flags"] = ";".join(flags)


def sequential_pool(records, score_top=2000, le_top=500, exclude_flags=False):
    if min(score_top, le_top) < 1 or not records:
        raise ValueError("Use a nonempty cohort and positive stage sizes")
    if len({r["entry_number"] for r in records}) != len(records):
        raise ValueError("entry_number must be unique")
    for r in records:
        if r["hac"] <= 0 or not all(math.isfinite(float(r[k])) for k in
                                   ("docking_score", "ligand_efficiency")):
            raise ValueError("Scores and LE must be finite; HAC must be positive")
        if not math.isclose(r["ligand_efficiency"], -r["docking_score"] / r["hac"], abs_tol=.0006):
            raise ValueError("LE must equal -DockingScore/HAC")
        r.update(pass_score_gate=0, pass_le_gate=0, pass_review_gate=0,
                 le_rank_within_score_pool="", eligible_cluster_size=0)
    by_score = sorted(records, key=lambda r: (r["docking_score"], r["entry_number"]))
    threshold_score = by_score[min(score_top, len(by_score)) - 1]["docking_score"]
    score_pool = [r for r in by_score if r["docking_score"] <= threshold_score]
    for r in score_pool:
        r["pass_score_gate"] = 1
    by_le = sorted(score_pool, key=lambda r: (-r["ligand_efficiency"],
                                            r["docking_score"], r["entry_number"]))
    threshold_le = by_le[min(le_top, len(by_le)) - 1]["ligand_efficiency"]
    for rank, r in enumerate(by_le, 1):
        r["le_rank_within_score_pool"] = rank
    le_pool = [r for r in by_le if r["ligand_efficiency"] >= threshold_le]
    for r in le_pool:
        r["pass_le_gate"] = 1
    eligible = [r for r in le_pool if not exclude_flags or not r.get("review_flags")]
    counts = Counter(r["cluster_id"] for r in eligible)
    for r in eligible:
        r["pass_review_gate"] = 1
    for r in records:
        r["eligible_cluster_size"] = counts.get(r["cluster_id"], 0)
    report = dict(input_parents=len(records), score_top_requested=score_top,
                  score_pool=len(score_pool), docking_score_cutoff=threshold_score,
                  le_top_requested=le_top, le_pool=len(le_pool), le_cutoff=threshold_le,
                  exclude_review_flags=exclude_flags, eligible=len(eligible),
                  excluded_by_review_flags=len(le_pool) - len(eligible),
                  eligible_clusters=len(counts), gate_boundary_ties="include all",
                  final_order="global cluster size descending; within cluster score ascending, LE descending",
                  weights=None)
    return eligible, report


def select_candidates(eligible, number, max_per_cluster=1):
    if min(number, max_per_cluster) < 1:
        raise ValueError("number and max_per_cluster must be positive")
    groups = {}
    for r in eligible:
        groups.setdefault(r["cluster_id"], []).append(r)
    ordered_groups = []
    for members in groups.values():
        members = sorted(members, key=lambda r: (r["docking_score"],
                                                 -r["ligand_efficiency"], r["entry_number"]))
        ordered_groups.append(members)
    ordered_groups.sort(key=lambda g: (-g[0]["cluster_size"], g[0]["docking_score"],
                                       -g[0]["ligand_efficiency"], g[0]["cluster_id"]))
    # Explicit rounds maintain broad coverage when an explicit cap >1 is used.
    result = []
    for round_index in range(max_per_cluster):
        for members in ordered_groups:
            if round_index < len(members):
                result.append(dict(members[round_index], selection_rank=len(result) + 1))
                if len(result) == number:
                    return result
    return result  # Never relax gates or per-cluster limits to fill a shortage.
