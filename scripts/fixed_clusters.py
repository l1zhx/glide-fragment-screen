"""Exactly K chemical groups using average-linkage Tanimoto distances."""
import numpy as np
from rdkit import DataStructs


def tanimoto_condensed(fingerprints):
    """SciPy upper-triangle row order (different from RDKit Butina order)."""
    n = len(fingerprints)
    result = np.empty(n * (n - 1) // 2, dtype=np.float64)
    start = 0
    for i in range(n - 1):
        values = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[i + 1:])
        result[start:start + len(values)] = 1.0 - np.asarray(values)
        start += len(values)
    return result


def fixed_clusters(records, fingerprints, n_clusters=100):
    from scipy.cluster.hierarchy import linkage, cut_tree
    if len(records) != len(fingerprints) or not 1 <= n_clusters <= len(records):
        raise ValueError("Require 1 <= n_clusters <= number of parent molecules")
    if len(records) == 1:
        tree = np.empty((0, 4))
        labels = np.zeros(1, dtype=int)
    else:
        distances = tanimoto_condensed(fingerprints)
        tree = linkage(distances, method="average", optimal_ordering=False)
        del distances
        labels = cut_tree(tree, n_clusters=[n_clusters]).ravel()
    groups = {}
    for index, label in enumerate(labels):
        groups.setdefault(int(label), []).append(index)
    # IDs summarize size, not quality; deterministic input-order tie break.
    clusters = sorted(groups.values(), key=lambda g: (-len(g), min(g)))
    assert len(clusters) == n_clusters
    for cluster_id, members in enumerate(clusters, 1):
        member_fps = [fingerprints[i] for i in members]
        sums = np.array([sum(DataStructs.BulkTanimotoSimilarity(fp, member_fps))
                         for fp in member_fps])
        local_medoid = int(np.argmax(sums))
        medoid_index = members[local_medoid]
        sims = DataStructs.BulkTanimotoSimilarity(fingerprints[medoid_index], member_fps)
        # Mean over distinct member pairs (singleton convention = 1).
        mean_similarity = ((float(sums.sum()) - len(members)) /
                           (len(members) * (len(members) - 1)) if len(members) > 1 else 1.)
        for index, similarity in zip(members, sims):
            records[index].update(cluster_id=cluster_id, cluster_size=len(members),
                                  leaf_index=index, is_cluster_medoid=int(index == medoid_index),
                                  similarity_to_medoid=float(similarity),
                                  cluster_mean_pair_similarity=mean_similarity)
    return clusters, tree


def make_hierarchy_plots(output_dir, records, fingerprints, clusters, tree):
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    k = len(clusters)
    medoids = [next(i for i in g if records[i]["is_cluster_medoid"]) for g in clusters]
    matrix = np.array([DataStructs.BulkTanimotoSimilarity(fingerprints[i],
                       [fingerprints[j] for j in medoids]) for i in medoids])
    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap="viridis", origin="lower",
                   extent=(.5, k+.5, .5, k+.5))
    ax.set(xlabel="Cluster ID (size order)", ylabel="Cluster ID (size order)",
           title=f"{k} cluster medoids: Morgan/Tanimoto similarity")
    fig.colorbar(im, ax=ax, label="Tanimoto similarity")
    fig.savefig(output_dir / "cluster_medoid_similarity.png", dpi=180)
    plt.close(fig)
    if len(tree):
        fig, ax = plt.subplots(figsize=(18, 6), layout="constrained")
        dendrogram(tree, truncate_mode="lastp", p=k, show_leaf_counts=True, ax=ax,
                   leaf_rotation=90, leaf_font_size=7, color_threshold=0)
        ax.set(xlabel="Tree leaves/subtrees (counts in parentheses; not size-sorted cluster IDs)",
               ylabel="Average Tanimoto distance",
               title=f"Average-linkage hierarchy, truncated to {k} leaves")
        fig.savefig(output_dir / "cluster_dendrogram.png", dpi=180)
        plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(16, 9), layout="constrained")
    ids = list(range(1, k + 1))
    axes[0].bar(ids, [len(g) for g in clusters], color="#487DA6")
    axes[0].set(xlabel="Cluster ID", ylabel="Parent molecules", title=f"All {k} clusters")
    axes[1].bar(ids, [records[g[0]]["cluster_mean_pair_similarity"] for g in clusters], color="#54A58C")
    axes[1].set(xlabel="Cluster ID", ylabel="Mean pairwise Tanimoto", ylim=(0, 1),
                title="Cluster cohesion (singleton convention = 1)")
    fig.savefig(output_dir / "all_clusters_size_and_cohesion.png", dpi=180)
    plt.close(fig)


def plot_selection_stages(output_dir, records, report, selected):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
    axes[0].scatter([r["docking_score"] for r in records],
                    [r["ligand_efficiency"] for r in records], color="#BBBBBB", s=10, alpha=.35,
                    label="All parents")
    for key, color, label in [("pass_score_gate", "#4477AA", "Score gate"),
                               ("pass_le_gate", "#228866", "Score + LE gates")]:
        subset = [r for r in records if r[key]]
        axes[0].scatter([r["docking_score"] for r in subset],
                        [r["ligand_efficiency"] for r in subset], s=12, color=color, alpha=.5, label=label)
    axes[0].scatter([r["docking_score"] for r in selected],
                    [r["ligand_efficiency"] for r in selected], color="#CC5544", marker="*", s=60,
                    label=f"Selected ({len(selected)})")
    axes[0].axvline(report["docking_score_cutoff"], linestyle="--", color="#4477AA")
    axes[0].axhline(report["le_cutoff"], linestyle="--", color="#228866")
    axes[0].set(xlabel="DockingScore", ylabel="LE = -DockingScore / HAC", title="Sequential gates")
    axes[0].legend(fontsize=8)
    names = ["Input", "Score gate", "LE gate", "Eligible", "Selected"]
    counts = [report["input_parents"], report["score_pool"], report["le_pool"],
              report["eligible"], len(selected)]
    bars = axes[1].bar(names, counts, color=["#AAAAAA", "#4477AA", "#228866", "#44AA99", "#CC5544"])
    axes[1].bar_label(bars, padding=3)
    axes[1].set(ylabel="Parent molecules", ylim=(0, max(counts)*1.2), title="Actual stage counts")
    fig.savefig(output_dir / "sequential_selection.png", dpi=180)
    plt.close(fig)
