# Public packaging changes

## 2026-10-27: complete docking analysis and GitHub profile integration

- Retained exact100 clustering and sequential score → LE → cluster-size ordering; no composite weights or required cross-docking.
- Integrated latest four-region EF with complete supplied background and explicit Other/Ambiguous/Unassigned cases.
- Added score/size/LE/Pareto figures, review sheets, cluster/pocket coverage summaries and heatmaps.
- Extended orchestration with configurable pose count, receptor-entry override, optional independent EF input and input hash checks.
- Added standalone annotation of original ARM4–5 volume/occupancy outputs, preserving zeros and omitting the static reference line; no parameter forcing or independent-frame p values.
- Updated public synthetic demo, regression tests, Chinese runbook and l1zhx profile files. No private campaign data or real new shortlist is bundled.
- Added an MIT license for project-authored code, matching the related public AF3 repository.

## 2026-10-21: fixed100 and sequential selection

- Replaced similarity-threshold Butina clustering with average-linkage Tanimoto clustering cut to exactly100 groups; impossible counts raise an error.
- Removed the weighted score and its percentile module from the current workflow. Selection now uses an explicit score gate, then an LE gate, then cluster-size ordering and per-cluster representatives.
- Boundary ties are included. Alerts annotate by default; exclusion is explicit. Gate/cap shortages are reported without automatic relaxation.
- Added medoids, cohesion metrics, dendrogram/heatmap/all-cluster plots, stage tables, gate visualization and a sensitivity comparison tool.
- Added saved-parent-CSV input for environments without Schrödinger. CSV mode recomputes descriptors/LE and does not export MAEGZ.
- Historical manual30 requires explicit archival mode. New ranks/pockets must be supplied after review.
- Updated public demo, tests and documentation. Original external working scripts and the2026-10-16 ZIP remain available as historical artifacts.

## 2026-10-16

- Organized the existing analysis, clustering, Top100 depiction and manual30 scripts into a standalone portable project; source working scripts outside this repository remain unchanged.
- Replaced personal server paths with relative/private-path examples.
- Added explicit pose-file audit, score sorting, stable source-entry traceability, orchestrated logging and file hashes.
- Kept historical percentile ranking as `legacy`; added opt-in average-tie ranking to the cluster script and made it explicit/default for the new orchestration flow.
- Preserved the historical30 selections and pocket labels; no automatic remapping to a newly generated Top60.
- Expanded density bounds to cover ligand extents and corrected the NumPy Gaussian fallback for grids shorter than the kernel.
- Added finite-score handling, a synthetic RDKit example, regression tests, CI, server documentation and GitHub profile text.
- Did not reconstruct missing original docking inputs, rerun private docking, or assert binding/activity validation.
