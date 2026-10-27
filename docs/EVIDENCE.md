# Evidence ledger

| Claim | Source | What this release supports |
| --- | --- | --- |
| β-catenin vs approximately60k molecules in Glide SP | Researcher account | Campaign context, not independently audited total success |
| Protein-centroid grid intended as global docking | Researcher account | Search intent; actual box coverage still needs archived grid settings |
| 4525 clusters at0.60;3729 at0.55 | Reported server results | Historical observations; private data and exact version lock not supplied |
| 30 candidates chosen after reviewing Top60 | Explicit rank list and pocket labels | Reproducible mapping to the ORIGINAL saved Top60 |
| Selected pocket labels16/11/2/1 | Direct count of supplied mapping | Annotation count only |
| Real molecular names, scores and poses for final30 | Original CSV/MAEGZ required | Not bundled; never reconstructed from guesses |
| Binding, inhibition, druggability or functional selectivity | No experiments in this package | Not established |
| Swimming / MD pocket analysis | Researcher-reported follow-up trajectories | Optional existing-volume annotation supplied; no trajectory validation or causal claim here |

The repository contains public synthetic demonstrations. Toy molecules are real chemical graphs but their scores and pocket labels are invented. They must not be described as β-catenin hits. The source scripts are adapted from the analysis code developed for this project; public-release improvements are documented in [CHANGELOG](../CHANGELOG.md).

The exact receptor structure, construct/chain numbering, preparation settings, ligand input manifest, grid dimensions and original docking job input/log were not supplied for this release. The code reproduces the postprocessing method when real inputs are provided; it does not promise exact reproduction of the original docking campaign from scratch.

The2026-10-21 update implements fixed100 groups and sequential score/LE gates. No private-data100-cluster count, new real Top30 or new pocket map has been computed locally. Historical4525/3729 counts and30 manual annotations remain historical records only.

The 2026-10-27 release adds four-region EF, chemical-cluster/pocket review summaries and the latest MED volume annotation. EF arithmetic and review joins are exercised with synthetic scores/labels. No real full-library EF values were reconstructed or published. Regions are operational residue neighborhoods, not experimentally established pocket boundaries.
