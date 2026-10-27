# Data and software availability

The public repository distributes project-authored analysis scripts, documentation, synthetic demo graphs and invented scores. It does not contain the private compound library, receptor coordinates, docking grid, commercial software or real docking poses. Do not infer molecule identities from the historical rank list alone.

Keep original compounds, MAEGZ, grid archives and full run manifests in private storage. Manifests and logs can include local paths. Even a CSV may disclose compound identities through SMILES. Before sharing derived figures, apply the permissions governing the underlying research data.

Schrödinger must be independently installed and licensed; it is not installable through this requirements.txt. Project-authored code is MIT licensed, consistent with the related AF3 project. RDKit and other dependencies retain their own licenses. The code license grants no rights to private compounds or vendor assets.

The public-demo requirements express compatible ranges, not a frozen research environment. Store the vendor release, Python/NumPy/RDKit versions and input hashes alongside each real campaign. `clustering_metadata.json` records available Python package versions; also archive the actual installed Schrödinger release separately.

The2026-10-21 release also records SciPy. Its fixed100/sequential workflow has not been run on the private campaign in this local environment. Historical Butina results and manual pocket labels must not be presented as newly computed100-group results.
