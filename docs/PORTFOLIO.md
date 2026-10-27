# 发布到l1zhx：项目与个人主页

目标账号：https://github.com/l1zhx

已读取该账号公开仓库列表及现有AF3项目说明；当前可见 `af3-ligand-screen`，没有读取到同名主页仓库。本交付准备两个独立目录：

- `glide-fragment-screen/`：研究分析项目，上传到已创建的 `l1zhx/glide-fragment-screen`。
- `l1zhx/README.md`：个人主页，上传到拟建的同名公开仓库 `l1zhx/l1zhx`。如果发布前主页已存在，先合并内容，不要覆盖现有介绍。

## 网页发布

1. 在GitHub创建公开空仓库 `glide-fragment-screen`，上传该目录**里面的内容**至仓库根目录，不要再套同名目录。包含隐藏文件 `.github/`、`.gitignore`、`.gitattributes`。
2. 提交后检查Actions中的public-demo任务。CI验证公开Python部分，不运行Schrödinger、真实docking或Slurm。
3. 创建账号同名公开仓库 `l1zhx`，将交付的主页README放在仓库根目录；主页项目链接在Glide仓库发布后生效。
4. 在个人主页的Pinned中并列固定 `af3-ligand-screen` 与 `glide-fragment-screen`。

上传源码目录，不上传ZIP本身、results/、私有CSV/SMILES、MAEGZ、轨迹、grid或本地环境。交付ZIP已按白名单打包；运行后新增结果不自动适合公开。

## 使用Git

当前Glide仓库已有初始README，应先clone，保留其初始提交。VS Code远程SSH的完整步骤见 [UPLOAD_VSCODE_SSH.md](UPLOAD_VSCODE_SSH.md)。以下假定解压文件位于release_20261027/，且克隆目录尚不存在：

```bash
git clone git@github.com:l1zhx/glide-fragment-screen.git
cd glide-fragment-screen
cp -a ../release_20261027/glide-fragment-screen/. .
sha256sum -c MANIFEST.sha256
git add README.md README.zh-CN.md LICENSE CHANGELOG.md requirements.txt VERSION MANIFEST.sha256 .gitignore .gitattributes .github scripts tests slurm config docs
git diff --cached --stat
git commit -m "Add reproducible beta-catenin Glide screening workflow"
git push -u origin main
```

个人主页仓库创建并初始化README后，单独clone；先检查已有内容，仅在只有初始化标题时使用下面的整页替换。已有介绍则合并profile_snippet.md的内容：

```bash
cd ..
git clone git@github.com:l1zhx/l1zhx.git profile_checkout
cd profile_checkout
cp ../release_20261027/l1zhx/README.md README.md
git add README.md
git commit -m "Introduce computational molecular screening projects"
git push -u origin main
```

已有仓库应clone后合并文件再提交；这些命令不使用强制推送。若只想添加Glide项目一节，可使用 `docs/profile_snippet.md`，保留原主页其他内容。

## 项目展示

About建议：

> Reproducible beta-catenin Glide SP analysis: chemical clustering, regional score enrichment and traceable candidate selection.

Topics建议：`cheminformatics`, `virtual-screening`, `molecular-docking`, `rdkit`, `beta-catenin`, `reproducible-research`。

可以描述已完成约6万分子docking、早期聚类与人工分析，并实现新版固定100类和逐级筛选。运行新版真实数据之前，不把合成demo写成新版真实Top30；没有活性数据时不写“发现抑制剂”或“验证可成药口袋”。

讲解项目时依次说明：研究问题与grid覆盖 → pose/母体身份 → 分数与LE → 化学分组与类内一致性 → 逐级门槛与多样性 → 四区域富集及背景选择 → 人工结构复核 → swimming中进一步观察什么。
