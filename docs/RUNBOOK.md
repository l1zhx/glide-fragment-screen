# 新版服务器流程：100类与逐级筛选

## 1. 环境和原始资料

归档受体/配体准备、化合物ID、原grid、Glide任务.in与日志、软件release、实际成功/失败数量。已有docking无需重跑。只有蛋白质心位置不足以证明grid覆盖全蛋白。

先验证当前Schrödinger Python依赖：

```bash
"$SCHRODINGER/run" python3 -c "import numpy, scipy, matplotlib; from rdkit import Chem; print(scipy.__version__)"
```

若供应商环境缺少SciPy，可在独立环境安装requirements.txt并使用保存的CSV重新聚类；不要盲目升级供应商Python中的NumPy/RDKit。完整MAEGZ分析仍要求相应环境兼容。少量配体试跑后再分析全量。

## 2. 完整结果分析

在项目根目录，使用尚不存在的输出目录：

```bash
"$SCHRODINGER/run" python3 scripts/run_analysis.py \
  --input /private/docking/truncated_hits..maegz \
  --out results/fixed100_sequential \
  --n-clusters 100 --score-top 2000 --le-top 500
```

总流程顺序：审计并排序输入→导出Top1万/Top100→分布图和DX→全体母体100类→分数门槛→LE门槛→簇大小与代表选择→Top100结构图→母体LE/尺寸图。主要输出位于prepared/、global/、cluster/、le/。run_manifest.json记录步骤状态与文件哈希；失败看相应日志，并使用新目录重跑。

主入口可用--number控制pose人数，用--receptor-entry明确受体条目。可选--with-ef --ef-input读取完整背景计算四区域EF，见[EF指南](REGIONAL_EF.md)。EF不会自动改用Top1万背景；若未提供--ef-input则使用主输入文件的全部可用条目。

prepare_top_hits.py只证明所给文件内的Top N；如需全库Top N，请输入完整pose文件。实际pose数、排序状态、受体条目见prepared/input_audit.json。受体识别属于启发式，必要时单独运行审计脚本并明确指定--receptor-entry。

## 3. 只重做100类与筛选

如果truncated文件已按分数排序：

```bash
"$SCHRODINGER/run" python3 scripts/cluster_glide_top1w.py \
  --input /private/docking/truncated_hits..maegz \
  --output-dir results/cluster100 \
  --n-clusters 100 --score-top 2000 --le-top 500
```

独立脚本读取前10000个有分数条目，不能把未排序全库直接当Top1万。若保存过旧cluster_assignments.csv：

```bash
python scripts/cluster_glide_top1w.py \
  --input-csv /private/old_cluster/cluster_assignments.csv \
  --output-dir results/cluster100_csv --n-clusters 100 \
  --score-top 2000 --le-top 500
```

CSV模式处理所给CSV的全部记录，重新标准化和去重，不受--max-molecules限制；输入应为希望分析的Top1万母体表。CSV中的原entry_number仍保留，但此模式不自动回写MAEGZ。

## 4. 筛选规则与输出检查

默认分数前2000（越负越好）→其中LE前500（越大越好）→簇大小排序，每簇最多1个。边界同分全保留，因此阶段人数可能略大于门槛。实际cutoff记录在selection_report.json中。所有筛选都基于化学母体记录，不是受体或重复状态。

首先打开：

- cluster_summary.csv：100类大小、medoid、类内平均Tanimoto、合格成员数量。
- cluster_assignments.csv：每个分子的类、LE与pass_score_gate/pass_le_gate/pass_review_gate。
- stage1_score_pool.csv与stage2_le_pool.csv：每一级的完整名单。
- eligible_cluster_representatives.csv：所有合格簇的候选，便于查看小簇。
- diverse_top60_for_review.csv/.maegz：默认最多60个；不足时不会补充未达标者。
- all_clusters_size_and_cohesion.png、cluster_dendrogram.png、cluster_medoid_similarity.png、cluster_medoids.svg：聚类分布及类内一致性。

cluster_size是全体母体中该类的大小；eligible_cluster_size是两级筛选及可选警示排除之后的大小。第三步按前者排序，后者供判断是否只有个别分子得分较好。大簇可能化学上很宽泛，需要查看cohesion图。

PAINS/理化警示默认只标注。若决定排除，使用--exclude-review-flags；这个选项在LE门槛之后执行，报告排除数量，不补录。

## 5. 不足60时和敏感性

可扩大门槛池，例如：

```bash
"$SCHRODINGER/run" python3 scripts/cluster_glide_top1w.py \
  --input /private/docking/truncated_hits..maegz --output-dir results/cluster100_wider \
  --n-clusters 100 --score-top 3000 --le-top 1000
```

或明确允许每簇2个：--max-per-cluster 2。程序先各类第1个，再各类第2个，不自动突破筛选门槛。若想先比较参数而不重算层次聚类：

```bash
python scripts/selection_sensitivity.py \
  --csv results/cluster100/cluster_assignments.csv --out results/gate_sensitivity \
  --score-tops 1000,2000,3000 --le-tops 300,500,1000 --number 30
```

查看sensitivity.json中的实际cutoff、候选身份、数量与overlap_jaccard.csv。对参数变化稳定和不稳定的分子分别检查，并记录最终选用门槛的原因。这一步不给出自动最优参数。

## 6. Maestro人工口袋审查与图片

新cluster ID和排名已经变化，旧30个行号不能转用。导入新的Top60，按i_user_selection_rank检查配体身份、质子化、碰撞、极性配对、埋藏/暴露程度、口袋位置与同类分子是否一致。这里不要求交叉docking。

重新选定后显式给出行号和口袋标签，例如：

```bash
"$SCHRODINGER/run" python3 scripts/export_manual_top30.py \
  --csv results/cluster100/diverse_top60_for_review.csv \
  --maegz results/cluster100/diverse_top60_for_review.maegz \
  --rows 1,3,6 --pockets ARM2_3,ARM4_5,N_terminal_bottom \
  --output-dir results/new_manual
```

真实选择30个时，--rows和--pockets都填写30项并保持一一对应。仅CSV可加--csv-only，画标准化母体而不是实际对接化学态。

只有导出**原先存档的旧Top60**时使用历史选项：

```bash
"$SCHRODINGER/run" python3 scripts/export_manual_top30.py \
  --csv /private/original_cluster/diverse_top60_for_review.csv \
  --maegz /private/original_cluster/diverse_top60_for_review.maegz \
  --use-archived-selection --output-dir results/old_manual30
```

输出PNG/SVG/CSV/MAEGZ及selection_provenance.json。图显示名称、DockingScore、LE、Cluster、Pocket；排名仅存于数据文件。

## 7. 其余图与DX

总流程保留DockingScore直方图、DockingScore–HAC、DockingScore–LE、LE–HAC和DX，score不标单位。查看density时使用原docking坐标系，阈值如0.25仅用于显示：

空间分析默认按title保留最佳pose，与聚类的标准化母体去重不同。若title不可靠，在总入口使用--density-id-property指定可信母体ID字段；如需按pose频次而非化合物频次显示，显式使用--density-keep-all-poses并报告这一口径。不得把title重复合并后的图误称为全部独立化合物分布。

```pml
load /private/receptor_in_docking_frame.pdb, receptor
load results/fixed100_sequential/global/top1w_centroid_density.dx, pose_density
isosurface pose_hotspots, pose_density, 0.25
show cartoon, receptor
```

Top100大图可独立生成：

```bash
"$SCHRODINGER/run" python3 scripts/draw_top100_2d.py \
  --input /private/docking/top100_SP_pv.maegz --output results/top100_2D.png
```

## 8. Slurm与交接

设置GLIDE_PROJECT、GLIDE_INPUT、GLIDE_OUTPUT三个绝对路径，加载Schrödinger后提交slurm/analyze.slurm；GLIDE_N_CLUSTERS、GLIDE_SCORE_TOP、GLIDE_LE_TOP默认100、2000、500。自行指定本集群partition/account。16GB是分析任务初始资源预算，尚无全量内存基准。

最终候选交接swimming时保留化合物ID、原化学态、pose、口袋、筛选理由和源文件哈希。本代码不提交新的docking或swimming。真实30分子名单和图片需在服务器运行后产生；公开demo不是研究结果。

## 9. 同类分子是否指向相同口袋

对全体母体的Maestro文件做口袋标注时，先生成匹配的复核表：

```bash
python scripts/review_cluster_pockets.py init \
  --csv results/glide100/cluster/cluster_assignments.csv \
  --out results/glide100/cluster/pocket_review.csv
```

按entry_number/名称核对原pose，填写pocket_label：ARM2_3、N_terminal_bottom、ARM4_5、C_terminal、Other或Ambiguous。未检查的留空，不能自动当作Other。keep与review_notes用于人工记录，不自动改变筛选规则。保持parent_smiles、cluster_id与source_csv_sha256不变。

```bash
python scripts/review_cluster_pockets.py summarize \
  --csv results/glide100/cluster/cluster_assignments.csv \
  --review results/glide100/cluster/pocket_review.csv \
  --out results/glide100/cluster_pockets
```

输出各簇已检查比例、主标签占比、标签熵和热图。未检查簇留空。主标签占比只以已检查成员为分母，因此必须同时看review_coverage。若只想检查Top60，将两个命令的--csv均换成diverse_top60_for_review.csv；这只代表候选子集，不能推断整个化学簇的空间一致性。

## 10. 独立LE图与可选MD

```bash
python scripts/plot_le_analysis.py \
  --csv results/glide100/cluster/cluster_assignments.csv --out results/le_extra
```

输出Score–HAC等LE线、LE–HAC、分口袋/分HAC描述图及Score–HAC Pareto候选。若使用带pocket列的CSV则按区域着色；聚类CSV本身没有自动口袋标签。默认不做理论亲和力换算。

已有口袋体积的MED标注见[swimming接口](SWIMMING_INTERFACE.md)。该模块不改MDpocket参数、不重跑体积、不显示static reference。

可将config/analysis.example.sh复制为被git忽略的analysis.local.sh，编辑输入/输出路径后执行。模板默认只分析已有docking；EF需明确打开。
