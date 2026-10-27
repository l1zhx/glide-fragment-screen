# β-catenin Glide SP：从约6万分子筛选到可追溯候选

该项目整理已有 Glide SP docking 的后处理：**输入审计 → Top1万空间/分数分析 → 固定100类 → DockingScore、LE逐级筛选 → Maestro人工复核 → 最终结构图与swimming交接**。另提供四个残基区域的 EF 分析和既有 MDpocket 体积的 MED 占据标注。

已有研究流程与新版软件实现需要区分：研究者已完成 docking、早期聚类及人工筛选；本公开版本没有私有全库文件，不能声称已重新计算出新版真实 Top30 或全库 EF。展示图均为明确标识的合成演示。

## 当前采用的规则

| 阶段 | 当前实现 |
| --- | --- |
| 数据审计 | 扫描所给文件全部有限分数，稳定排序，导出最多10000 poses和Top100 |
| 化学身份 | 标准化母体、去重复、保留母体最低分状态；原始对接态保留用于结构导出 |
| 聚类 | Morgan半径2、2048位、Tanimoto距离、平均连接层次聚类，切成恰好100类 |
| 第一级 | 默认取DockingScore最低的2000个母体，保留边界同分 |
| 第二级 | 在第一级内取LE_score最高的500个，保留边界同值 |
| 第三级 | 合格簇按整个聚类群体的cluster size排序，簇内取最低分代表，默认每簇1个 |
| 输出 | 最多60个人工复核、最多30个自动候选；不足时报告，不自动突破门槛 |
| 人工复核 | 检查化学态、碰撞、极性配对、暴露情况、口袋与类内空间一致性 |

**没有60/30/10权重，也不要求交叉docking。** 2000/500是可调审查预算，不是活性阈值；大簇也不等于好结合。固定100类不保证类内相似度均超过0.60。LE_score与DockingScore存在数学耦合，不能认为用LE彻底校正了尺寸偏差。

## 服务器主命令

在项目根目录运行，输出目录使用新名称：

```bash
"$SCHRODINGER/run" python3 scripts/run_analysis.py \
  --input /private/docking/truncated_hits..maegz \
  --out results/glide100 \
  --n-clusters 100 --score-top 2000 --le-top 500
```

路径替换为实际服务器路径。文件名不能证明输入恰为全库前1万个；查看 `prepared/input_audit.json`。如果完整pose文件可用，主流程也可直接输入完整文件，再由程序排序截取。全流程需要兼容的Schrödinger Python、RDKit、NumPy、SciPy、Matplotlib；不要为安装依赖直接升级供应商环境。

四区域EF可单独计算，也可在总命令增加 `--with-ef --ef-input /private/docking/full_poses_pv.maegz`。没有完整结果时可以算截断样本内EF，但不能将其解释为全库EF。

## 文件和阅读顺序

1. [RUNBOOK.md](docs/RUNBOOK.md)：环境、完整运行、只重聚类、阈值敏感性、人工导出、Slurm。
2. [REGIONAL_EF.md](docs/REGIONAL_EF.md)：四口袋定义、EF分母、完整/截断背景与结果解释。
3. [METHODS.md](docs/METHODS.md)：聚类、LE、pose density与筛选的具体实现及局限。
4. [SWIMMING_INTERFACE.md](docs/SWIMMING_INTERFACE.md)：最终分子交接、MED标注图、零体积及统计解释。
5. [VALIDATION.md](docs/VALIDATION.md)：实际测试范围与未完成的服务器验证。
6. [PORTFOLIO.md](docs/PORTFOLIO.md)：发布到l1zhx账号与个人主页的方法。

主要结果分布于 `prepared/`、`global/`、`cluster/`、`le/`、可选的 `regional_ef/`。DockingScore图不标单位。最终人工结构大图保留分子名称、DockingScore、LE、Cluster、Pocket，不显示旧Top60排名或Selected编号。

## 无商业软件的演示

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/demo.py --out demo_output
```

演示仅20个化学图、5类，分数与口袋标签为人为构造。另有125个受控指纹切成100类的回归检查。真实Top1万的内存、运行时间和MAEGZ接口仍需在服务器验证。

Docking优选说明计算候选值得复核；区域EF说明高分姿态的区域富集。两者均不能直接证明抑制活性、亲和力或口袋可成药性。
