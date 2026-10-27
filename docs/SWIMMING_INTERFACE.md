# 从docking候选到swimming：交接和体积标注

此模块是可选后续步骤。Glide项目不生成或提交MD轨迹，也不把docking分数、pose density或EF当作动态结合验证。当前保留最初MDpocket参数及输出，**不为强制连续识别口袋而改参数**。

## 候选交接

最终30个分子应保留：原化合物ID、对接化学态、母体SMILES、最佳pose原条目、DockingScore、HAC及其定义、LE_score、cluster ID/size、口袋、人工审查理由、源文件哈希。原化学态与标准化母体可能不同；swimming体系准备需要重新核对质子化、立体化学与参数化。

新的100类和Top60排序不能复用旧Top60行号。口袋标签以当前四区域为准；旧Ser246/Thr150等标注仅在原存档中保留。候选可能在动态过程中访问其他区域，这不是自动判定原docking“对/错”的充分依据。

## 输入已有ARM4–5体积和MED状态

`scripts/plot_volume_occupancy.py` 是最新标注逻辑的独立公开入口，仅读取已有结果，不调用MDpocket：

```bash
python scripts/plot_volume_occupancy.py \
  --volume /private/pocket_growth_r1/ARM4_5/volume_timeseries.csv \
  --occupancy /private/site_analysis/stable_occupancy_by_frame.csv \
  --descriptors /private/pocket_growth_r1/ARM4_5/measure/mdpout_descriptors.txt \
  --out results/ARM45_MED_original
```

体积CSV必须只有MD采样行，含 `time_ns,volume_A3`；如有 `nb_AS` 则直接使用。若缺少nb_AS，需原描述符文件，其第一行数据是static reference，之后逐行对应MD，数量和体积必须完全匹配。原描述符保持不变，不显示static-reference线。

状态CSV含严格递增的 `time_ns` 和默认 `MED_ARM4-5` 列，值为0/1或true/false。可用 `--state-column` 指定另一列。两表必须来自同一replica和时间原点；默认最近时间差不超过0.001 ns。不匹配会报错，不插值，不把未知状态当作未结合。

标注应使用逐帧接触状态，而不是手工把“前1 μs”全部标为结合。已有stable状态依赖其原接触阈值、持续时间与gap规则；报告这些定义。若只有MED列，其他fragment是否占据记unknown。只有18种fragment的对应列均存在时，才进一步区分MED alone、MED with others等。

## 输出与解释

- `volume_with_MED_occupancy.png/svg`：原始体积随时间、MED阴影区间、占据条带；零值保留灰色标记，缺口不插值。
- `MED_present_absent_comparison.png/svg`：检测比例与检测到且体积为正时的条件分布。
- `volume_occupancy_joined.csv`、`conditional_volume_summary.csv`、`MED_intervals.csv`：对齐后的数据、计数与区间。
- `MED_present_minus_absent.json`：描述性效应差，未计算独立帧p值。
- `analysis_provenance.json`：来源哈希、时间容差、状态列及局限。

MDpocket的0可能表示没有匹配到通过算法筛选的alpha spheres，不足以证明物理口袋消失。应同时报告检测比例、检测阳性条件下的体积、原始含0均值；三者回答不同问题。只画正体积箱线图会隐去检测差异。

同一轨迹帧相关，不能把数百帧当数百个独立重复做普通t检验。早期结合、晚期不结合还混杂时间漂移；分块bootstrap也不能从一次解离事件中消除这个混杂。多个独立replica与多次结合/解离事件才更适合估计不确定性。当前图只能表述条件分布/关联，不能推出“MED离开导致口袋变小”。MED absent的混合fragment体系也不是apo对照。

统计背景参见 [Grossfield等：模拟数据的不确定性和采样质量](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286151/)。
