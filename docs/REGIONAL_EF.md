# 四区域的高分富集：regional EF

本分析回答“某区域在全局高分分子中是否占比更高”，不是根据已知活性/诱饵计算的传统验证EF，也不是经过实验校准的druggability评分。默认不做理化性质校正。

## 1. 输入与背景

尽可能输入**完整docking pose文件**。程序扫描全部有效条目，标准化化学母体，每个母体仅保留全局最低DockingScore姿态。该姿态决定该分子的区域归属，不统计它在各口袋的所有替代姿态。

如果只输入 `truncated_hits..maegz`，输出仅代表这个截断群体的高分富集。不能用Top1万的口袋比例代替全库背景比例，也不能把同一个Top1万同时称为全库Top1%。MAEGZ里的条目数、状态数和母体数不同。

## 2. 四区域

| 程序名 | 锚定残基 |
| --- | --- |
| ARM2_3 | Met202、Thr205、Lys242 |
| N_terminal_bottom | Met189、Glu226、Lys233 |
| ARM4_5 | Leu286、Asn326 |
| C_terminal | Leu644、Ser681 |

默认：各组Cα坐标的几何中心为球心，半径8 Å；配体重原子质心落在球内则为候选归属。若落入多个球且最近与次近中心距离差≤2 Å，记为Ambiguous；所有球外记Other。已有标签缺失记Unassigned。残基类型、编号或链缺失/不唯一时直接报错，不省略缺失锚点。

这些球是探索性区域，不是由SiteMap/fpocket确定的口袋边界。两个相邻区域可能重叠。必须在**原docking坐标系**查看 `regions_used.json` 和分子归属；对齐或重新平移后的受体不能与原pose混用。

## 3. 计算公式

设N为全背景母体数，n_p为其中归于区域p的母体数；K为全局高分子集实际大小，k_p为其中归于p的母体数：

`EF_p = (k_p / K) / (n_p / N)`

- 默认分析全局前1%、2%、5%、10%，分数越低越好。
- 先向上取整人数，再把边界相同分数全部纳入，因此实际K可能大于请求人数。
- N与K均包括Other、Ambiguous、Unassigned，不能只保留四口袋后重新归一化。
- n_p=0时EF为空，不能记为0；k_p=0且n_p>0时EF=0。
- EF>1表示相对所给背景富集，EF≈1表示接近背景，EF<1表示相对缺失；小n_p可能产生不稳定的高EF。

至少同时报告N、K、n_p、k_p、实际DockingScore cutoff与EF。不要只展示EF最大者。

## 4. 运行

```bash
"$SCHRODINGER/run" python3 scripts/calculate_pocket_ef.py \
  --input /private/docking/full_poses_pv.maegz \
  --output-dir results/ef_full --chain A \
  --cohort-label 'complete supplied SP pose cohort'
```

链A只是示例，须先确认。受体不在文件里时增加 `--receptor /private/prepared_receptor.mae`。程序不会自动对齐受体。输入文件保持不变。

若球心需要人工复核，用JSON列表提供四个对象，每个包含 `name`、`center:[x,y,z]`、`radius`，然后使用 `--definitions /private/reviewed_regions.json`。不要把分子质心均值当作已经验证的口袋中心。

若已在Maestro完整标注每个pose的区域，使用：

```bash
"$SCHRODINGER/run" python3 scripts/calculate_pocket_ef.py \
  --input /private/annotated_full_poses.maegz \
  --label-property s_user_pocket --output-dir results/ef_labels
```

标签使用上表程序名或Other/Ambiguous/Unassigned，不能混用Ser246、Asn326等旧命名。仅标注Top60再计算不是全库EF。自定义原化合物ID可用 `--id-property`；缺失ID会记录为拒绝条目，不自动用title补齐。

## 5. 输出检查与展示

- `audit.json`：输入哈希、有效条目、母体数、拒绝数、背景定义、各区域数。
- `rejected_entries.csv`：转换失败与无效分数等条目；这些不在可分析背景中，须报告。
- `molecule_pocket_assignments.csv`：母体、最佳pose源条目、名称、分数、质心与区域距离。
- `pocket_EF.csv/png`：四区域多个高分截点的计数、cutoff与EF曲线。
- `regions_used.json`：实际球心、半径、锚点及链信息。

汇报建议：一张受体表面/区域示意图 + 一张四区域EF曲线 + 计数表，说明输入是全库还是截断样本。总流程的 `le/` 使用聚类母体HAC；EF表的HAC来自原pose重原子，含盐/多组分时两者可能不同，不能无检查混合。

EF受grid覆盖、口袋暴露/大小、姿态采样、库组成及打分偏差影响。可以提出“区域p表现出高分姿态富集，值得后续研究”，不能仅凭EF宣称“区域p可成药”或“小分子具有抑制活性”。
