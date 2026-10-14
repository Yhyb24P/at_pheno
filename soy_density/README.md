# 植物全基因组 SNP 连续表型回归：BiMSGP 对标

本项目以 BiMSGP 为主要对标方法，使用大豆 SoyBase 数据研究全基因组 SNP 的连续表型回归预测。核心问题是相同数据下 BiMSGP 与 GBLUP/Ridge 的预测差异；SNP 密度分析用于描述预测规律，计算成本和高维表征的不足为后续算法研究提供依据。目标是评价方法，不以重现论文原始分数为完成条件。当前先使用 protein。

## 两套数据与比较范围

| 用途 | 数据与材料 | SNP | 划分与参数 |
| --- | --- | --- | --- |
| 当前 BiMSGP 历史成绩比较 | 本机已拉取的 `fprior/data/soybase/geno-big.parquet`、`pheno_selected.parquet`；protein 15,878 个材料 | 原列序 42,195 个，包含 scaffold | 沿用既有 protein 十折文件；[benchmark.json](benchmark.json) |
| 已完成加性密度分析 | `fprior/preprocessing/processed/soybase_official_v2/`；14,020 个开发材料，3,506 个封存材料 | 官方 VCF 的 41,346 个核染色体 SNP，Wm82.gnm2.DTC4 | 已确定近重复组十折；[config.json](config.json) |

本次按用户要求不重训 BiMSGP、不运行其推理，复用已有论文逐折指标，在旧数据和旧外折上重新计算加性基准。BiMSGP 原成绩按测试折选 epoch，本轮基准仅在内层选参；因此只能报告有明确协议限制的历史成绩比较，不能称为统一无测试泄漏协议下的完整方法对标。旧数据 protein 十折 checkpoint 已有，逐材料 BiMSGP 预测尚未找到。

两套矩阵与标签不可混用：旧矩阵覆盖当前开发集中的 12,576 个材料；按染色体与位置对应有 41,306 个共同 SNP，当前目录另有 40 个位点不在旧矩阵中。旧队列中有 3,139 个材料属于官方队列的封存名单。本次旧队列运行只读取指定旧表型，不读取官方封存标签文件；这项比较不作为官方队列的独立封存验证。

实验要求见 [最小实现与实验验收规范](最小实现与实验验收规范.md)。fprior 的输入与原实验保持不变，拟南芥研究独立保存。

## 运行

使用根项目的 Miniconda 环境与本地 `.venv`。加性计算需要 NumPy、PyTorch、threadpoolctl，绘图使用 matplotlib，旧数据读取使用 PyArrow。当前 `benchmark.py` 复用本地 `fprior/.deps` 的 Parquet 依赖，不安装或加载 Mamba。

若本地 `.venv` 尚不存在，可用 Miniconda 环境的 Python 执行 `python -m venv --system-site-packages .venv`，复用该环境的 NumPy、PyTorch 和 pytest；随后运行 `.venv/bin/python -m pip install matplotlib threadpoolctl`。

在根项目目录运行：

```bash
.venv/bin/python -m pytest -q soy_density/tests/test_numerics.py
.venv/bin/python soy_density/benchmark.py evaluate
```

`benchmark.py` 一次将旧基因型载入内存，沿用原十折；每个外折以四个内折选择正则化参数，先测量第一个完整外折，再继续剩余折。两个加性流程分别为训练均值中心化的 `GBLUP-centered` 和按训练 SNP 方差缩放的 `Ridge-standardized`，属于同一模型类的两种编码，不计作独立算法。所有 42,195 个原标记列保留，不套用官方 VCF 的 MAF/QC 规则。输出为配置中的独立目录。

此前的官方 VCF 密度入口保留：

```bash
.venv/bin/python soy_density/experiment.py profile
.venv/bin/python soy_density/experiment.py evaluate
.venv/bin/python soy_density/experiment.py report
```

`profile` 完成第 0 外折的四内折搜索及外层重拟合；检查实际成本后，`evaluate` 在配置、输入身份和相关源码一致时复用完整折，完成剩余外折并汇总。未完成折从头计算，不将部分预测视为完成；更改身份后须使用新的输出目录。`report` 只汇总十个已完成折。

官方 VCF 密度实验的开发基因型一次载入 CPU 内存，约 2.32 GB（float32）；训练计算转换为 float64。原始数组按位点 × 材料排列，模型内部按材料 × 位点排列。分块累加关系矩阵以控制 GPU 临时占用。每次训练只用实际训练行计算 QC、均值、LD 与核尺度。完整面板两种策略共用同一拟合，避免重复计算及人为产生数值差异。

## 输出

旧队列 protein 两种加性编码的十折已完成，结果见 [报告](results/legacy_protein_v1/analysis.md)、[pooled 指标](results/legacy_protein_v1/pooled_metrics.tsv)与[历史逐折对照](results/legacy_protein_v1/historical_comparisons.json)。上述结果已公开；不将不存在的 BiMSGP 新训练或逐材料预测写作已完成。

旧队列比较保存本轮逐材料预测、pooled PCC/MSE/RMSE/MAE/R²、逐折指标、选参记录和成本，并另存既有 BiMSGP、GBLUP 的历史逐折指标对照。历史 RMSE 可按折人数重构近似 pooled MSE，历史 pooled PCC 和逐材料配对误差不能从汇总指标恢复。

官方 VCF 密度实验每折保存完整训练规模的内层分组、QC 均值与保留掩码、两种嵌套位点排序、外层面板索引、物理与 LD 统计、关系矩阵诊断、搜索结果、逐样本预测和成本。`marker_catalog.tsv` 将索引对应到原始 SNP ID 和坐标。

十折完成后生成 pooled 指标表、配对误差与组 bootstrap 区间、密度曲线 PNG/PDF、成本汇总及完成标记。`complete.json` 只在相关计算全部通过时生成。运行结果与未完成事项以实际产物为准，本文不把已实现入口等同于实验已完成。

## 已发布结果

- [BiMSGP 原数据加性基准报告](results/legacy_protein_v1/analysis.md)：15,878 个材料、42,195 SNP 的本轮十折结果与历史成绩对照。[完整结果包](results/legacy_protein_v1/soybase_legacy_protein_v1_artifacts.zip)包含逐材料预测、原折分配、预处理统计、源码快照及对应历史指标；[SHA256 摘要](results/legacy_protein_v1/checksums.sha256)用于核对下载文件。
- [官方 VCF 加性密度报告](results/analysis.md)：此前 14,020 个开发材料、41,346 核 SNP 的五档密度结果。[完整结果包](results/soy_density_v1_artifacts.zip)继续保留。

结果包不包含源基因型矩阵、BiMSGP checkpoint 或官方封存标签。两项实验的材料、标签、面板和评价协议分别说明，不能将两套成绩拼作同一实验。
