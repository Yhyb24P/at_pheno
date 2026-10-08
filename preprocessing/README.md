# 数据预处理

本目录保存基因型与表型原始数据、处理产物及其格式说明。`raw/` 按基因型和表型分为两类；`processed/` 保存可供模型读取的数据文件。

## 数据来源

| 数据 | 来源 | 本地位置 |
| --- | --- | --- |
| 基因型 VCF | [1001 Genomes v3.1](https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz) | `raw/genotypes/` |
| 预插补 SNP 矩阵与亲缘矩阵 | [1001 Genomes HDF5 数据包](https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz) | `raw/genotypes/` |
| 1135 个材料的样本信息 | [1001 Genomes 样本表](https://1001genomes.org/json/1135_accessions_table.json) | `raw/genotypes/1135_accessions_table.json` |
| 群体标签 | [ADMIXTURE 数据表](https://docs.google.com/spreadsheets/d/1_jp6KKfUC0z1WteS9HUIp4f3ssFZlCV3KXPIvrTW0Fk/export?format=csv) | `raw/genotypes/admixture.csv` |
| 性状描述与观测值 | [AraPheno](https://arapheno.1001genomes.org/) | `raw/phenotypes/` |

AraPheno 原始文件以 `phenotype_<编号>.json` 保存性状说明，以 `values_<编号>.json` 保存观测值。原始目录中保留下载记录和发布者校验信息，便于确定所使用的数据来源。

## 基因型处理

参考基因组为 TAIR10。由 VCF 保留核染色体 1–5 上的双等位 SNP，参考和替代碱基均限于 A、C、G、T；仅接受 `FILTER` 为 `PASS` 或 `.` 的记录。多等位记录整体排除，不拆分为多个位点。

矩阵按“材料 × 位点”排列，形状为 `1135 × 10707430`，使用 `int8` 存储。基因型按 VCF 替代等位基因的拷贝数编码：

| 数值 | 含义 |
| --- | --- |
| 0 | 两个参考等位基因 |
| 1 | 一个参考、一个替代等位基因 |
| 2 | 两个替代等位基因 |
| -1 | 缺失基因型 |

矩阵行与 `samples.tsv` 一一对应，列与 `variants.tsv` 一一对应。位点按染色体、位置、参考碱基、替代碱基排序；位置采用从 1 开始的坐标。`variants_compact.npy` 保存位点信息的紧凑数组表示。

原始 HDF5 是预插补的二态矩阵，用于近似相同材料的识别；模型的基因型输入使用 VCF 派生矩阵，以保留替代等位基因方向和原始缺失信息。

## 表型处理

以材料编号 `accession_id` 将表型与基因型对应，使用 AraPheno 发布的材料级数值，不在本项目中重新汇总重复测量或进行数值变换。每项性状中，每个材料对应一条记录；不同性状可共享同一材料。

| AraPheno 编号 | 发布名称 | 研究用途 | 对应材料数 |
| --- | --- | --- | --- |
| 703 | DTF1 | 主分析 | 936 |
| 748 | stomata density | 独立研究数据复核 | 330 |
| 704 | RL | 同一研究的补充分析 | 850 |
| 705 | CL | 同一研究的补充分析 | 904 |

这四项性状共形成 3020 条记录。其他候选性状的原始文件保留在 `raw/phenotypes/`，不纳入当前分析表。编号 390、261、262 的重复观测对应关系尚未确定，不将其自动平均为材料级数值。

## 处理产物

以下文件均位于 `processed/`：

| 文件 | 内容 |
| --- | --- |
| `genotypes.npy` | 基因型矩阵，约 12.2 GB |
| `samples.tsv` | 矩阵行对应的材料编号与群体标签 |
| `variants.tsv` | 矩阵列对应的染色体、位置、参考及替代碱基 |
| `variants_compact.npy` | 位点信息的紧凑数组 |
| `phenotypes.tsv` | `accession_id`、`trait_id`、`value` 三列的表型长表 |
| `traits.tsv` | 性状名称、研究用途、材料数与原始文件路径 |
| `blocks.tsv` | 材料与近似相同材料组的对应关系 |
| `splits.tsv` | 各性状、划分方案和内外层交叉验证的样本分配，共 36,240 条记录 |
| `logo_splits.tsv` | 每次留出一个群体的样本分配 |
| `marker_rank_at_pheno_common_v1*.npy` | 主排序与两套备选位点排序，元素为矩阵列索引 |
| `dataset.json` | 矩阵形状、编码、来源、性状和文件摘要 |

`splits.tsv` 中，`protocol` 区分分组平衡划分（`blocked`）和随机划分（`iid`）；`salt` 区分主划分（`canonical`）与备选划分（`s2`、`s3`）；`level` 区分外层（`outer`）与内层（`inner`）。内层记录的 `outer_fold` 指明其所属的外层测试折；`fold` 为该层的折编号。选择一份划分时须同时指定性状、方案、划分标识和层级，内层还须指定外层测试折。

位点排序按 1 Mb 物理区间分箱，箱内采用固定散列排序，再依次从各箱取位点。备选排序使用不同的固定散列参数，以评估位点选择的敏感性；排序不使用表型值。

## 训练时的数据处理

处理产物保留原始基因型缺失值。检出率、次要等位基因频率（MAF）和填补均值在每次实际训练的样本集合内估计，包括内层交叉验证的训练部分。

研究设计要求检出率不低于 0.95、MAF 不低于 0.05，并排除单态位点。缺失值用该训练集合的位点均值填补。测试样本和验证样本不参与这些参数的估计，也不预先用全体材料计算一套统一筛选结果。相应计算由 `src/at_pheno/core.py` 中的 `fit_qc` 和 `additive_kernel` 实现。

原始文件和大型数组不纳入 Git。读取基因型时可使用 NumPy 的只读内存映射；不应将整个矩阵一次转换成浮点数组。
