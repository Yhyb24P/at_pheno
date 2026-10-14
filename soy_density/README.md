# 大豆芯片 SNP 密度与 protein 回归

本项目研究同一 SoySNP50K 芯片母集内，增加 SNP 数量是否改善未参与拟合材料的 protein 预测，以及物理覆盖和 LD 代表性选择是否改变这一规律。只使用加性 Ridge/GBLUP；二者按一致的编码和惩罚定义视为同一模型。

本地基因型产物为 41,346 个核染色体 SNP × 20,087 个分型材料，参考组装为 Wm82.gnm2.DTC4。表型匹配队列为 17,526 个材料，其中开发集 14,020、封存集 3,506。上述规模已核对本地产物的数组头、位点表与划分表；来源及此前完整数据核验记录见只读输入目录的 `source_manifest.json`、`gt_audit.json` 和 `partition_validation.json`。不重新读取全量源表型或封存标签。

独立输入为 `../fprior/preprocessing/processed/soybase_official_v2/`，参数唯一来源为 [config.json](config.json)。完整实验要求见 [最小实现与实验验收规范](最小实现与实验验收规范.md)。fprior 原有算法、配置和产物不修改；根项目的拟南芥设计继续独立保存。

## 运行

使用根项目的 Miniconda 环境与本地 `.venv`。本项目另外需要 matplotlib、threadpoolctl 和 pytest；PyTorch 仅承担常规双精度矩阵运算，不使用神经网络。

若本地 `.venv` 尚不存在，可用 Miniconda 环境的 Python 执行 `python -m venv --system-site-packages .venv`，复用该环境的 NumPy、PyTorch 和 pytest；随后运行 `.venv/bin/python -m pip install matplotlib threadpoolctl`。

在根项目目录运行：

```bash
.venv/bin/python -m pytest -q soy_density/tests/test_numerics.py
.venv/bin/python soy_density/experiment.py profile
.venv/bin/python soy_density/experiment.py evaluate
.venv/bin/python soy_density/experiment.py report
```

`profile` 完成第 0 外折的四内折搜索及外层重拟合；检查实际成本后，`evaluate` 在配置、输入身份和相关源码一致时复用完整折，完成剩余外折并汇总。未完成折从头计算，不将部分预测视为完成；更改身份后须使用新的输出目录。`report` 只汇总十个已完成折。

开发基因型一次载入 CPU 内存，约 2.32 GB（float32）；训练计算转换为 float64。原始数组按位点 × 材料排列，模型内部按材料 × 位点排列。分块累加关系矩阵以控制 GPU 临时占用。每次训练只用实际训练行计算 QC、均值、LD 与核尺度。完整面板两种策略共用同一拟合，避免重复计算及人为产生数值差异。

## 输出

独立运行目录由配置确定，不入 Git。每折保存完整训练规模的内层分组、QC 均值与保留掩码、两种嵌套位点排序、外层面板索引、物理与 LD 统计、关系矩阵诊断、搜索结果、逐样本预测和成本。`marker_catalog.tsv` 将索引对应到原始 SNP ID 和坐标。

十折完成后生成 pooled 指标表、配对误差与组 bootstrap 区间、密度曲线 PNG/PDF、成本汇总及完成标记。`complete.json` 只在相关计算全部通过时生成。运行结果与未完成事项以实际产物为准，本文不把已实现入口等同于实验已完成。

## 已完成实验结果

完整开发集十折结果见[研究报告](results/analysis.md)、[指标表](results/metrics.tsv)和[密度曲线](results/density_performance.png)。[完整结果压缩包](results/soy_density_v1_artifacts.zip)包含全部折外预测、面板清单、训练侧统计、计算成本和对应源码快照；解压后保留原项目相对目录。数据源基因型及封存表型不在包内。文件完整性摘要见[checksums.sha256](results/checksums.sha256)。
