# HDF5 方向不变预备分支

此分支使用官方 v3.1 `imputed_snps_binary.hdf5`，服务于数据可行性、密度和群体结构预备分析；不是 VCF-derived ALT-dosage 的确认性分支。

## 已核验

完整文件扫描的输入 SHA256 为 `f3713e66f017187fb3d30b8d6b854794393ab259e6b5bae4b18ca2cafcbb1862`。`snps` 为10,709,949×1135的binary矩阵；五条染色体分别为2,597,825、1,868,869、2,194,365、1,767,088、2,281,802个位点。全文件扫描未发现非0/1值；全体样本下各染色体MAC≥5的位点数分别为1,177,150、833,790、988,820、798,347、1,048,941，合计4,846,048。完整原始报告在忽略的`data/raw/hdf5_binary_scan_2026-09-22/`；运行同一脚本可重建。

这里的1只表示该文件的二态编码状态，**不解释为ALT allele**。对训练样本中心化的二态矩阵 Z，逐位点全体翻转使 Z 变为−Z，因此关系矩阵 ZZᵀ 不变。代码以独立的`BinaryQC`、`fit_binary_qc()`和`binary_additive_kernel()`实现这一约束；测试验证任意三个位点的0↔1翻转不会改变训练或测试核。

## 能回答的预备问题

- 在固定binary编码和训练内频率筛选下，密度梯度的 additive relationship 是否已饱和。
- PCA、IBS/kinship和近克隆候选的样本结构；这些输出用于冻结划分，不把它们当作表型模型的训练内特征。
- 物理位置及染色体边界的工程可行性、marker数量、读取时间和显存/CPU成本。
- 对方向不敏感的 MAF/MAC、carrier count 和 hash-density 组成审计。

## 不能回答的问题

- REF→ALT效应方向，或包含等位基因方向的回归系数。
- 原始missingness、异合原调用、calling uncertainty，或fold-specific imputation。
- REF/ALT依赖的功能注释、sequence-context variant effect和等位基因特异解释。
- 需要原始VCF语义的确认性结论。

因此，HDF5的pilot必须在产物中标记`orientation_unknown_binary`和`pilot_not_confirmatory`。当前通用CLI拒绝将此表示用于`--mode formal`；不能为绕过门禁把0/1转换成伪ALT dosage。

## 计算边界

全marker关系矩阵虽只存储约1135²个元素，但直接计算仍是O(n²M)。在CPU上对10.7M位点反复执行外/内CV不属于廉价操作。先以hash嵌套10K、50K、250K、1M的分块测试测量耗时；ALL只能在预算已实测且不破坏隔离规则时运行。预计算文件中的`kinship_ibs_mac5.hdf5`可作为外部一致性检查，不替代每个训练折重新估计的关系矩阵。

正式 VCF 分支需要单独核验样本、位点主键、REF/ALT、原始缺失和HDF5位置映射。两分支的结果只能在明确输入语义后并列，不能相互替代。
