# 外部实现参考

`BiMSGP/` 是独立 Git checkout，不纳入本项目源码，也未修改上游代码。

- 上游：https://github.com/liuxumail-create/BiMSGP
- 固定提交：`71084c8fa5444fbc747b0291ace627b1f2d6ef33`
- 获取日期：2026-09-22
- 用途：课题组已有方法的架构、训练协议和可复现性参照。

重建参考目录：

```bash
git clone https://github.com/liuxumail-create/BiMSGP reference/BiMSGP
git -C reference/BiMSGP checkout --detach 71084c8fa5444fbc747b0291ace627b1f2d6ef33
```

## 本次代码核验

| 位置（固定提交） | 观察 | 对新项目的影响 |
| --- | --- | --- |
| `bimsgp/model.py:82` | wide head 第一层为 `Linear(fused_dim*n_markers, 4*d_model)` | 不能直接扩展到百万级 SNP；预测头必须重新设计 |
| `bimsgp/model.py:133` 附近 | 融合特征 reshape 后进入 decoder | 这是位点数相关的 flatten head，不是固定维 pooling |
| `bimsgp/train.py:290` 附近 | 默认 selection set = outer test；`val_frac>0` 才从训练集划分验证集 | 论文复现与本研究公平比较应分开命名；正式比较强制训练内选择 |
| `bimsgp/train.py:309` | scaler 在真正用于拟合的 `X_tr` 上拟合 | 开启 val_frac 后该处预处理隔离是正确的，不应笼统说整个仓库存在全数据标准化 |
| `bimsgp/train.py:379`、`:386` | scheduler 和 checkpoint 使用 selection set | 默认路径有测试集驱动选择；不能只修 early stopping 而保留 scheduler 路径 |
| README 的 Data / protocol | 明示源矩阵列序、非 nested CV 和原始数据构建的范围 | 原文已披露的限制不等于隐瞒或结果无效；限制的是可作何种泛化结论 |

`--val-frac 0.1` 提供训练内 holdout，可以隔离 epoch selection；它不是完整内层 K-fold 超参数搜索，也没有自动保证分群、亲缘分块、原始 QC 和本项目其他规则。这里只静态检查源码，未安装 Mamba 或复跑论文实验。

以 wide、双向、d=64 计算，第一层有 `32768*M + 256` 个参数。M=33,709 时约 1.105B；M=10,707,430 时约 350.86B，仅 FP16 权重约 701.7 GB（十进制），不含优化器、梯度和激活。该公式是代码推导，不是本机实测。

公开预印本：https://www.preprints.org/manuscript/202609.1051 。其页面链接正式版 https://www.mdpi.com/2073-4425/17/9/1156 ，本次正式版正文请求受限，因此版本审计采用上述固定代码和可访问预印本，不以搜索未命中否认正式发表状态。
