# at_pheno

本项目研究拟南芥基因组标记密度对表型预测精度的影响，比较加性基因组模型在不同单核苷酸多态性（SNP）数量下的预测表现。主要比较为 25 万个常见 SNP 与训练样本中全部符合筛选条件的常见 SNP。

基因型来自 1001 Genomes v3.1，表型来自 AraPheno。处理后的基因型矩阵包含 1135 个材料和 10,707,430 个位点；四项研究性状合计有 3020 条材料与性状对应的观测值。

## 目录

```text
preprocessing/
  README.md       数据来源、处理方法与文件格式
  raw/            原始基因型、样本信息和表型数据
  processed/      基因型矩阵、表型表、位点信息和样本划分
src/at_pheno/     数据读取与加性核模型
configs/          运行参数和研究分析规则
tests/            数值计算与训练、测试数据隔离的测试
docs/研究设计.md  研究问题、实验设计和统计方法
environment.yml   Miniconda 环境定义
```

原始数据和处理产物保存在本地，不随 Git 仓库分发。克隆仓库后，需将数据放入上述对应目录。数据来源及文件说明见 [数据预处理](preprocessing/README.md)，实验方法见 [研究设计](docs/研究设计.md)。

## 安装

在仓库根目录使用 Miniconda 创建并激活环境：

```bash
conda env create --prefix ./.venv --file environment.yml
conda activate "$PWD/.venv"
```

已有环境时直接执行激活命令。环境定义中的可编辑安装会提供 `at-pheno` 命令。

## 使用

以下示例生成小规模合成数据，并运行加性核模型：

```bash
at-pheno demo --out /tmp/at-pheno-demo
at-pheno run --data /tmp/at-pheno-demo --trait synthetic --config configs/smoke.toml --out /tmp/at-pheno-run
python -m pytest -q
```

运行真实数据时，`--data` 指向 `preprocessing/processed`，`--trait` 使用 AraPheno 性状编号。默认参数文件为 `configs/run.toml`；输出目录必须是尚不存在的新目录。

`run` 实现了训练样本内的位点筛选、缺失值填补、内层交叉验证选参和外层预测。输出包括样本划分、所选位点、参数选择结果、逐样本预测和评价指标。

该命令使用固定正则化参数网格，并自行生成内层划分；它尚未实现研究设计中的物理分箱位点排序、参数搜索范围扩展及统计区间估计。`--splits` 接受单项性状的一份外层划分表，不能直接读取包含多项性状和内层划分的合并表。现有命令适用于模型开发与计算验证，不能据此报告完整研究设计的结论。

## 开发约定

数据处理内容集中在 `preprocessing/`，模型实现集中在 `src/at_pheno/`。修改优先在现有模块中完成；同一用途只维护一份说明和一套参数来源。

检查围绕具体改动执行。已确认的数据不因目录移动而重新下载或全量校验；输入内容或相关计算方法发生变化时，检查受影响的部分。运行结果写入独立目录，避免覆盖原始数据与已有产物。
