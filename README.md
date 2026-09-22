# at_pheno

拟南芥自然群体的超高维基因组预测：研究标记密度、低频变异与基因组结构何时带来可泛化的表型预测收益。

2026-09-22重新初始化。当前交付为**研究审计、修订工作协议、可运行CPU基线，以及官方HDF5/七项表型的文件级样本审计**；尚未开始真实全基因组模型比较，不含生物学结论。原稿已原样保留在 [archive](docs/archive/2026-09-22/)。

## 从这里阅读

- [研究审计与重启决策](docs/研究审计与重启决策.md)：原方案的问题、证据与修改理由。
- [修订综述](docs/拟南芥超高维基因组预测研究综述.md)、[实验协议v0.1](docs/拟南芥超高维基因组预测_实验设计方案.md)：研究主线与执行边界。
- [文献证据表](docs/文献证据表.md)：已核验、待核验和来源链接。
- [数据契约](docs/数据契约.md)、[实施路线](docs/实施路线与验收.md)：数据准备与阶段验收。
- [BiMSGP参考](reference/README.md)：课题组开源仓库固定版本与代码审计。
- [初始化验证与数据初审](docs/初始化验证与数据初审.md)：历史初审、真实表型计数及Kas-2重复记录。
- [方向不变 HDF5 分支](docs/方向不变HDF5分支.md)：binary HDF5 能做与不能做的预备分析边界。
- [正式实验门禁](docs/正式实验门禁.md)：全库表型可行性审计、近亲块标定和确认性统计的未完成项。

## 安装与自检

使用Python≥3.11与uv；本机裸`python`指向无法启动的旧环境，因此使用`uv run`。

```bash
uv sync --locked --extra dev
uv run pytest -q
uv run at-pheno --help
```

合成数据端到端验证（输出目录必须不存在；已运行时换一个run_id）：

```bash
uv run at-pheno demo --out data/interim/demo
uv run at-pheno audit --data data/interim/demo --trait synthetic
uv run at-pheno run --data data/interim/demo --trait synthetic --config configs/smoke.toml --out runs/demo_iid
uv run at-pheno run --data data/interim/demo --trait synthetic --config configs/smoke.toml --protocol group --out runs/demo_group
```

真实数据按契约准备后，可用`configs/pilot.toml`运行。当前实现均值基线、分块additive kernel ridge、每个内外训练折独立QC/填补、哈希嵌套密度和内层选参。输出原尺度OOF、划分、QC、调参及输入/代码校验值。它是pilot，不等于已实施完整实验协议。

## 小型公开数据快照

```bash
uv run python scripts/snapshot_arapheno.py --out data/raw/arapheno_snapshot --phenotype 261 262
```

工具下载表型目录及候选原始记录并写URL/时间/SHA256；不下载基因组大文件、不自动汇总重复或选正式性状。初始快照保存在`data/raw/arapheno_2026-09-22/`，大文件和原始数据不入Git。

HDF5审计需额外依赖：`uv sync --locked --extra dev --extra data`。真实数据准备进展与来源见 [样本交集审计](docs/真实数据准备与候选面板.md)。

## 目录

```text
configs/              pilot和软件自检配置
src/at_pheno/         数据审计、训练内QC、分块kernel、嵌套评估
tests/                泄漏隔离与数值正确性测试
scripts/              小型公开数据快照
docs/                 工作稿、证据、原稿归档、数据/运行审计
data/raw/             原始下载（忽略）
data/interim/         中间数据与合成测试（忽略）
data/processed/       按契约处理的真实输入（忽略）
data/manifests/       可入库的小型来源记录
reference/BiMSGP/     独立上游checkout（忽略；固定SHA见reference/README）
runs/                 每次运行的原始预测与完整记录（忽略）
```

已核验HDF5为1135×10,709,949，且全文件只含0/1；全库 AraPheno feasibility registry 已生成并计入非有限值过滤，候选性状的来源级重复/汇总语义（390/261/262 组）仍按 registry TSV 状态保持未决。运行时有 106 项本地测试（CI 侧 105 passed + 1 skipped，被跳过的检查依赖 61 MB 的 catalog 工件，该工件不入 Git），含物理排序、样本列序、group 嵌套隔离、二态全位点翻转不变性、minority-state 边界扫描、formal v2 provenance binding（6 项运行重算 + record-only 交叉绑定）与 int8 ALT-dosage 存储契约（GT round-trip fixture，-1 不进入任何 dosage 统计）。**剩余 P0-1（formal VCF→int8 ALT-dosage conversion 尚未实施）**：int8、值域 `{-1,0,1,2}`、shape `sample×variant`，禁止先经 float32 全量矩阵生成；raw-file 审计由 build 阶段 `at_pheno.formal_provenance.verify_source_vcf` 承担，执行 gate 不重读 19 GB 源 VCF。HDF5 尚无完整 density prediction runner；PC-only、LD/RKHS/多核、标准GBLUP外部一致性及结构化SSM也尚未实现；不以合成自检或理论复杂度宣称真实全密度预测成功。
