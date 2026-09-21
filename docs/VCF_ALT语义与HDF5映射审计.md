# Population VCF ALT 语义与 HDF5 坐标/状态映射审计（TASK02/TASK03）

本文是「数据语义收口阶段」TASK02/TASK03 的验收与证据链文档。只做审计
与有界 spot-check，**不生成全量训练矩阵**、**不建立 10M×1135 全量 NPY**、
**不把 HDF5 的 0/1 值称作 ALT、不做 dosage 转换、不写回 HDF5**。
全部数值引自 repo 内 JSON 与本轮命令实跑输出；证据文件清单见 §8。

## 1. TASK02 — 来源与校验（§2.1/2.2）

| 项 | 值 | 证据 |
| --- | --- | --- |
| URL | `https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz` | `1001g_v31_population_vcf_source.json` → `source_url` |
| 官方 MD5 | `a29dec1389cf31c2d18ecd76947d4e54` | 官方 `.md5` 侧文件 `data/raw/1001g_v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz.md5` |
| 本地文件 | 19,230,910,695 B（2026-09-22 04:30 mtime） | `stat` |
| 本地 MD5 | `a29dec1389cf31c2d18ecd76947d4e54`（== 官方，`md5_matches_official: true`） | source JSON |
| 本地 SHA256 | `7f825afb784b2424e35501fd8b88c32a4798013eb9fe762724addfb249007e7a` | source JSON `local_sha256` |
| 流编码 | gzip（读满 1 遍不解除压到明文临时文件） | source JSON `stream_encoding` + 2.4/3.2 流程 |
| bcftools | 1.24 `view -h` 可读：样本头行为 8 个固定列 + 字面量 `FORMAT` 列 + 1135 个样本列（88、108、139…，末 19949/19950/19951） | 2026-09-22 06:57 实跑：`bcftools --version` = 1.24 + `view -h` 成功 |

下载时点如实记录：`download_time = "aria2c session completed 2026-09-21T20:30:30Z (file mtime UTC; mis-typed 1970 session stamp discarded)"`
（先前误记的 1970 时戳为错误时间戳，弃用；`download_tool` 记录 aria2c 续传方式）。
本轮主扫描（06:50 运行）未传 `--bcftools`，故 source JSON `bcftools_probe` 为 null；
2.2 的 bcftools 可读性以 §1 末行的一次性实跑为准，不重烧全量扫描。

## 2. TASK02 — Header 与样本（§2.3）

- `##FORMAT=<ID=GT,...>` 在 1001G 官方文件中用**逗号**分隔属性；审计模块
  `scripts/audit_1001g_vcf.py` 的 meta 解析按 `;`/`,` 分离（fixture：
  `test_vcf_meta_comma_separators_parse_gt_and_filter`），修后 `gt_column_index = 0`。
- `#CHROM` 行为 1144 列，第 8 列是**字面量 token `FORMAT`**（bcftools `-h`
  输出与 header JSON `format_column_in_chrom: true` 双重佐证），样本列 =
  第 9 列起，`n_samples = 1135`。此数由审计得出，**非假定**。
  全量数据行样本字段宽度一致：`sample_fields_width_mismatch: null`。
- 与官方 1135 表的映射（候选 3 键，结果在 `1001g_v31_vcf_header_sample_audit.json` →
  `sample_id_map`）：
  - `accession_id`（表字段 `Accession_ID`）：双侧 set 相等、双侧唯一、各 1135 →
    **resolved**；`order_match: false`，`order_divergences: 761` —— 行序
    不收入 gate，以名称键 1:1 对应为映射成立的标准（§2.3 语义）；
  - `accession_name` / `cs_number` 候选键：unresolved，差集已记录（不 fuzzy、
    不猜生态型名）。

## 3. TASK02 §2.4 窗口 + TASK03 §3.1/3.2 目录

窗口：5 染色体 × 每 chr 头 400 条核 SNP records = 2,000 条、2,270,000 个 GT 单元格。
GT 判定全程 as-is：GT 单元格等位基因数直方图 {1 等位基因: 0, 2 等位基因: 2,270,000, >2: 0}
（无 haploid 单等位基因 GT，无 >2）；allele-index 全部在 0/1，`out_of_range_or_unexpected: 0`。
（交付包 §2.4 连带：不得仅凭此把 1 改成二倍性 dosage 2；GT ploidy/dosage 判断保持
UNRESOLVED，不自动推断。）
`format_column_assumed: true` 已记录。

§3.1/§3.2 全量 12,883,854 条记录的 7-bucket 计数（见 catalog summary JSON）：

| bucket | 数量 |
| --- | --- |
| nuclear_biallelic_snp | 10,707,430 |
| multiallelic（不分解，整体排除） | 1,113,934 |
| indel_or_complex | 1,062,490 |
| n_or_noncanonical | 0 |
| malformed | 0 |
| non_nuclear | 0 |
| ref_equals_alt | 0 |

duplicate `(chrom,pos)` = 0、duplicate `(chrom,pos,ref,alt)` 四元组 = 0
（流式一次扫描生成 catalog gz，再用一次回读校验重复键值）。
位点池：交付包规定的固定 hash/salt `md5(chr:pos:at_pheno_task03_v1)`，cap = 2,000
（1,000–2,000 区间上限），在 1:1 唯一匹配坐标集上选出 2,000 站；mapping 侧复核
`pool_sites_missing_in_vcf = 0` —— 2,000 站全部落进 10,707,430 行 catalog。

## 4. TASK03 §3.3 坐标映射

HDF5 `data/raw/1001g_v3.1/extracted/imputed_snps_binary.hdf5`（`sha256 f3713e66f017187fb3d30b8d6b854794393ab259e6b5bae4b18ca2cafcbb1862`）：
`positions` 轴 10,709,949 行（int32），按 5 个连续 `[start, stop)` region 段
组织（轴属性 `chr_regions`），chrom 标签 = `['1','2','3','4','5']`，轴上 0 重复坐标；
序列化 manifest：`data/manifests/1001g_hdf5_2026-09-22.json`。

streaming merge / sorted join（`scripts/audit_hdf5_vcf_mapping.py`，§3.3 六项）：

| 数量 | 值 |
| --- | --- |
| `hdf5_coordinate_rows` | 10,709,949 |
| `unique_matches` | 10,700,811 |
| `multiple_candidates` | 0 |
| `hdf5_only_coordinates` | 9,138 |
| `vcf_only_coordinates` | 6,619 |
| `coord_hits` | 10,700,811 |

§3.3 不变量：坐标唯一匹配**不推出** HDF5 0/1 的 REF/ALT 方向（summary JSON
`no_orientation_claimed: true`）。

## 5. TASK03 §3.4 Spot（2,000 站 × 1,135 共享样本）

共同样本事实（记录，非推断）：HDF5 `accessions` 行序 == VCF `#CHROM` 样本头序
（排除第 8 列字面量 `FORMAT`） == accession 升序，1,135/1,135 一一对应；官方
1135 表 `Accession_ID` 列的集合相同（表的行序不同，因此表行序不作序证据）。
样本重复/共享上下文（`--sample-file` JSON）记录此事，`shared_fraction = 1.0`；
SAMPLE-ORDERING 仍按规范保持人工 cross-check 事实，不作推断。

两侧 as-is 计数（2,270,000 单元格 = 2,000 × 1,135）：

| 侧 | 构成 |
| --- | --- |
| HDF5 int8 观测 | `{0: 2,133,944, 1: 136,056}`，无 0/1 以外取值（`out_of_0_1_cells = 0`） |
| VCF 状态码 | homozygous_0: 1,793,153；homozygous_1: 94,663；heterozygous: **0**；missing: 382,184；out_of_range_or_unexpected: 0 |

配对（坐标两侧均有值、且 VCF 侧给出逐样本行）共 1,887,816 对
（= VCF 侧全部 homozygous-非缺失样本；异合 0，故未出现异合对）；
总体读取（具体见 `hdf5_vcf_concordance_spot_check.json` → `overall`）：

```text
concordance_if_hdf5_0_is_REF = 0.9999629201150959（≈ 1,887,746 / 1,887,816，失配 ≈ 70 对）
concordance_if_hdf5_1_is_REF = 3.7079884904037255e-05
best_orientation = "hdf5_0_is_ref"（仅为标签；不写回 HDF5、不生成 ALT-dosage）
n_sites_with_pairs = 2000
per-site 分布：1,999 站 0-is-ref / 1 站 1-is-ref / 0 平局（无整站全缺失）
suspect_sites = []（无悬案位点；合计为 0 的两类逐站枚举：全缺失站 / 平局站）
```

如实申报（不做阈值、不做强行）：2,000 站窗口内 VCF 侧 **heterozygous 状态为
0**，missing 计 382,184 / 2,270,000（≈16.8%）。这是原始观测记录；GT 策略与
imputation 覆盖的解释由下一轮数据评审承担，本阶段**不**据此升降级 gate。

## 6. Gate 判定（§3.4 验收 6 项，左右对齐）

| 验收项 | 本轮实跑结果 |
| --- | --- |
| 官方 source checksum 正确 | `md5_matches_official: true`；SHA256 `7f825a…`（§1） |
| sample mapping 无歧义 | `accession_id` resolved（set 1:1、双侧唯一；行序差 761 已记录） |
| GT ploidy / ALT-dosage 语义明确 | 窗口全 2 等位基因几何；as-is 记录，无 dosage、无 haploid 误判 |
| variant filter policy 可复算 | 审计脚本路径 + summary 7-bucket + duplicate 计数 0 + commit 完整 |
| HDF5↔VCF 坐标映射统计完成 | §4 六项全部计算；`hdf5_only`/`vcf_only` 已如实列出 |
| spot-check 未暴露系统性错位 | `suspect_sites: []`；已配对 1,887,816 对；方向标签 0.99996；无阈值强通 |

→ 两侧 gate 源同时给出：
- VCF 侧简行 JSON：`{"mapping": "resolved", "dup_chr_pos": 0, "dup_quads": 0, "catalog_rows": 10707430, "n_picked": 2000, "run": "VCF_SEMANTICS_READY_FOR_CONVERSION", "exit_code": 0}`
- mapping 侧 gate（`hdf5_vcf_position_mapping_summary.json` → `gate`）：`{"status": "PASS", "gate": "VCF_SEMANTICS_READY_FOR_CONVERSION", "reasons": []}`

依据：TASK01/03 的可用证据与本阶段完成项归入后再收紧 TASK04（见
`docs/数据契约.md` / `docs/正式实验门禁.md` 的 `at_pheno_formal_v2` 节）。
748（可用的 panel 候选）仍保持**不冻结**；
390 / 261 / 262 不自动冻结。

## 7. 未过 / 延迟

- VCF 侧 2,000 站窗口 hetero=0 与 missing≈16.8%：记录在案；GT 策略 /
  imputation 覆盖的下一步解释由下一轮数据评审承担，本 gate **不**强制通过/不通过。
- multiallelic（1,113,934 条）按交付包 §3.1「先排除，不分解」策略处理；包文本同时注明
  「这是最小、可核验的 ALT 语义基线，不代表永久排除 multiallelic」。
- 方向标签 `hdf5_0_is_ref` 仅为 spot-check 结论标签：不写回、不生成 12 GB
  ALT-dosage `genotypes.npy`；跨版本对齐 / multiallelic 属性声明、SESOI、正式 trait panel
  仍 UNRESOLVED（P0-1 未冻）。
- kinship threshold / 显著性方案本阶段不冻结。
- 大文件不入 git：18 GB VCF、877 MB HDF5、61.6 MB catalog gz、53.7 MB VCF 侧
  spot JSON（`data/manifests/1001g_v31_concordance_spot_check.json`，per_sample
  逐样本行 2000×1135）与全量逐点 mapping 均在 `data/raw` / 本轮工作区；
  Git 只提交小的 JSON/TSV/文档/脚本（1.5 MB 的 H5 侧映射/spot 摘要 JSON 入仓库）。

## 8. 证据文件与复现

仓库内（本轮提交）：

```text
docs/VCF_ALT语义与HDF5映射审计.md   ← 本文
docs/性状来源与目标定义审计.md       ← TASK01（426de76 已提交；SESOI 措辞入提交③）
data/manifests/1001g_v31_population_vcf_source.json
data/manifests/1001g_v31_vcf_header_sample_audit.json
data/manifests/1001g_v31_gt_audit.json
data/manifests/1001g_v31_snp_catalog_summary.json
data/manifests/hdf5_vcf_position_mapping_summary.json
data/manifests/hdf5_vcf_concordance_spot_check.json
data/manifests/1001g_hdf5_2026-09-22.json
scripts/audit_1001g_vcf.py / scripts/audit_hdf5_vcf_mapping.py
tests/test_audit_1001g_vcf.py / tests/test_audit_hdf5_vcf_mapping.py
```

制品 sha256（2026-09-22）：

```text
VCF  1001genomes_snp-short-indel_only_ACGTN.vcf.gz
  md5      a29dec1389cf31c2d18ecd76947d4e54（== 官方）
  sha256   7f825afb784b2424e35501fd8b88c32a4798013eb9fe762724addfb249007e7a
HDF5 extracted/imputed_snps_binary.hdf5
  sha256   f3713e66f017187fb3d30b8d6b854794393ab259e6b5bae4b18ca2cafcbb1862
表   1135_accessions_table.json
  sha256   653c6b16c3562f34264d77f48a4f7555c5618011b72c14bd58c1ee44a5f2a3db
catalog gz（不入 git；mtime-dependent sha，见下）
  sha256   31649a02b053928bb4eb54065e551414f68e1747d185933801ffff5f206d3b59（本轮）
小 JSON（data/manifests，sha256）：
  1001g_v31_population_vcf_source.json       15f8bb935abf4cc94716d14553c1dedd80b7ea767a61e827bc5d53c35eb763c5
  1001g_v31_vcf_header_sample_audit.json     70694c5ef0db6ddd9ffd0ec60fcf66069f5ab0137994293428ebaf42a04cd211
  1001g_v31_gt_audit.json                    33a804a095cedc5468f604d61e1641d494132985671c8cd19c19b90b92ff2cc5
  1001g_v31_snp_catalog_summary.json         c84a99efed1cf7edcc12df2f12beb57d5a693527ab17dfbd8a91d7a343ed5729
  （仓库内小 JSON；53.7 MB 的 VCF 侧 spot JSON 未入 Git，其制品值同样录制：）
  1001g_v31_concordance_spot_check.json      5bdf639d298b2d9569e8f7876c51e65c8b237176b32c97ce9837bb3033dd67a5 [工作区制品，不入 Git]
  hdf5_vcf_position_mapping_summary.json      5cbed77e48040d3f5658c69435bbea5341dbcbc944e7d3701f36ccb5f5f4513d
  hdf5_vcf_concordance_spot_check.json        03d14a62b7a1d5514a83caa72db06b6f7a2ab27084d5340ef11c005c35062dbf
  1001g_hdf5_2026-09-22.json                  ef52a70c48eabbc46c26ae51372c880ec7de5430ea2a69c6ccbd84c3ab360af1
```

catalog gz 确定性说明（两遍同命令实跑，/tmp 侧 07:08 与 07:18）：
**展开内容逐字节一致**（展开 md5 均为 `829ccfcab58107f57fce4683858c10f9`）；
gz 整包 sha 仅相差 gzip 头 mtime 字段的 2 字节（全流字节 diff 即第 4-5 字节），
因此 catalog gz 整包 sha256（上文 `3164…`，含本轮 mtime）**不是**稳定完整性标记；
git 只提交小 JSON 与脚本；如需内容一致性复核，重跑脚本后比较展开内容（上条 md5）即可。

复现（同机重跑）：

```bash
# TASK02/03-part1 VCF 审计（全程约 11 min）
uv run python scripts/audit_1001g_vcf.py \
  --vcf data/raw/1001g_v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
  --out data/manifests \
  --accession-table data/raw/1001g_v3.1/1135_accessions_table.json \
  --official-md5 a29dec1389cf31c2d18ecd76947d4e54 \
  --source-url "https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz" \
  --download-time "aria2c session completed 2026-09-21T20:30:30Z (file mtime UTC; mis-typed 1970 session stamp discarded)" \
  --download-tool "aria2c -c -x8 -s8 (resumed via control file)" \
  --bcftools /home/yhshy/miniconda3/envs/audit-1001g/bin/bcftools \
  --picked-sites <mapping#1 summary JSON: position_mapping.site_pool_sample>

# TASK03 part2 坐标映射
uv run python scripts/audit_hdf5_vcf_mapping.py \
  --hdf5 data/raw/1001g_v3.1/extracted/imputed_snps_binary.hdf5 \
  --catalog data/manifests/1001g_v31_biallelic_snp_catalog.tsv.gz \
  --per-site-vcf data/manifests/1001g_v31_concordance_spot_check.json \
  --hdf5-values <2000 站 × 1135 列的 per-sample 0/1 状态 JSON（嵌套 {chrom: {pos: [...]}}）> \
  --sample-file <共享样本上下文档
（shared_accessions / shared_n / n_samples_vcf / n_samples_hdf5 / note，note 含序证据说明）> \
  --hdf5-manifest data/manifests/1001g_hdf5_2026-09-22.json \
  --out data/manifests
```

（两个 `<...>` 输入是 `/tmp` 侧一次性生成物，非仓库文件；生成方法与
`scripts/audit_hdf5_vcf_mapping.py` 的纯函数对应。）
注：本条命令里的 `--bcftools` 属 2.2 可选探针入参；本轮 06:50 主扫描**未**传该参数
（故 source JSON `bcftools_probe` 为 null），可读性已由 §1 的一次性 `view -h` 实跑佐证，
未为它重烧全量扫描。
