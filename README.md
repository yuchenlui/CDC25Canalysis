# CDC25C Pan-Cancer Analysis

CDC25C 单基因泛癌全流程分析（单文件 R 脚本）：表达差异 → 预后（OS/DSS/DFI/PFI）→ TMB/MSI → 免疫微环境（ESTIMATE / CIBERSORT / TIMER2 / 免疫功能）→ Hallmark GSVA → 基因相关 / GSEA → 药敏（pRRophetic）→ IMvigor210 免疫治疗响应。

## 使用

```bash
Rscript CDC25C.pancancer.analysis.R
```

在包含各参考数据目录的项目根下执行（见脚本头部注释）。

## 数据下载

唯一表达量输入（约 9.5 GB）：UCSC Xena Toil Recomputed 泛癌 TPM 大表

- 下载页面：<https://xenabrowser.net/datapages/?cohort=TCGA%20Target-GTEx>
- 选择 dataset：`TcgaTargetGtex_rsem_gene_tpm`（Hugo symbol 版本，文件名 `TcgaTargetGtex_rsem_gene_tpm.symbol.tsv`）
- 放置到 `data/TcgaTargetGtex_rsem_gene_tpm.symbol.tsv`

脚本自动从该大表完成：
1. `data/CDC25C.tpm.denseDataOnlyDownload.tsv`（表达差异分析输入，含 TCGA/GTEx、样本类型注释）
2. `data/CDC25CexprTime.txt`（预后及下游各模块输入，仅 TCGA 样本）

## 输入数据布局

```
<data root>/
├── data/TcgaTargetGtex_rsem_gene_tpm.symbol.tsv
├── 07.preSurvival/ 13.cliCor/ 15.TMBcor/ 17.MSIcor/
├── 20.estimateCor/ 21.CIBERSORT/ 23.geneCor/ 25.GSEA/
├── 29.drugsensitivityIC50_pRRophetic/symbol.txt
└── 33.IMvigor210/IMvigor210/
```

## 输出

- 各步骤结果：`<步骤目录>/<输出目录>/`
- 汇总副本：`report/CDC25C/`

## 说明

- 表达口径为 Toil `log2(TPM+1)`；高低分组按中位数，与 FPKM 口径分组等价
- 步骤逐个 `tryCatch`：单步失败记录后继续，不中断整体流程
