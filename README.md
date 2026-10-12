# CDC25C Pan-Cancer Analysis

Single-file R pipeline for single-gene pan-cancer analysis of CDC25C: differential expression -> survival (OS/DSS/DFI/PFI) -> TMB/MSI -> immune microenvironment (ESTIMATE / CIBERSORT / TIMER2 / immune functions) -> Hallmark GSVA -> gene correlations / GSEA -> drug sensitivity (pRRophetic) -> IMvigor210 immunotherapy response.

## Usage

```bash
Rscript CDC25C.pancancer.analysis.R
```

Run from the project root that contains the per-step reference directories (see the script header).

## Data download

Expression input (single source, ~9.5 GB): UCSC Xena Toil Recomputed pan-cancer TPM matrix

- Download page: <https://xenabrowser.net/datapages/?cohort=TCGA%20Target-GTEx>
- Dataset: `TcgaTargetGtex_rsem_gene_tpm` (Hugo symbol version, file name `TcgaTargetGtex_rsem_gene_tpm.symbol.tsv`)
- Place at `data/TcgaTargetGtex_rsem_gene_tpm.symbol.tsv`

The script derives from this matrix:

1. `data/CDC25C.tpm.denseDataOnlyDownload.tsv` (input to differential expression; TCGA/GTEx samples with sample-type annotation)
2. `data/CDC25CexprTime.txt` (expression input to steps 08-29, TCGA samples only)

Survival data (OS/DSS/DFI/PFI) is not derived from the expression matrix. It is downloaded from UCSC Xena as `07.preSurvival/TCGA-<CANCER>.survival.tsv` and merged with expression by sample at runtime for steps 08-12.

## Input layout

```
<data root>/
├── data/TcgaTargetGtex_rsem_gene_tpm.symbol.tsv
├── 07.preSurvival/ 13.cliCor/ 15.TMBcor/ 17.MSIcor/
├── 20.estimateCor/ 21.CIBERSORT/ 23.geneCor/ 25.GSEA/
├── 29.drugsensitivityIC50_pRRophetic/symbol.txt
└── 33.IMvigor210/IMvigor210/
```

## Output

- Per-step results: `<step dir>/<output dir>/`
- Collected copy: `report/CDC25C/`

## Notes

- Expression is Toil `log2(TPM+1)`; High/Low groups are split at the median, equivalent to the FPKM-based grouping
- Each step runs inside `tryCatch`; a failed step is logged and the pipeline continues
