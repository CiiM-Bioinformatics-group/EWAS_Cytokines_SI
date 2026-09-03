# EWAS of stimulation-induced cytokine responses

Analysis code for a study of how DNA methylation shapes inter-individual
variation in cytokine responses to immune stimulation.

Cytokine responses were profiled in the **SI cohort of older adults (n = 516)**
as discovery and the independent **BCG Prime cohort (n = 384)** for replication.
Genetic variation explained the largest share of response variance, particularly
for antiviral responses, while DNA methylation contributed a smaller but
consistent fraction and improved model performance in innate immune contexts.
Epigenome-wide association analyses identified **1,560 cytokine-associated CpG
sites (cCpGs)**, of which 45% replicated in BCG Prime and 28% remained
significant after multiple-testing correction.

## Repository layout

```
01_cytokines_QC_and_explained_variance_Figure2/
  01_QC_cytokines.R                     cytokine QC, filtering and rank-based normalisation
  01_explained_variance.R               variance partitioned between methylation and genetics

02_EWAS_run/
  EWAS_run.R                            per-cytokine epigenome-wide association analysis
  Validaiton_EWAS_results_BCG_PRIME.R   replication of discovery cCpGs in BCG Prime
```

Run the folders in order. Each script has a settings block at the top defining
its input and output paths.

## Requirements

R with `ggplot2`, `dplyr`, `tidyr`, `data.table`, `MASS`, `sandwich`, `lmtest`,
`foreach`, `doParallel`, `pheatmap`, `Hmisc`, `openxlsx`.

## Data

Individual-level methylation, genotype and cytokine data are not included in
this repository and are available from the cohort studies under their
respective access procedures.
