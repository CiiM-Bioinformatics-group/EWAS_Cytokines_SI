# EWAS of stimulation-induced cytokine responses

Analysis code for a study of the role of DNA methylation in inter-individual
variation in cytokine responses to immune stimulation, using the SI cohort of
older adults as discovery and the independent BCG Prime cohort for replication.


## Repository layout

```
01_cytokines_QC_and_explained_variance_Figure2/
  01_QC_cytokines.R                     cytokine QC, filtering and rank-based normalisation
  01_explained_variance.R               variance partitioned between methylation and genetics

02_EWAS_run/
  EWAS_run.R                            per-cytokine epigenome-wide association analysis
  Validaiton_EWAS_results_BCG_PRIME.R   replication of discovery cCpGs in BCG Prime

03_longituidal_methatlion_change_AND_GWAS_EWAS_PUBLIC_DATA_ENRICHEMT_FIGURE3/
  delta_longitudinal_methylation_change.R
                                        methylation change at cCpGs across three time points
                                        in 300BCG, tested against matched background
```

The cCpGs are identified in the SI cohort. To ask whether these sites change following immune perturbation, 
they are carried into the independent 300BCG cohort, which received BCG
vaccination and was sampled before vaccination, 14 days after and 3 months
after.

Run the folders in order. Each script has a settings block at the top defining
its input and output paths.

## Requirements

R with `ggplot2`, `dplyr`, `tidyr`, `data.table`, `MASS`, `sandwich`, `lmtest`,
`foreach`, `doParallel`, `pheatmap`, `Hmisc`, `openxlsx`.

## Data

Individual-level methylation, genotype and cytokine data are not included in
this repository and are available from the cohort studies under their
respective access procedures.
