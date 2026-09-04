# EWAS of stimulation-induced cytokine responses

Analysis code for a study of the role of DNA methylation in inter-individual
variation in cytokine responses to immune stimulation, using the SI cohort of
older adults as discovery and the independent BCG Prime cohort for replication.


## Repository layout

```
01_cytokines_QC_and_explained_variance_Figure2/
    cytokine QC and normalisation, and variance partitioned between
    methylation and genetics

02_EWAS_run/
    per-cytokine epigenome-wide association analysis, and replication
    of the discovery cCpGs in BCG Prime

step3_Figure3_analysis/
    longitudinal methylation change at cCpGs in 300BCG, and enrichment
    of cCpGs in public GWAS and EWAS trait data

Step4_Figure4_genomic_features_enrichment_by_direction/
    genomic feature annotation of cCpGs and enrichment by direction of
    effect, against the non-significant array background

Step5_Mendelian_randomization/
    two-sample MR of methylation on cytokine responses and on disease
    and immune traits
```

The cCpGs are identified in the SI cohort. To ask whether these sites change following immune perturbation, 
they are carried into the independent 300BCG cohort, which received BCG
vaccination and was sampled before vaccination, 14 days after and 3 months
after.

Run the folders in order. Each script has a settings block at the top defining
its input and output paths.

## Requirements

R with `ggplot2`, `dplyr`, `tidyr`, `data.table`, `MASS`, `sandwich`, `lmtest`,
`foreach`, `doParallel`, `pheatmap`, `Hmisc`, `openxlsx`, `TwoSampleMR` and
`GenomicRanges`.

## Data

Individual-level methylation, genotype and cytokine data are not included in
this repository and are available from the cohort studies under their
respective access procedures.
