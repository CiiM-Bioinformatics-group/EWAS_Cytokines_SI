#!/usr/bin/env Rscript
# Two-sample MR: methylation at a CpG as exposure, disease or immune trait as
# outcome. Instruments are the cis-mQTL SNPs for that CpG.
#
# Pairs are kept whatever the instrument count, with the method matched to it:
#   1 SNP    Wald ratio, no sensitivity tests
#   2 SNPs   IVW, no sensitivity tests
#   3+ SNPs  IVW plus MR-Egger intercept, Cochran's Q and leave-one-out

#
# Usage:
#   Rscript 02_MR_methylation_to_disease.R <mqtl_file> <gwas_file> <out.txt> [prune_file|NONE]

library(dplyr)
library(stringr)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02_MR_methylation_to_disease.R <mqtl_file> <gwas_file> <out.txt> [prune_file|NONE]")
}

mqtl_file   <- args[1]
GWAS_file   <- args[2]
output_file <- args[3]
prune_file  <- if (length(args) >= 4) args[4] else "NONE"

SKIP_PRUNING <- toupper(prune_file) == "NONE"

# Alleles are the last two fields of the SNP identifier. Splitting on "_" rather
# than matching single characters keeps indels intact.
split_alleles <- function(snp_ids) {
  parts_list <- strsplit(snp_ids, "_")
  n <- lengths(parts_list)
  other  <- mapply(function(p, ni) if (ni >= 4) p[ni - 1] else NA_character_, parts_list, n)
  effect <- mapply(function(p, ni) if (ni >= 4) p[ni]     else NA_character_, parts_list, n)
  data.frame(other_allele = other, effect_allele = effect, stringsAsFactors = FALSE)
}


# --- exposure: mQTL -----------------------------------------------------------
# Layout: SNP  CpG  beta  se  pvalue   (whitespace separated, no header)

cat("Reading mQTL file:", mqtl_file, "\n")

mqtl <- read.table(mqtl_file, header = FALSE, sep = "", stringsAsFactors = FALSE, quote = "\"")
if (ncol(mqtl) != 5) {
  stop("mQTL file parsed into ", ncol(mqtl), " columns instead of 5. Check the delimiter.")
}
colnames(mqtl) <- c("snp", "cpg", "beta", "se", "pvalue")
mqtl$cpg <- gsub('"', '', mqtl$cpg)

mqtl_alleles       <- split_alleles(mqtl$snp)
mqtl$other_allele  <- mqtl_alleles$other_allele
mqtl$effect_allele <- mqtl_alleles$effect_allele

n_before <- nrow(mqtl)
mqtl <- mqtl[!is.na(mqtl$other_allele) & !is.na(mqtl$effect_allele), ]
cat("  mQTL rows:", nrow(mqtl), "(dropped", n_before - nrow(mqtl), "with malformed SNP IDs)\n")
cat("  Unique CpGs:", length(unique(mqtl$cpg)), "\n\n")


# --- outcome: disease GWAS ----------------------------------------------------
# Layout: SNP  trait  beta  t.stat  p.value   (tab separated, no header)

cat("Reading GWAS file:", GWAS_file, "\n")

GWAS <- read.table(GWAS_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(GWAS) != 5) {
  stop("GWAS file parsed into ", ncol(GWAS), " columns instead of 5. Check the delimiter.")
}
colnames(GWAS) <- c("SNP", "gene", "beta", "t.stat", "p.value")

GWAS_alleles       <- split_alleles(GWAS$SNP)
GWAS$other_allele  <- GWAS_alleles$other_allele
GWAS$effect_allele <- GWAS_alleles$effect_allele
GWAS <- GWAS[!is.na(GWAS$other_allele) & !is.na(GWAS$effect_allele), ]
GWAS$se <- abs(GWAS$beta / GWAS$t.stat)

cat("  GWAS rows:", nrow(GWAS), "\n")

if (SKIP_PRUNING) {
  cat("  No LD pruning applied.\n\n")
} else {
  cat("Reading LD pruning list:", prune_file, "\n")
  if (!file.exists(prune_file)) stop("Pruning file not found: ", prune_file)

  # Entries look like "chr17:130664:T:C;rs9674640"
  pruned_snp_ids <- unique(gsub(":", "_", sub(";.*$", "", readLines(prune_file))))
  cat("  Independent SNPs in list:", length(pruned_snp_ids), "\n")

  n_before_prune <- nrow(GWAS)
  GWAS <- GWAS[GWAS$SNP %in% pruned_snp_ids, ]
  cat("  GWAS rows after pruning:", nrow(GWAS),
      "(dropped", n_before_prune - nrow(GWAS), ")\n\n")
}

if (nrow(GWAS) == 0) {
  cat("No GWAS rows remain. Exiting.\n")
  quit(save = "no", status = 0)
}


# --- harmonise ----------------------------------------------------------------
# action = 2 infers the strand for palindromic SNPs from allele frequency and
# drops the ones that stay ambiguous.

exposure_dat <- format_data(
  mqtl, type = "exposure",
  snp_col = "snp", beta_col = "beta", se_col = "se", pval_col = "pvalue",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  phenotype_col = "cpg")

outcome_dat <- format_data(
  GWAS, type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "p.value",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  phenotype_col = "gene")

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
cat("Harmonised rows:", nrow(dat), "\n")


# --- run MR -------------------------------------------------------------------

snp_counts <- dat %>%
  filter(mr_keep == TRUE) %>%
  group_by(id.exposure, id.outcome) %>%
  summarise(n_snp = n(), .groups = "drop")

cat("CpG-outcome pairs by instrument count:\n")
cat("  1 SNP  :", sum(snp_counts$n_snp == 1), "\n")
cat("  2 SNPs :", sum(snp_counts$n_snp == 2), "\n")
cat("  3+ SNPs:", sum(snp_counts$n_snp >= 3), "\n\n")

if (nrow(snp_counts) == 0) {
  cat("No CpG had valid instruments after harmonisation. Exiting.\n")
  quit(save = "no", status = 0)
}

mr_results <- mr(dat)
print(mr_results)

mr_results$method <- trimws(as.character(mr_results$method))
mr_results <- mr_results %>% left_join(snp_counts, by = c("id.exposure", "id.outcome"))

mr_selected <- mr_results %>%
  filter(
    (n_snp == 1 & method == "Wald ratio") |
      (n_snp >= 2 & method == "Inverse variance weighted")
  ) %>%
  mutate(
    instrument_tier = case_when(
      n_snp == 1 ~ "1 SNP (Wald ratio, no sensitivity)",
      n_snp == 2 ~ "2 SNPs (IVW, no sensitivity)",
      TRUE       ~ "3+ SNPs (IVW + full sensitivity)"
    )
  )

if (nrow(mr_selected) == 0) {
  cat("No Wald ratio or IVW results found. Exiting.\n")
  quit(save = "no", status = 0)
}


# --- sensitivity tests, pairs with 3 or more instruments ----------------------

dat_multi <- dat %>%
  semi_join(snp_counts %>% filter(n_snp >= 3), by = c("id.exposure", "id.outcome"))

if (nrow(dat_multi) > 0) {
  pleio <- tryCatch(mr_pleiotropy_test(dat_multi), error = function(e) NULL)
  het   <- tryCatch(mr_heterogeneity(dat_multi),   error = function(e) NULL)
  loo   <- tryCatch(mr_leaveoneout(dat_multi),     error = function(e) NULL)
} else {
  pleio <- het <- loo <- NULL
}

mr_selected$pleiotropy_pval    <- NA
mr_selected$pleiotropy_note    <- NA
mr_selected$heterogeneity_pval <- NA
mr_selected$heterogeneity_note <- NA
mr_selected$loo_min_pval       <- NA
mr_selected$loo_max_pval       <- NA
mr_selected$loo_note           <- NA

if (!is.null(pleio) && nrow(pleio) > 0) {
  mr_selected <- mr_selected %>%
    left_join(pleio %>% select(id.exposure, id.outcome, pval) %>%
                rename(pleiotropy_pval_new = pval),
              by = c("id.exposure", "id.outcome")) %>%
    mutate(
      pleiotropy_pval = coalesce(pleiotropy_pval_new, pleiotropy_pval),
      pleiotropy_note = case_when(
        n_snp < 3                                        ~ NA_character_,
        !is.na(pleiotropy_pval) & pleiotropy_pval < 0.05 ~ "Significant",
        TRUE                                             ~ "Not significant"
      )
    ) %>%
    select(-pleiotropy_pval_new)
}

if (!is.null(het) && nrow(het) > 0) {
  het_ivw <- het %>% filter(method == "Inverse variance weighted")
  if (nrow(het_ivw) > 0) {
    mr_selected <- mr_selected %>%
      left_join(het_ivw %>% select(id.exposure, id.outcome, Q_pval) %>%
                  rename(het_pval_new = Q_pval),
                by = c("id.exposure", "id.outcome")) %>%
      mutate(
        heterogeneity_pval = coalesce(het_pval_new, heterogeneity_pval),
        heterogeneity_note = case_when(
          n_snp < 3                                              ~ NA_character_,
          !is.na(heterogeneity_pval) & heterogeneity_pval < 0.05 ~ "Significant",
          TRUE                                                   ~ "Not significant"
        )
      ) %>%
      select(-het_pval_new)
  }
}

if (!is.null(loo) && nrow(loo) > 0) {
  loo_summary <- loo %>%
    filter(SNP != "All") %>%
    group_by(id.exposure, id.outcome) %>%
    summarise(loo_min_pval_new = min(p, na.rm = TRUE),
              loo_max_pval_new = max(p, na.rm = TRUE), .groups = "drop")

  mr_selected <- mr_selected %>%
    left_join(loo_summary, by = c("id.exposure", "id.outcome")) %>%
    mutate(
      loo_min_pval = coalesce(loo_min_pval_new, loo_min_pval),
      loo_max_pval = coalesce(loo_max_pval_new, loo_max_pval),
      loo_note = case_when(
        n_snp < 3                                  ~ NA_character_,
        !is.na(loo_max_pval) & loo_max_pval < 0.05 ~ "Stable",
        TRUE                                       ~ "Driven by single SNP"
      )
    ) %>%
    select(-loo_min_pval_new, -loo_max_pval_new)
}


# --- keep the nominally significant results -----------------------------------

mr_filtered <- mr_selected %>%
  filter(
    pval < 0.05,
    n_snp < 3 | (
      (is.na(pleiotropy_note)    | pleiotropy_note    == "Not significant") &
        (is.na(heterogeneity_note) | heterogeneity_note == "Not significant")
    )
  )

cat("\nResults passing filters:\n")
cat("  1-2 instruments (p<0.05):", sum(mr_filtered$n_snp < 3), "\n")
cat("  3+ instruments (p<0.05, no pleiotropy or heterogeneity):",
    sum(mr_filtered$n_snp >= 3), "\n")
cat("  Total:", nrow(mr_filtered), "\n")
cat("  Pruning applied:", !SKIP_PRUNING, "\n")

if (nrow(mr_filtered) > 0) {
  write.table(mr_filtered, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved to:", output_file, "\n")
} else {
  cat("No results passed the filters. No output written.\n")
}
