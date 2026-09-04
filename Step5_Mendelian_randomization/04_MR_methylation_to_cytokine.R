#!/usr/bin/env Rscript
# Two-sample MR: methylation at a CpG as exposure, cytokine response as outcome.
# One results file per cytokine.
#
# Each significant result is compared against the observational EWAS effect for
# the same CpG and cytokine, and labelled Concordant or Discordant on the sign.
#
# Usage:
#   Rscript 04_MR_methylation_to_cytokine.R <input_dir> <output_dir> <ewas_beta_file> [cytokine|ALL]

suppressMessages({
  library(dplyr)
  library(stringr)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 04_MR_methylation_to_cytokine.R <input_dir> <output_dir> <ewas_beta_file> [cytokine|ALL]")
}

cyto_dir       <- args[1]
out_dir        <- args[2]
ewas_beta_file <- args[3]
cytokine_arg   <- if (length(args) >= 4) args[4] else "ALL"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Alleles are the last two fields of the SNP identifier. Splitting on "_" rather
# than matching single characters keeps indels intact.
split_alleles <- function(snp_ids) {
  parts_list <- strsplit(snp_ids, "_")
  n <- lengths(parts_list)
  other  <- mapply(function(p, ni) if (ni >= 4) p[ni - 1] else NA_character_, parts_list, n)
  effect <- mapply(function(p, ni) if (ni >= 4) p[ni]     else NA_character_, parts_list, n)
  data.frame(other_allele = other, effect_allele = effect, stringsAsFactors = FALSE)
}


# --- observational EWAS effects, for the direction check ----------------------

cat("Reading EWAS betas:", ewas_beta_file, "\n")
if (!file.exists(ewas_beta_file)) stop("EWAS betas file not found: ", ewas_beta_file)

ewas_db <- read.table(ewas_beta_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(ewas_db) <- c("cpg", "ewas_beta", "file")

# Cytokine labels are spelled inconsistently across files, so normalise before
# matching on them
ewas_db <- ewas_db %>%
  mutate(
    file    = tolower(file),
    file    = gsub("cov[-\\.](ctrl|c)", "cov.ctrl", file),
    file    = gsub("cov[-\\.](n)", "cov.n", file),
    outcome = gsub("\\.cpgs\\.ewas\\.direction\\.txt$", "", file, ignore.case = TRUE)
  ) %>%
  select(cpg, outcome, ewas_beta)

cat("  Rows:", nrow(ewas_db),
    "| CpGs:", length(unique(ewas_db$cpg)),
    "| cytokines:", length(unique(ewas_db$outcome)), "\n\n")


# --- one cytokine -------------------------------------------------------------

run_mr_for_cytokine <- function(mqtl_file, gwas_file, cpgs_file, trait_label,
                                output_file, ewas_db) {

  cpg_list <- tryCatch(readLines(cpgs_file), error = function(e) NULL)
  if (is.null(cpg_list) || length(cpg_list) == 0) {
    cat("  Skipped, CpG list empty or unreadable\n"); return(invisible(NULL))
  }
  cpg_list <- trimws(cpg_list)
  cpg_list <- cpg_list[cpg_list != ""]

  # Exposure
  mqtl <- tryCatch(
    read.table(mqtl_file, header = FALSE, sep = "", stringsAsFactors = FALSE, quote = "\""),
    error = function(e) NULL)
  if (is.null(mqtl) || ncol(mqtl) != 5) {
    cat("  Skipped, mQTL file not in 5 columns\n"); return(invisible(NULL))
  }
  colnames(mqtl) <- c("cpg", "snp", "beta", "se", "pvalue")
  mqtl$cpg <- gsub('"', '', mqtl$cpg)

  mqtl <- mqtl[mqtl$cpg %in% cpg_list, ]
  if (nrow(mqtl) == 0) {
    cat("  Skipped, no mQTL rows for this cytokine's CpGs\n"); return(invisible(NULL))
  }

  mqtl_alleles       <- split_alleles(mqtl$snp)
  mqtl$other_allele  <- mqtl_alleles$other_allele
  mqtl$effect_allele <- mqtl_alleles$effect_allele
  mqtl <- mqtl[!is.na(mqtl$other_allele) & !is.na(mqtl$effect_allele), ]
  if (nrow(mqtl) == 0) {
    cat("  Skipped, no valid mQTL rows after allele parsing\n"); return(invisible(NULL))
  }

  exposure_dat <- format_data(
    mqtl, type = "exposure",
    snp_col = "snp", beta_col = "beta", se_col = "se", pval_col = "pvalue",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    phenotype_col = "cpg")

  # Outcome
  GWAS <- tryCatch(
    read.table(gwas_file, header = FALSE, sep = "", stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(GWAS) || ncol(GWAS) != 5) {
    cat("  Skipped, GWAS file not in 5 columns\n"); return(invisible(NULL))
  }
  colnames(GWAS) <- c("SNP", "gene", "beta", "t.stat", "p.value")

  GWAS_alleles       <- split_alleles(GWAS$SNP)
  GWAS$other_allele  <- GWAS_alleles$other_allele
  GWAS$effect_allele <- GWAS_alleles$effect_allele
  GWAS <- GWAS[!is.na(GWAS$other_allele) & !is.na(GWAS$effect_allele), ]
  GWAS$se <- abs(GWAS$beta / GWAS$t.stat)
  if (nrow(GWAS) == 0) {
    cat("  Skipped, no valid GWAS rows after allele parsing\n"); return(invisible(NULL))
  }

  outcome_dat <- format_data(
    GWAS, type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "p.value",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    phenotype_col = "gene")

  # action = 2 infers the strand for palindromic SNPs from allele frequency and
  # drops the ones that stay ambiguous
  dat <- tryCatch(harmonise_data(exposure_dat, outcome_dat, action = 2),
                  error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) {
    cat("  Skipped, harmonisation returned nothing\n"); return(invisible(NULL))
  }

  snp_counts <- dat %>%
    filter(mr_keep == TRUE) %>%
    group_by(id.exposure, id.outcome) %>%
    summarise(n_snp = n(), .groups = "drop")

  if (nrow(snp_counts) == 0) {
    cat("  Skipped, no valid instruments after harmonisation\n"); return(invisible(NULL))
  }

  mr_results <- tryCatch(mr(dat), error = function(e) NULL)
  if (is.null(mr_results) || nrow(mr_results) == 0) {
    cat("  Skipped, mr() returned nothing\n"); return(invisible(NULL))
  }

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
    cat("  Skipped, no Wald ratio or IVW results\n"); return(invisible(NULL))
  }

  # Sensitivity tests, pairs with 3 or more instruments
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

  mr_filtered <- mr_selected %>%
    filter(
      pval < 0.05,
      n_snp < 3 | (
        (is.na(pleiotropy_note)    | pleiotropy_note    == "Not significant") &
          (is.na(heterogeneity_note) | heterogeneity_note == "Not significant")
      )
    )

  if (nrow(mr_filtered) == 0) {
    cat("  No results cleared the filters\n"); return(invisible(NULL))
  }

  # Does the causal estimate point the same way as the observational EWAS effect
  mr_filtered <- mr_filtered %>%
    mutate(
      outcome_norm = tolower(outcome),
      outcome_norm = gsub("cov[-\\.](ctrl|c)", "cov.ctrl", outcome_norm),
      outcome_norm = gsub("cov[-\\.](n)", "cov.n", outcome_norm)
    ) %>%
    left_join(ewas_db, by = c("exposure" = "cpg", "outcome_norm" = "outcome")) %>%
    mutate(
      direction_consistent = case_when(
        is.na(ewas_beta)           ~ NA_character_,
        sign(b) == sign(ewas_beta) ~ "Concordant",
        TRUE                       ~ "Discordant"
      )
    ) %>%
    select(-outcome_norm)

  cat("  EWAS direction:",
      sum(mr_filtered$direction_consistent == "Concordant", na.rm = TRUE), "concordant,",
      sum(mr_filtered$direction_consistent == "Discordant", na.rm = TRUE), "discordant,",
      sum(is.na(mr_filtered$direction_consistent)), "no EWAS match\n")

  write.table(mr_filtered, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved", nrow(mr_filtered), "result(s)\n")

  invisible(mr_filtered)
}


# --- which cytokines to run ---------------------------------------------------

if (toupper(cytokine_arg) == "ALL") {
  cpgs_files   <- list.files(cyto_dir, pattern = "\\.cpgs\\.txt$", full.names = TRUE)
  trait_labels <- sub("\\.cpgs\\.txt$", "", basename(cpgs_files))
  cat("Running all", length(trait_labels), "cytokines\n\n")
} else {
  trait_labels <- cytokine_arg
  cat("Running one cytokine:", trait_labels, "\n\n")
}

n_done    <- 0
n_skip    <- 0
n_missing <- 0

for (trait_label in trait_labels) {

  out_file <- file.path(out_dir, paste0(trait_label, ".MR_results.txt"))

  # Existing results are left alone, so an interrupted batch can be resumed
  if (file.exists(out_file)) {
    cat("Already done:", trait_label, "\n")
    n_skip <- n_skip + 1
    next
  }

  mqtl_file <- file.path(cyto_dir, paste0(trait_label, ".mqtl.rebuilt.txt"))
  gwas_file <- file.path(cyto_dir, paste0(trait_label, ".GWAS.rebuilt.txt"))
  cpgs_file <- file.path(cyto_dir, paste0(trait_label, ".cpgs.txt"))

  if (!file.exists(mqtl_file) || !file.exists(gwas_file) || !file.exists(cpgs_file)) {
    cat("Missing input files:", trait_label, "\n")
    n_missing <- n_missing + 1
    next
  }

  cat("Processing:", trait_label, "\n")

  result <- tryCatch(
    run_mr_for_cytokine(mqtl_file, gwas_file, cpgs_file, trait_label, out_file, ewas_db),
    error = function(e) {
      cat("  Error:", conditionMessage(e), "\n")
      NULL
    })

  if (!is.null(result) && nrow(result) > 0) n_done <- n_done + 1
}

cat("\nCytokines requested:", length(trait_labels), "\n")
cat("Already done:       ", n_skip, "\n")
cat("Missing inputs:     ", n_missing, "\n")
cat("Produced results:   ", n_done, "\n")
cat("Results in:", out_dir, "\n")
