#!/usr/bin/env Rscript
# Build the exposure and outcome files for the cytokine MR, one pair per cytokine.
#
#   exposure  the GoDMC cis-mQTLs for that cytokine's CpGs, taken from the
#             master mQTL file
#   outcome   that cytokine's own genome-wide cQTL results, restricted to the
#             SNPs the exposure needs
#
# The GoDMC instruments are already a clumped lead-SNP set, so no further LD
# pruning is applied by default. Pass a PLINK .prune.in file as the last
# argument to restrict the outcome side to an independent subset.
#
# Usage:
#   Rscript 03_build_cytokine_MR_inputs.R <cytokine|ALL> <cpg_list_dir> \
#     <master_mqtl_file> <cQTL_dir> <output_dir> [prune_file|NONE]

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 03_build_cytokine_MR_inputs.R <cytokine|ALL> <cpg_list_dir> <master_mqtl_file> <cQTL_dir> <output_dir> [prune_file|NONE]")
}

cytokine_arg     <- args[1]   # one trait label, or ALL
cyto_dir         <- args[2]   # holds <trait>.cpgs.txt
shared_mqtl_file <- args[3]
cqtl_dir         <- args[4]   # holds <trait>.snps.input.cyto.all.tsv
out_dir          <- args[5]
prune_file_arg   <- if (length(args) >= 6) args[6] else "NONE"

SKIP_PRUNING <- toupper(prune_file_arg) == "NONE"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- master mQTL file ---------------------------------------------------------
# Layout: SNP  CpG  beta  se  pvalue   (whitespace separated, CpG quoted)

cat("Reading master mQTL file:", shared_mqtl_file, "\n")

master_mqtl <- read.table(shared_mqtl_file, header = FALSE, sep = "",
                          stringsAsFactors = FALSE, quote = "\"")
if (ncol(master_mqtl) != 5) {
  stop("Master mQTL file parsed into ", ncol(master_mqtl), " columns instead of 5.")
}
colnames(master_mqtl) <- c("snp", "cpg", "beta", "se", "pvalue")
master_mqtl$cpg <- gsub('"', '', master_mqtl$cpg)

cat("  Rows:", nrow(master_mqtl),
    "| unique CpGs:", length(unique(master_mqtl$cpg)), "\n\n")


# --- optional pruning list ----------------------------------------------------

if (SKIP_PRUNING) {
  cat("No LD pruning applied.\n\n")
  pruned_snp_ids <- NULL
} else {
  cat("Reading LD pruning list:", prune_file_arg, "\n")
  if (!file.exists(prune_file_arg)) stop("Pruning file not found: ", prune_file_arg)
  pruned_snp_ids <- unique(gsub(":", "_", sub(";.*$", "", readLines(prune_file_arg))))
  cat("  Independent SNPs in list:", length(pruned_snp_ids), "\n\n")
}


# --- which cytokines to build -------------------------------------------------

if (toupper(cytokine_arg) == "ALL") {
  cpgs_files   <- list.files(cyto_dir, pattern = "\\.cpgs\\.txt$", full.names = TRUE)
  trait_labels <- sub("\\.cpgs\\.txt$", "", basename(cpgs_files))
  cat("Building all", length(trait_labels), "cytokines\n\n")
} else {
  trait_labels <- cytokine_arg
  cat("Building one cytokine:", trait_labels, "\n\n")
}


# --- one cytokine -------------------------------------------------------------

build_one_cytokine <- function(trait_label) {

  cat("----------------------------------------------------------\n")
  cat("Building inputs for:", trait_label, "\n")

  cpgs_file <- file.path(cyto_dir,  paste0(trait_label, ".cpgs.txt"))
  cqtl_file <- file.path(cqtl_dir,  paste0(trait_label, ".snps.input.cyto.all.tsv"))
  out_mqtl  <- file.path(out_dir,   paste0(trait_label, ".mqtl.rebuilt.txt"))
  out_gwas  <- file.path(out_dir,   paste0(trait_label, ".GWAS.rebuilt.txt"))

  if (!file.exists(cpgs_file)) {
    cat("  Skipped, CpG list not found:", cpgs_file, "\n\n"); return(invisible(NULL))
  }
  if (!file.exists(cqtl_file)) {
    cat("  Skipped, cQTL file not found:", cqtl_file, "\n\n"); return(invisible(NULL))
  }

  cpg_list <- trimws(readLines(cpgs_file))
  cpg_list <- cpg_list[cpg_list != ""]
  cat("  CpGs:", length(cpg_list), "\n")

  # Exposure
  exposure <- master_mqtl[master_mqtl$cpg %in% cpg_list, ]
  cat("  mQTL instrument rows:", nrow(exposure),
      "| CpGs matched:", length(unique(exposure$cpg)), "of", length(cpg_list), "\n")

  if (nrow(exposure) == 0) {
    cat("  Skipped, no mQTL rows for these CpGs\n\n"); return(invisible(NULL))
  }

  missing_cpgs <- setdiff(cpg_list, unique(exposure$cpg))
  if (length(missing_cpgs) > 0) {
    cat("  ", length(missing_cpgs), "CpG(s) have no instruments in the master file\n")
  }

  exposure_out <- exposure %>% select(cpg, snp, beta, se, pvalue)
  write.table(exposure_out, out_mqtl, sep = " ", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat("  Saved exposure:", basename(out_mqtl), "\n")

  needed_snps <- unique(exposure$snp)

  # Outcome
  cqtl <- tryCatch(
    read.table(cqtl_file, header = TRUE, sep = "", stringsAsFactors = FALSE,
               quote = "", comment.char = ""),
    error = function(e) NULL)

  if (is.null(cqtl) || ncol(cqtl) != 5) {
    cat("  Skipped, could not parse the cQTL file\n\n"); return(invisible(NULL))
  }
  colnames(cqtl) <- c("SNP_raw", "gene", "beta", "t.stat", "p.value")

  # cQTL SNP identifiers are colon or semicolon separated; rewrite to match the
  # underscore form used by the mQTL file
  parts <- strsplit(cqtl$SNP_raw, "[:;]")
  cqtl$SNP <- vapply(parts, function(x) paste(x[1], x[2], x[3], x[4], sep = "_"),
                     character(1))

  hits <- cqtl %>% filter(SNP %in% needed_snps)
  if (!SKIP_PRUNING) hits <- hits %>% filter(SNP %in% pruned_snp_ids)

  hits <- hits %>%
    select(SNP, gene, beta, t.stat, p.value) %>%
    distinct(SNP, .keep_all = TRUE)

  cat("  Outcome rows matched:", nrow(hits), "of", length(needed_snps), "needed SNPs\n")

  if (nrow(hits) == 0) {
    cat("  No outcome SNPs remain; exposure written but no GWAS file\n\n")
    return(invisible(NULL))
  }

  write.table(hits, out_gwas, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat("  Saved outcome:", basename(out_gwas), "\n\n")

  invisible(list(exposure = nrow(exposure_out), outcome = nrow(hits)))
}


# --- run ----------------------------------------------------------------------

n_done <- 0
n_skip <- 0

for (trait_label in trait_labels) {
  result <- tryCatch(build_one_cytokine(trait_label),
                     error = function(e) {
                       cat("  Error on", trait_label, "-", conditionMessage(e), "\n\n")
                       NULL
                     })
  if (!is.null(result)) n_done <- n_done + 1 else n_skip <- n_skip + 1
}

cat("Cytokines requested:", length(trait_labels), "\n")
cat("Built:              ", n_done, "\n")
cat("Skipped or failed:  ", n_skip, "\n")
cat("Pruning applied:    ", !SKIP_PRUNING, "\n")
cat("Files in:", out_dir, "\n")
