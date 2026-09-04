# =============================================================================
# Matched background CpG sets for the immune-trait enrichment analysis
#
# Generates 100 random CpG sets, each the same size as the cCpG list and each
# reproducing the cCpGs' chromosome by genomic feature composition. Sets are
# built stratum by stratum from the full array annotation rather than sampled
# uniformly at random, so the background carries the same genomic composition
# as the target.
#
# One set of 100 draws is generated here and reused by the enrichment script,
# so the EWAS Catalog and GWAS Catalog results rest on the same draws rather
# than on separate sources of sampling noise.
#
# The cCpGs themselves are removed from the pool, so a background set can never
# contain a real cCpG.
#
# Input
#   chrom_category_counts_sig.table.rds   strata table from the previous step
#   array_annotation.rds                  full array annotation, the pool to
#                                         sample from
#   cpgs_all_SI.txt                       the cCpGs, excluded from the pool
#
# Output
#   iteration_cpg_sets_LONG.tsv           one row per iteration and CpG
# =============================================================================

library(dplyr)
library(data.table)

data_dir <- "data"
out_dir  <- "results/EWAS_GWAS_enrichment"

STRATA_TABLE_PATH          <- file.path(out_dir,  "chrom_category_counts_sig.table.rds")
background_annotation_path <- file.path(data_dir, "array_annotation.rds")
CCPG_LIST_PATH             <- file.path(data_dir, "cpgs_all_SI.txt")

N_ITER <- 100
SEED   <- 42

set.seed(SEED)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strip <- function(x) sub("_[A-Z]{2}[0-9]{2}$", "", x)


# --- strata table ------------------------------------------------------------

sig <- readRDS(STRATA_TABLE_PATH)

target_total <- sig %>% distinct(seqnames, total_CpGs) %>% pull(total_CpGs) %>% sum()
cat("Strata table loaded. Target CpGs per iteration:", target_total, "\n")


# --- background pool ---------------------------------------------------------

back_ann <- as.data.frame(readRDS(background_annotation_path))
back_ann$name_clean <- strip(back_ann$name)

cat("Background annotation rows loaded:", nrow(back_ann), "\n")
cat("Unique background CpGs:", length(unique(back_ann$name_clean)), "\n")

ccpgs_clean <- unique(strip(read.table(CCPG_LIST_PATH, header = FALSE,
                                       stringsAsFactors = FALSE)$V1))

bg <- back_ann %>%
  distinct(seqnames, annot.type, name_clean) %>%
  filter(!name_clean %in% ccpgs_clean)

cat("Background pool after excluding cCpGs:", length(unique(bg$name_clean)),
    "CpGs across", nrow(bg), "strata rows\n")


# --- one matched draw --------------------------------------------------------

sample_one_iteration <- function(bg_pool, sig_table) {

  sampled_total <- data.frame()

  for (chr in unique(sig_table$seqnames)) {
    sig_chr <- sig_table %>% filter(seqnames == chr)

    for (annot in unique(sig_chr$annot.type)) {

      n_to_sample <- sig_chr %>% filter(annot.type == annot) %>% pull(n_CpGs)
      if (length(n_to_sample) == 0 || n_to_sample == 0) next

      bg_chr_annot <- bg_pool %>% filter(seqnames == chr, annot.type == annot)

      if (nrow(bg_chr_annot) < n_to_sample) {
        warning(paste("Not enough background CpGs for", chr, annot,
                      "- needed", n_to_sample, "found", nrow(bg_chr_annot)))
        next
      }

      sampled <- bg_chr_annot[sample(seq_len(nrow(bg_chr_annot)), n_to_sample,
                                     replace = FALSE), ]
      sampled_total <- rbind(sampled_total, sampled)
    }
  }

  sampled_total
}


# --- run the iterations ------------------------------------------------------

cat("Generating", N_ITER, "matched sets...\n")

all_iterations_long <- list()

for (iter in seq_len(N_ITER)) {

  sampled <- sample_one_iteration(bg, sig)

  if (is.null(sampled) || nrow(sampled) == 0) {
    cat("Iteration", iter, "skipped: no CpGs drawn\n")
    next
  }

  # A CpG belonging to several annotation strata is drawn more than once, so
  # the stratum draws sum to more than the target size
  if (nrow(sampled) > target_total) {
    sampled <- sampled[seq_len(target_total), ]
  }

  iter_cpgs <- unique(sampled$name_clean)

  all_iterations_long[[iter]] <- data.table(Iteration = iter,
                                            slot      = seq_along(iter_cpgs),
                                            cpg       = iter_cpgs)

  cat("Iteration", iter, "done -", length(iter_cpgs), "CpGs (target:", target_total, ")\n")
}

iterations_long <- rbindlist(all_iterations_long)

cat("Iterations generated:", length(unique(iterations_long$Iteration)), "\n")
print(iterations_long[, .N, by = Iteration][1:5])

fwrite(iterations_long, file.path(out_dir, "iteration_cpg_sets_LONG.tsv"), sep = "\t")

cat("\nWritten to ", file.path(out_dir, "iteration_cpg_sets_LONG.tsv"), "\n")
