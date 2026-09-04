# =============================================================================
# Chromosome by genomic feature composition of the cCpGs
#
# Builds the strata table that defines how the matched background sets are
# drawn: for every chromosome and annotation type, how many cCpGs fall in that
# combination. Each random background set must reproduce these counts.
#
# A CpG can carry several annotation rows in the source file, for example when
# overlapping transcripts annotate the same position. Counts are deduplicated
# per chromosome, annotation type and CpG before counting, so a CpG contributes
# at most one to any single stratum. It can still appear in several different
# strata, for example both a CpG shore and an enhancer, which is intended.
#
# Input
#   cpgs_all_SI.txt          the cCpG identifiers
#   ccpg_annotation.rds      annotation rows for those CpGs
#
# Output
#   chrom_category_counts_sig.table.rds        strata table for the sampling step
#   cCpG_chrom_feature_distribution_summary.tsv
#   cCpG_feature_type_totals.tsv
#   cCpG_chrom_feature_distribution_plot.pdf / .png
# =============================================================================

library(dplyr)
library(data.table)
library(ggplot2)

data_dir <- "data"
out_dir  <- "results/EWAS_GWAS_enrichment"

CCPG_LIST_PATH         <- file.path(data_dir, "cpgs_all_SI.txt")
ANNOTATION_SOURCE_PATH <- file.path(data_dir, "ccpg_annotation.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Annotation names carry a trailing array-design suffix that the CpG lists do not
strip <- function(x) sub("_[A-Z]{2}[0-9]{2}$", "", x)


# --- cCpGs and their annotation ----------------------------------------------

ccpgs_clean <- unique(strip(read.table(CCPG_LIST_PATH, header = FALSE,
                                       stringsAsFactors = FALSE)$V1))
cat("Target cCpGs:", length(ccpgs_clean), "\n")

data <- readRDS(ANNOTATION_SOURCE_PATH)
data$name_clean <- strip(data$name)
data_filtered <- data[data$name_clean %in% ccpgs_clean, ]

n_found <- length(unique(data_filtered$name_clean))
cat("Unique cCpGs found in the annotation source:", n_found, "\n")

if (n_found != length(ccpgs_clean)) {
  warning("Found ", n_found, " annotated cCpGs but expected ", length(ccpgs_clean),
          ". The remainder have no annotation row and cannot be matched.")
}


# --- strata table ------------------------------------------------------------

sig_table <- data_filtered %>%
  distinct(seqnames, annot.type, name_clean) %>%
  group_by(seqnames, annot.type) %>%
  summarise(n_CpGs = n(), .groups = "drop")

# Per chromosome totals count each CpG once, rather than summing across
# annotation types, since one CpG can belong to several feature classes
total_per_chr <- data_filtered %>%
  distinct(seqnames, name_clean) %>%
  group_by(seqnames) %>%
  summarise(total_CpGs = n(), .groups = "drop")

sig_table <- sig_table %>% left_join(total_per_chr, by = "seqnames")

cat("\nSum of unique per-chromosome CpG counts:", sum(total_per_chr$total_CpGs),
    "(expected", length(ccpgs_clean), ")\n")
cat("Chromosomes represented:", nrow(total_per_chr), "\n")
cat("Distinct annotation types:", length(unique(sig_table$annot.type)), "\n")

saveRDS(sig_table, file.path(out_dir, "chrom_category_counts_sig.table.rds"))

fwrite(sig_table %>% arrange(seqnames, desc(n_CpGs)),
       file.path(out_dir, "cCpG_chrom_feature_distribution_summary.tsv"), sep = "\t")

feature_totals <- sig_table %>%
  group_by(annot.type) %>%
  summarise(total_n_CpGs = sum(n_CpGs)) %>%
  arrange(desc(total_n_CpGs))

cat("\nPer-chromosome totals:\n")
print(total_per_chr %>% arrange(desc(total_CpGs)))

# This total exceeds the number of cCpGs because a CpG can belong to several
# feature types
cat("\nPer-feature-type totals:\n")
print(feature_totals)
fwrite(feature_totals, file.path(out_dir, "cCpG_feature_type_totals.tsv"), sep = "\t")


# --- composition heatmap -----------------------------------------------------

chr_order <- paste0("chr", c(1:22, "X", "Y"))
sig_table$seqnames <- factor(
  sig_table$seqnames,
  levels = rev(chr_order[chr_order %in% unique(sig_table$seqnames)])
)

p <- ggplot(sig_table, aes(x = annot.type, y = seqnames, fill = n_CpGs)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n_CpGs), size = 2.8, color = "black") +
  scale_fill_gradient(low = "#F0F4F8", high = "#2E6C7E", name = "cCpGs") +
  labs(x = "Genomic feature type",
       y = "Chromosome",
       title = "Chromosome by genomic feature distribution of cytokine-associated CpGs") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(out_dir, "cCpG_chrom_feature_distribution_plot.pdf"), p, width = 10, height = 8)
ggsave(file.path(out_dir, "cCpG_chrom_feature_distribution_plot.png"), p, width = 10, height = 8, dpi = 300)

cat("\nWritten to ", out_dir, "\n")
