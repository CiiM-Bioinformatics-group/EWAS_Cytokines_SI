#!/usr/bin/env Rscript
# =============================================================================
# Validation of the SI EWAS findings in the BCG-PRIME cohort
#
# Every CpG that reached FDR significance for any cytokine in the discovery
# cohort is looked up in the BCG-PRIME results. A CpG counts as consistent when
# BCG-PRIME shows an effect in the same direction with p < 0.05 for at least one
# cytokine.
#
# Outputs: the result tables as text, CSV and an annotated Excel workbook, and
# the signed p-value scatter plot used as the replication figure.
#
# Usage:
#   Rscript Validaiton_EWAS_results_BCG_PRIME.R
#
# Expected inputs
#   discovery_dir    one file per cytokine, whitespace separated, no header,
#                    first three columns: CpG, effect size, p-value
#   replication_dir  one file per cytokine, tab separated, with a header
#                    containing CpGsite, Estimate and Pr(>|z|)
# =============================================================================

library(data.table)
library(ggplot2)
library(openxlsx)


# --- settings ----------------------------------------------------------------

discovery_dir   <- "results/EWAS_discovery_SI"
replication_dir <- "results/EWAS_replication_BCG_PRIME"
out_dir         <- "results/validation_SI_in_BCG_PRIME"

# Patterns identifying the result files, also used to derive the cytokine label
discovery_pattern   <- "\\.sig\\.FDR\\.txt$"
replication_pattern <- "resultsmodel1\\.txt$"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", out_dir, "\n\n")


# =============================================================================
# STEP 1 — Discovery hits, one row per CpG
# A CpG can be significant for several cytokines. The row kept is the one with
# the lowest discovery p-value.
# =============================================================================
cat("STEP 1: Loading discovery FDR files...\n")

si_files <- list.files(discovery_dir, pattern = discovery_pattern, full.names = TRUE)
if (length(si_files) == 0) stop("No discovery result files found in: ", discovery_dir)

si_all <- rbindlist(lapply(si_files, function(f) {
  if (file.size(f) == 0) { message("  Skipping empty: ", basename(f)); return(NULL) }
  dt <- fread(f, header = FALSE, sep = " ")
  dt <- dt[, 1:3]
  setnames(dt, c("CpG", "SI_Effect", "SI_Pvalue"))
  dt[, SI_file := sub(discovery_pattern, "", basename(f))]
  dt
}), fill = TRUE)

cat(sprintf("  Total SI entries (all files) : %s\n", format(nrow(si_all),        big.mark = ",")))
cat(sprintf("  Unique CpG sites             : %s\n", format(uniqueN(si_all$CpG), big.mark = ",")))

si_unique <- si_all[order(SI_Pvalue)][!duplicated(CpG)]
cat(sprintf("  After deduplication          : %s\n\n", format(nrow(si_unique), big.mark = ",")))


# =============================================================================
# STEP 2 — Replication results
# =============================================================================
cat("STEP 2: Loading BCG-PRIME files...\n")

prime_files <- list.files(replication_dir, pattern = replication_pattern, full.names = TRUE)
if (length(prime_files) == 0) stop("No replication result files found in: ", replication_dir)

n_prime_files <- length(prime_files)

prime_all <- rbindlist(lapply(prime_files, function(f) {
  dt <- fread(f, header = TRUE, sep = "\t")
  dt <- dt[, .(CpG            = CpGsite,
               Prime_Estimate = Estimate,
               Prime_Pvalue   = `Pr(>|z|)`)]
  dt[, Prime_file := sub(replication_pattern, "", basename(f))]
  dt
}))

cat(sprintf("  Replication files loaded: %s\n", n_prime_files))
cat(sprintf("  Total rows loaded       : %s\n\n", format(nrow(prime_all), big.mark = ",")))


# =============================================================================
# STEP 3 — Best replication hit per CpG
#   A: best same-direction hit with p<0.05  -> Consistent = TRUE
#   B: otherwise the best hit in any direction -> Consistent = FALSE
#
# The p-value kept is a minimum taken over all replication files, so it is a
# selected value rather than a single test. N_Prime_SameDir_p05 records how many
# files qualified, so the selection can be accounted for when reporting.
# =============================================================================
cat("STEP 3: Extracting best replication hit per CpG...\n")

merged <- merge(si_unique[, .(CpG, SI_Effect, SI_Pvalue, SI_file)],
                prime_all, by = "CpG", allow.cartesian = TRUE)
merged[, same_dir := sign(SI_Effect) == sign(Prime_Estimate)]

best_same <- merged[same_dir == TRUE & Prime_Pvalue < 0.05][order(Prime_Pvalue),
  .(Prime_Estimate      = Prime_Estimate[1],
    Prime_Pvalue        = Prime_Pvalue[1],
    Prime_file          = Prime_file[1],
    N_Prime_SameDir_p05 = .N,
    Consistent          = TRUE),
  by = CpG]

cpgs_no_rep <- setdiff(si_unique$CpG, best_same$CpG)
best_any <- merged[CpG %in% cpgs_no_rep][order(Prime_Pvalue),
  .(Prime_Estimate      = Prime_Estimate[1],
    Prime_Pvalue        = Prime_Pvalue[1],
    Prime_file          = Prime_file[1],
    N_Prime_SameDir_p05 = 0L,
    Consistent          = FALSE),
  by = CpG]

results <- merge(si_unique, rbind(best_same, best_any), by = "CpG")

cat(sprintf("  Consistent (same dir + p<0.05)  : %s\n", format(sum(results$Consistent),  big.mark = ",")))
cat(sprintf("  Non-consistent                  : %s\n", format(sum(!results$Consistent), big.mark = ",")))
cat(sprintf("  Total CpGs                      : %s\n\n", format(nrow(results),          big.mark = ",")))


# =============================================================================
# STEP 4 — BH correction and plot coordinates
# =============================================================================
cat("STEP 4: BH FDR correction across", nrow(results), "CpGs...\n")

results[, FDR_BH := p.adjust(Prime_Pvalue, method = "BH")]

# Signed axis values: significance on the scale, direction of effect in the sign
results[, SI_signed    := -log10(SI_Pvalue)    * sign(SI_Effect)]
results[, Prime_signed := -log10(Prime_Pvalue) * sign(Prime_Estimate)]

# The threshold is read off the data: the largest p-value still passing FDR<0.05
fdr_sig         <- results[FDR_BH < 0.05]
fdr_cutoff_pval <- if (nrow(fdr_sig) > 0) max(fdr_sig$Prime_Pvalue) else NA_real_
fdr_line_pos    <- if (!is.na(fdr_cutoff_pval)) -log10(fdr_cutoff_pval) else NA_real_
fdr_line_neg    <- if (!is.na(fdr_cutoff_pval))  log10(fdr_cutoff_pval) else NA_real_

n_p05 <- sum(results$Consistent)
n_fdr <- sum(results$Consistent & results$FDR_BH < 0.05)
n_nc  <- sum(!results$Consistent)

cat(sprintf("  FDR<0.05 p-value cutoff         : %.4e\n", fdr_cutoff_pval))
cat(sprintf("  -log10 threshold                : %.2f\n",  fdr_line_pos))
cat(sprintf("  Consistent p<0.05               : %s\n",    format(n_p05, big.mark = ",")))
cat(sprintf("  Consistent FDR<0.05             : %s\n\n",  format(n_fdr, big.mark = ",")))

lab_nc  <- sprintf("SI, Non-consistent CpGs  (n = %s)",                 format(n_nc,  big.mark = ","))
lab_con <- sprintf("SI, Consistent CpGs (p < 0.05 in Prime)  (n = %s)", format(n_p05, big.mark = ","))

results[, Category := ifelse(Consistent, lab_con, lab_nc)]
results[, Category := factor(Category, levels = c(lab_nc, lab_con))]


# =============================================================================
# STEP 5 — Text and CSV output
# =============================================================================
cat("STEP 5: Saving text/CSV files...\n")

col_order <- c("CpG", "SI_file", "SI_Effect", "SI_Pvalue",
               "Prime_file", "Prime_Estimate", "Prime_Pvalue", "FDR_BH",
               "N_Prime_SameDir_p05", "Consistent", "Category")

all_cpgs        <- results[order(Prime_Pvalue), ..col_order]
consistent_p05  <- results[Consistent == TRUE][order(Prime_Pvalue), ..col_order]
consistent_fdr  <- results[Consistent == TRUE & FDR_BH < 0.05][order(FDR_BH), ..col_order]

fwrite(all_cpgs,       file.path(out_dir, "all_cpgs_with_prime_info.txt"),     sep = "\t", quote = FALSE, na = "NA")
fwrite(all_cpgs,       file.path(out_dir, "all_cpgs_with_prime_info.csv"),                 quote = FALSE, na = "NA")
fwrite(consistent_p05, file.path(out_dir, "replicated_samedirection_p05.txt"),  sep = "\t", quote = FALSE, na = "NA")
fwrite(consistent_p05, file.path(out_dir, "replicated_samedirection_p05.csv"),              quote = FALSE, na = "NA")
fwrite(consistent_fdr, file.path(out_dir, "replicated_samedirection_FDR05.txt"), sep = "\t", quote = FALSE, na = "NA")
fwrite(consistent_fdr, file.path(out_dir, "replicated_samedirection_FDR05.csv"),             quote = FALSE, na = "NA")

cat("  Text and CSV files saved.\n\n")


# =============================================================================
# STEP 6 — Excel workbook with a legend describing every column
# =============================================================================
cat("STEP 6: Building Excel workbook...\n")

wb <- createWorkbook()

write_sheet <- function(wb, sheet_name, dt) {
  addWorksheet(wb, sheet_name)
  hs <- createStyle(fontColour = "#FFFFFF", fgFill = "#2E6B3E",
                    halign = "CENTER", textDecoration = "Bold",
                    border = "Bottom", wrapText = TRUE)
  writeDataTable(wb, sheet_name, as.data.frame(dt),
                 tableStyle = "TableStyleMedium7", withFilter = TRUE)
  addStyle(wb, sheet_name, hs, rows = 1, cols = seq_len(ncol(dt)), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(dt)), widths = 20)
}

addWorksheet(wb, "Legend")
legend_df <- data.frame(
  Sheet = c(
    rep("All_CpGs", 11), "---",
    "Consistent_p05", "Consistent_FDR05", "---",
    "Plot", "Plot", "Plot"
  ),
  Column_or_Item = c(
    "CpG", "SI_file", "SI_Effect", "SI_Pvalue",
    "Prime_file", "Prime_Estimate", "Prime_Pvalue", "FDR_BH",
    "N_Prime_SameDir_p05", "Consistent", "Category", "---",
    "Content", "Content", "---",
    "Grey circles", "Navy triangles", "Amber horizontal line"
  ),
  Description = c(
    "CpG site identifier",
    "Discovery cytokine of origin; row kept = lowest discovery p-value across all files",
    "Effect size (beta) from the discovery EWAS",
    "Raw p-value from the discovery EWAS (FDR-significant in discovery)",
    "Replication cytokine with the lowest p-value for this CpG, same direction where possible",
    "Effect size (beta) from the replication EWAS, best matching file",
    "Raw p-value from the replication EWAS, best matching file",
    paste0("BH-adjusted replication p-value, corrected across all ", nrow(results), " unique CpGs"),
    paste0("Number of replication files (out of ", n_prime_files,
           ") with same-direction effect and p<0.05 for this CpG"),
    "TRUE = at least one same-direction hit with p<0.05; FALSE = no such hit",
    "Plot label: Consistent or Non-consistent, with counts",
    "---",
    paste0("All consistent CpGs (same direction + p<0.05 in any replication file), n=",
           n_p05, ", ordered by replication p-value"),
    paste0("Consistent CpGs with FDR_BH < 0.05, n=", n_fdr, ", ordered by FDR"),
    "---",
    paste0("Non-consistent CpGs: opposite direction or no p<0.05 replication hit (n=", n_nc, ")"),
    paste0("Consistent CpGs: same direction + p<0.05 in at least one replication file (n=", n_p05, ")"),
    paste0("FDR=0.05 threshold on the replication -log10(p) axis (p=",
           formatC(fdr_cutoff_pval, format = "e", digits = 2),
           "); CpGs beyond this line have FDR<0.05 (n=", n_fdr, ")")
  ),
  stringsAsFactors = FALSE
)
writeDataTable(wb, "Legend", legend_df, tableStyle = "TableStyleMedium7", withFilter = FALSE)
hs_leg <- createStyle(fontColour = "#FFFFFF", fgFill = "#1A4A6B", textDecoration = "Bold", halign = "CENTER")
addStyle(wb, "Legend", hs_leg, rows = 1, cols = 1:3, gridExpand = TRUE)
setColWidths(wb, "Legend", cols = 1:3, widths = c(22, 28, 90))

addWorksheet(wb, "Summary")
summary_df <- data.frame(
  Metric = c(
    "Total discovery entries loaded (all files)",
    "Unique CpG sites in discovery (best p kept)",
    "Replication files checked",
    "BH FDR correction: number of tests",
    "Non-consistent CpGs",
    "Consistent CpGs (same dir + p<0.05 in replication)",
    "Consistent CpGs with FDR < 0.05 (BH)",
    "FDR=0.05 replication p-value cutoff",
    "FDR=0.05 -log10(p) threshold on plot"
  ),
  Value = c(
    nrow(si_all), nrow(si_unique), n_prime_files, nrow(results),
    n_nc, n_p05, n_fdr,
    formatC(fdr_cutoff_pval, format = "e", digits = 3),
    round(fdr_line_pos, 3)
  ),
  stringsAsFactors = FALSE
)
writeDataTable(wb, "Summary", summary_df, tableStyle = "TableStyleMedium7", withFilter = FALSE)
hs_sum <- createStyle(fontColour = "#FFFFFF", fgFill = "#2E6B3E", textDecoration = "Bold", halign = "CENTER")
addStyle(wb, "Summary", hs_sum, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Summary", cols = 1:2, widths = c(50, 20))

write_sheet(wb, "All_CpGs",         all_cpgs)
write_sheet(wb, "Consistent_p05",   consistent_p05)
write_sheet(wb, "Consistent_FDR05", consistent_fdr)

xl_file <- file.path(out_dir, "replication_results.xlsx")
saveWorkbook(wb, xl_file, overwrite = TRUE)
cat(sprintf("  Excel --> %s\n\n", xl_file))


# =============================================================================
# STEP 7 — Signed p-value scatter
# Both axes carry significance in the magnitude and direction of effect in the
# sign, so agreement between cohorts falls in the lower-left and upper-right
# quadrants.
# =============================================================================
cat("STEP 7: Drawing replication plot...\n")

col_nc  <- "#adb5bd"   # grey  - non-replicated background
col_con <- "#1a5276"   # navy  - replicated signal
fdr_col <- "#c9a227"   # amber - FDR threshold

plot_nc  <- results[Consistent == FALSE]
plot_con <- results[Consistent == TRUE]

p <- ggplot() +

  # Non-consistent drawn first so the replicated points sit on top
  geom_point(data = plot_nc,
             aes(x = SI_signed, y = Prime_signed,
                 colour = Category, shape = Category),
             size = 1.6, alpha = 0.35) +

  geom_point(data = plot_con,
             aes(x = SI_signed, y = Prime_signed,
                 colour = Category, shape = Category),
             size = 2.6, alpha = 0.90) +

  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.4) +

  # Mirrored, because the axis is signed
  geom_hline(yintercept = fdr_line_pos, linetype = "longdash",
             colour = fdr_col, linewidth = 0.85) +
  geom_hline(yintercept = fdr_line_neg, linetype = "longdash",
             colour = fdr_col, linewidth = 0.85) +

  annotate("text",
           x      = max(abs(results$SI_signed), na.rm = TRUE) * 0.97,
           y      = fdr_line_pos + 0.22,
           label  = sprintf("FDR = 0.05  (n = %s)", format(n_fdr, big.mark = ",")),
           colour = fdr_col, size = 4, hjust = 1, fontface = "italic") +

  scale_colour_manual(values = setNames(c(col_nc, col_con), c(lab_nc, lab_con))) +
  scale_shape_manual(values  = setNames(c(16L, 17L),        c(lab_nc, lab_con))) +

  labs(
    title    = "Replication of SI CpGs in BCG-PRIME cohort",
    subtitle = "All cytokines - consistent: same direction, p < 0.05",
    x        = expression(-log[10](p) %.% sign(beta) ~~ "(Discovery: SI)"),
    y        = expression(-log[10](p) %.% sign(beta) ~~ "(Replication: BCG-PRIME)"),
    colour   = NULL,
    shape    = NULL
  ) +

  theme_classic(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 12, hjust = 0.5, colour = "grey40",
                                    margin = margin(b = 10)),
    axis.title       = element_text(size = 13, face = "bold"),
    axis.text        = element_text(size = 12, colour = "black"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key       = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    plot.margin      = margin(12, 16, 8, 10)
  ) +
  guides(
    colour = guide_legend(ncol = 1, override.aes = list(size = 4.5, alpha = 1)),
    shape  = guide_legend(ncol = 1, override.aes = list(size = 4.5, alpha = 1))
  )

pdf_file <- file.path(out_dir, "replication_SI_in_BCG_PRIME.pdf")
png_file <- file.path(out_dir, "replication_SI_in_BCG_PRIME.png")
ggsave(pdf_file, plot = p, width = 7, height = 7.5, device = "pdf")
ggsave(png_file, plot = p, width = 7, height = 7.5, dpi = 300)
cat(sprintf("  PDF --> %s\n", pdf_file))
cat(sprintf("  PNG --> %s\n\n", png_file))


# =============================================================================
# SUMMARY
# =============================================================================
cat(strrep("=", 65), "\n")
cat("SUMMARY\n")
cat(strrep("=", 65), "\n")
cat(sprintf("Discovery entries loaded                   : %s\n", format(nrow(si_all),    big.mark = ",")))
cat(sprintf("Unique CpG sites                           : %s\n", format(nrow(si_unique), big.mark = ",")))
cat(sprintf("Replication files checked                  : %s\n", n_prime_files))
cat(sprintf("BH FDR n tests                             : %s\n", format(nrow(results),   big.mark = ",")))
cat(sprintf("Non-consistent CpGs                        : %s\n", format(n_nc,            big.mark = ",")))
cat(sprintf("Consistent (same dir + p<0.05)             : %s\n", format(n_p05,           big.mark = ",")))
cat(sprintf("  of which FDR < 0.05 (BH)                 : %s\n", format(n_fdr,           big.mark = ",")))
cat(sprintf("FDR=0.05 threshold on -log10(p) axis       : %.3f\n", fdr_line_pos))
cat(strrep("=", 65), "\n")
cat(sprintf("All outputs in: %s\n", out_dir))
