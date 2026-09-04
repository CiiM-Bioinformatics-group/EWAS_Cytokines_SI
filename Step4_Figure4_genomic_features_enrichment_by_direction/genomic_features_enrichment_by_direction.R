# =============================================================================
# Genomic feature annotation of cCpGs, and enrichment by direction of effect
#
# cCpGs are annotated to CpG context categories (islands, shores, shelves, open
# sea) and to functional elements (promoters, FANTOM5 enhancers, introns,
# exons, UTRs, CDS, intergenic regions, lncRNA loci and exon-intron
# boundaries), fifteen categories in total. The annotation itself is produced
# with the annotatr package on the hg38 build and is read in here.
#
# Each cCpG is classified by the direction of its association with cytokine
# responses across stimuli:
#   Positive-Beta   positive in every association
#   Negative-Beta   negative in every association
#   Mixed           both directions present, reported but not tested
#
# Each group is compared independently against the non-significant array
# background, meaning all array CpGs except the cCpGs. Because that background
# holds hundreds of thousands of CpGs, a Fisher test would be significant on
# negligible differences, so enrichment is assessed by permutation instead:
# sets of the same size as the group are drawn from the background and the
# percentage in each annotation category is recorded, giving an empirical p
# value for the observed percentage.
#
# Input
#   SI.txt                          the cCpG identifiers
#   all.merged.standardized.tsv     EWAS results: CpG, condition, beta, p
#   si_annotated_df.rds             annotatr output for the cCpGs
#   array_annotated.rds             annotatr output for the full array
#
# Output
#   SI_direction_pie.png
#   SI_cpg_context_{all,positive,negative}.png
#   SI_features_{all,positive,negative}.png
#   SI_annotation_distribution.png
#   SI_{all,pos,neg}_vs_nonsig_permutation_table.csv
#   SI_{all,pos,neg}_vs_nonsig_permutation.png
#   SI_pos_neg_diverging_enrichment.png
# =============================================================================

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)


# --- settings ----------------------------------------------------------------

DATA_DIR <- "data/Figure4"
OUT_DIR  <- "results/Figure4"

N_PERM <- 1000
SEED   <- 42

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# --- cCpGs and direction of effect -------------------------------------------

si_cpgs <- readLines(file.path(DATA_DIR, "SI.txt"))
cat("Total cCpGs:", length(si_cpgs), "\n")

df <- fread(file.path(DATA_DIR, "all.merged.standardized.tsv"), header = FALSE)
colnames(df) <- c("CpG", "Condition", "Beta", "Pvalue")

df_si <- df[CpG %in% si_cpgs]
cat("EWAS rows for cCpGs:", nrow(df_si), "\n")

cpg_dir <- df_si %>%
  group_by(CpG) %>%
  summarise(
    direction = case_when(
      all(Beta > 0) ~ "Always Positive",
      all(Beta < 0) ~ "Always Negative",
      TRUE          ~ "Mixed"
    ),
    .groups = "drop"
  )

cat("\nDirection distribution:\n")
print(table(cpg_dir$direction))

cpgs_pos <- cpg_dir$CpG[cpg_dir$direction == "Always Positive"]
cpgs_neg <- cpg_dir$CpG[cpg_dir$direction == "Always Negative"]
n_pos <- length(cpgs_pos)
n_neg <- length(cpgs_neg)
cat("Positive-Beta cCpGs:", n_pos, "\n")
cat("Negative-Beta cCpGs:", n_neg, "\n\n")


# --- annotation --------------------------------------------------------------

si_ann <- readRDS(file.path(DATA_DIR, "si_annotated_df.rds"))
cat("Annotation rows:", nrow(si_ann), " | unique cCpGs:", length(unique(si_ann$name)), "\n\n")

si_ann <- si_ann %>%
  left_join(cpg_dir %>% select(CpG, direction), by = c("name" = "CpG")) %>%
  mutate(direction = ifelse(is.na(direction), "Mixed", direction))

# The fifteen categories the enrichment is corrected across
annotation_categories <- c(
  "hg38_genes_promoters", "hg38_enhancers_fantom",
  "hg38_cpg_islands", "hg38_cpg_shores", "hg38_cpg_shelves", "hg38_cpg_inter",
  "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_5UTRs",
  "hg38_genes_3UTRs", "hg38_genes_cds", "hg38_genes_intergenic",
  "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries",
  "hg38_lncrna_gencode"
)
cat("Annotation categories present:", length(unique(si_ann$annot.type)),
    "of", length(annotation_categories), "expected\n\n")


# --- palettes ----------------------------------------------------------------

cpg_context_types <- c("hg38_cpg_islands", "hg38_cpg_shores",
                       "hg38_cpg_shelves", "hg38_cpg_inter")

cpg_colors <- c("Islands" = "#335C67", "Shores" = "#E09F3E",
                "Shelves" = "#BB9457", "Open Sea" = "#D36582")

feature_colors <- c(
  "promoters"            = "#1f77b4",
  "enhancers_fantom"     = "#ff7f0e",
  "exons"                = "#2ca02c",
  "introns"              = "#9467bd",
  "5UTRs"                = "#8c564b",
  "3UTRs"                = "#e377c2",
  "intergenic"           = "#17becf",
  "intronexonboundaries" = "#7f7f7f",
  "exonintronboundaries" = "#bcbd22",
  "lncrna_gencode"       = "#d62728",
  "cds"                  = "#aec7e8"
)

dir_colors <- c("Always Negative" = "#344E41",
                "Always Positive" = "#BB9457",
                "Mixed"           = "#A3B18A")

green_palette <- c("#DAD7CD", "#A3B18A", "#588157", "#3A5A40", "#344E41")


# --- direction pie chart -----------------------------------------------------

dir_counts <- as.data.frame(table(cpg_dir$direction))
colnames(dir_counts) <- c("Direction", "Count")
dir_counts$Direction <- factor(dir_counts$Direction,
                               levels = c("Always Negative", "Always Positive", "Mixed"))
dir_counts$Pct   <- round(dir_counts$Count / sum(dir_counts$Count) * 100, 1)
dir_counts$Label <- paste0(dir_counts$Pct, "%")

pie_dir <- ggplot(dir_counts, aes(x = "", y = Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.8) +
  coord_polar("y") +
  scale_fill_manual(values = dir_colors) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
            size = 6, fontface = "bold", color = "black") +
  labs(title = "cCpGs: Methylation Direction", fill = NULL) +
  theme_void(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
    legend.text     = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1.4, "lines"),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(OUT_DIR, "SI_direction_pie.png"), pie_dir, width = 7, height = 6, dpi = 300)
cat("Saved: SI_direction_pie.png\n")


# --- CpG context composition -------------------------------------------------

make_cpg_pie <- function(ann_df, total_n, title, out_file) {

  summary_df <- ann_df %>%
    filter(annot.type %in% cpg_context_types) %>%
    group_by(annot.type) %>%
    summarise(Count = n_distinct(name), .groups = "drop") %>%
    mutate(
      Percentage = round(Count / total_n * 100, 1),
      Category = case_when(
        annot.type == "hg38_cpg_islands" ~ "Islands",
        annot.type == "hg38_cpg_shores"  ~ "Shores",
        annot.type == "hg38_cpg_shelves" ~ "Shelves",
        annot.type == "hg38_cpg_inter"   ~ "Open Sea"
      ),
      Label = paste0(Category, "\n", Percentage, "%")
    )

  colors_use <- cpg_colors[summary_df$Category]
  names(colors_use) <- summary_df$Label

  p <- ggplot(summary_df, aes(x = "", y = Percentage, fill = Label)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.6) +
    coord_polar("y") +
    scale_fill_manual(values = colors_use) +
    labs(title = title, fill = "CpG Context") +
    theme_void(base_size = 16) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
      legend.title    = element_text(face = "bold", size = 14),
      legend.text     = element_text(size = 13),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(out_file, p, width = 7, height = 6, dpi = 300)
  cat("Saved:", basename(out_file), "\n")
}

make_cpg_pie(si_ann, length(si_cpgs),
             paste0("CpG Context: All cCpGs (n=", length(si_cpgs), ")"),
             file.path(OUT_DIR, "SI_cpg_context_all.png"))

make_cpg_pie(si_ann %>% filter(direction == "Always Positive"), n_pos,
             paste0("CpG Context: Positive-Beta cCpGs (n=", n_pos, ")"),
             file.path(OUT_DIR, "SI_cpg_context_positive.png"))

make_cpg_pie(si_ann %>% filter(direction == "Always Negative"), n_neg,
             paste0("CpG Context: Negative-Beta cCpGs (n=", n_neg, ")"),
             file.path(OUT_DIR, "SI_cpg_context_negative.png"))


# --- genomic feature composition ---------------------------------------------

make_feature_bar <- function(ann_df, total_n, title, out_file) {

  df_plot <- ann_df %>%
    filter(!annot.type %in% cpg_context_types) %>%
    group_by(annot.type) %>%
    summarise(Count = n_distinct(name), .groups = "drop") %>%
    mutate(
      Pct       = round(Count / total_n * 100, 1),
      Clean     = gsub("^hg38_genes_|^hg38_", "", annot.type),
      Label_txt = paste0(Count, " (", Pct, "%)")
    ) %>%
    arrange(desc(Pct))

  fill_vals <- feature_colors[df_plot$Clean]
  fill_vals[is.na(fill_vals)] <- "grey70"
  names(fill_vals) <- df_plot$Clean

  p <- ggplot(df_plot, aes(x = reorder(Clean, Pct), y = Pct, fill = Clean)) +
    geom_bar(stat = "identity", width = 0.78, color = "black") +
    geom_text(aes(label = Label_txt), hjust = -0.1, size = 5, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = fill_vals) +
    labs(title = title, x = NULL, y = "% of CpGs") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 17),
      axis.text       = element_text(size = 13, face = "bold"),
      axis.title.x    = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin     = margin(10, 90, 10, 10)
    ) +
    ylim(0, max(df_plot$Pct) * 1.45)

  ggsave(out_file, p, width = 10, height = 9, dpi = 300)
  cat("Saved:", basename(out_file), "\n")
}

make_feature_bar(si_ann, length(si_cpgs), "Genomic Features: All cCpGs",
                 file.path(OUT_DIR, "SI_features_all.png"))

make_feature_bar(si_ann %>% filter(direction == "Always Positive"), n_pos,
                 "Genomic Features: Positive-Beta cCpGs",
                 file.path(OUT_DIR, "SI_features_positive.png"))

make_feature_bar(si_ann %>% filter(direction == "Always Negative"), n_neg,
                 "Genomic Features: Negative-Beta cCpGs",
                 file.path(OUT_DIR, "SI_features_negative.png"))


# --- annotation distribution across all categories ---------------------------

ann_for_plot <- si_ann %>%
  mutate(
    Category = gsub("hg38_genes_|hg38_", "", annot.type),
    Category = ifelse(Category == "cpg_inter", "Open Sea", Category),
    Class    = ifelse(annot.type %in% cpg_context_types, "CpG Context", "Genomic Feature")
  ) %>%
  group_by(Category, Class) %>%
  summarise(Count = n_distinct(name), .groups = "drop") %>%
  mutate(Pct = round(Count / length(si_cpgs) * 100, 1)) %>%
  arrange(Class, desc(Pct)) %>%
  mutate(Category = factor(Category, levels = rev(Category)))

class_colors <- c("CpG Context" = "#588157", "Genomic Feature" = "#BB9457")

pA <- ggplot(ann_for_plot, aes(x = Category, y = Pct, fill = Class)) +
  geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.35) +
  geom_text(aes(label = paste0(Pct, "%")), hjust = -0.1, size = 4.5, fontface = "bold") +
  coord_flip(ylim = c(0, max(ann_for_plot$Pct) * 1.38)) +
  scale_fill_manual(values = class_colors) +
  labs(x = NULL, y = "% of cCpGs", fill = "Annotation Class",
       title = "Annotation Distribution: cCpGs") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y        = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 13, face = "bold"),
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.margin        = margin(10, 80, 10, 10)
  )

ggsave(file.path(OUT_DIR, "SI_annotation_distribution.png"), pA, width = 9, height = 10, dpi = 300)
cat("Saved: SI_annotation_distribution.png\n")


# --- non-significant background ----------------------------------------------

bg_cache <- file.path(DATA_DIR, "bg_nonsig_annotated.rds")

if (!file.exists(bg_cache)) {
  suppressMessages(library(GenomicRanges))
  d_array   <- readRDS(file.path(DATA_DIR, "array_annotated.rds"))
  ann_array <- mcols(d_array)$annot
  df_array  <- data.frame(
    name         = mcols(d_array)$name,
    annot.type   = mcols(ann_array)$type,
    annot.symbol = mcols(ann_array)$symbol,
    stringsAsFactors = FALSE
  )
  bg_nonsig <- df_array %>% filter(!name %in% si_cpgs)
  saveRDS(bg_nonsig, bg_cache)
  cat("Built and cached non-significant background:", nrow(bg_nonsig), "rows\n")
} else {
  bg_nonsig <- readRDS(bg_cache)
  cat("Loaded non-significant background:", nrow(bg_nonsig), "rows\n")
}

bg_names <- unique(bg_nonsig$name)
n_bg     <- length(bg_names)
n_bg_ns  <- n_bg
cat("Non-significant background CpGs:", n_bg, "\n\n")

bg_ann_counts <- bg_nonsig %>%
  group_by(annot.type) %>%
  summarise(BG_n = n_distinct(name), .groups = "drop") %>%
  mutate(BG_pct = BG_n / n_bg * 100)


# --- permutation enrichment --------------------------------------------------

compute_perm_res <- function(target_ann) {

  n_target   <- length(unique(target_ann$name))
  categories <- unique(target_ann$annot.type)

  obs <- target_ann %>%
    group_by(annot.type) %>%
    summarise(SI_n = n_distinct(name), .groups = "drop") %>%
    mutate(SI_pct = SI_n / n_target * 100)

  perm_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(categories),
                     dimnames = list(NULL, categories))

  set.seed(SEED)
  for (i in seq_len(N_PERM)) {
    samp <- sample(bg_names, n_target, replace = FALSE)
    samp_a <- bg_nonsig %>%
      filter(name %in% samp) %>%
      group_by(annot.type) %>%
      summarise(nn = n_distinct(name), .groups = "drop")
    for (ct in categories) {
      cnt <- samp_a$nn[samp_a$annot.type == ct]
      perm_mat[i, ct] <- ifelse(length(cnt) == 0, 0, cnt / n_target * 100)
    }
  }

  # Tail probability on the side the observation falls
  emp_pval <- sapply(categories, function(ct) {
    ob   <- obs$SI_pct[obs$annot.type == ct]
    perm <- perm_mat[, ct]
    if (ob >= mean(perm)) mean(perm >= ob) else mean(perm <= ob)
  })

  obs %>%
    left_join(bg_ann_counts, by = "annot.type") %>%
    mutate(
      OR        = round((SI_n / (n_target - SI_n)) / (BG_n / (n_bg - BG_n)), 3),
      perm_pval = emp_pval[annot.type],
      FDR_BH    = p.adjust(perm_pval, method = "BH"),
      Category  = gsub("hg38_genes_|hg38_", "", annot.type),
      Category  = ifelse(Category == "cpg_inter", "Open Sea", Category),
      sig       = ifelse(perm_pval < 0.05 & OR > 1, "*",
                         ifelse(perm_pval < 0.05 & OR < 1, "**", "")),
      label     = paste0(OR, ifelse(sig != "", " ", ""), sig),
      nlog10p   = -log10(perm_pval + 1e-4),
      n_target  = n_target
    ) %>%
    arrange(desc(OR))
}

cat("Permutation: all cCpGs\n")
res_all_perm <- compute_perm_res(si_ann)
cat("Permutation: Positive-Beta cCpGs\n")
res_pos_perm <- compute_perm_res(si_ann %>% filter(direction == "Always Positive"))
cat("Permutation: Negative-Beta cCpGs\n")
res_neg_perm <- compute_perm_res(si_ann %>% filter(direction == "Always Negative"))

# Shared colour scale so the positive and negative panels are comparable
global_max_nlog <- ceiling(max(c(res_pos_perm$nlog10p, res_neg_perm$nlog10p)) * 10) / 10
global_breaks   <- pretty(c(0, global_max_nlog), n = 4)
cat("Shared colour scale maximum:", global_max_nlog, "\n\n")

save_perm_table <- function(res, out_prefix) {
  tbl <- res %>%
    select(Category, CpGs_in_group = SI_n, Group_pct = SI_pct,
           Background_CpGs = BG_n, Background_pct = BG_pct,
           OR, perm_pvalue = perm_pval, FDR_BH, Significance = sig) %>%
    arrange(perm_pvalue)
  write.csv(tbl, file.path(OUT_DIR, paste0(out_prefix, "_permutation_table.csv")),
            row.names = FALSE)
  cat("Saved table:", out_prefix, "\n")
}

save_perm_table(res_all_perm, "SI_all_vs_nonsig")
save_perm_table(res_pos_perm, "SI_pos_vs_nonsig")
save_perm_table(res_neg_perm, "SI_neg_vs_nonsig")


# --- per-group enrichment panels ---------------------------------------------

make_perm_plot <- function(res, group_label, out_file, scale_max, scale_breaks) {

  n_target     <- res$n_target[1]
  res$Category <- factor(res$Category, levels = res$Category[order(res$OR)])
  y_max        <- max(res$OR) * 1.55

  p <- ggplot(res, aes(x = Category, y = OR, fill = nlog10p)) +
    geom_bar(stat = "identity", width = 0.72, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1.1) +
    geom_text(aes(label = label), hjust = -0.12, size = 6.5, fontface = "bold") +
    coord_flip(ylim = c(0, y_max)) +
    scale_fill_gradientn(colors = green_palette,
                         name = expression(-log[10](p[perm])),
                         limits = c(0, scale_max), breaks = scale_breaks) +
    annotate("text", x = 3.4, y = y_max * 0.985, label = "*   Enriched vs background",
             hjust = 1, size = 5, color = "grey25", fontface = "italic") +
    annotate("text", x = 2.0, y = y_max * 0.985, label = "**  Depleted vs background",
             hjust = 1, size = 5, color = "grey25", fontface = "italic") +
    labs(
      title    = paste0(group_label, " cCpGs (n=", n_target,
                        ") vs Non-significant Background (n=", n_bg_ns, ")"),
      subtitle = paste0("Permutation test - ", N_PERM, " permutations"),
      x = NULL, y = "Odds Ratio (Enrichment)",
      caption = "* nominal empirical permutation p < 0.05. FDR values are in the accompanying table."
    ) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle      = element_text(size = 12, hjust = 0.5, color = "grey40"),
      plot.caption       = element_text(size = 10, color = "grey40", hjust = 0),
      axis.text.y        = element_text(size = 15, face = "bold", color = "black"),
      axis.text.x        = element_text(size = 14),
      axis.title.x       = element_text(size = 15, face = "bold", margin = margin(t = 8)),
      legend.title       = element_text(size = 14, face = "bold"),
      legend.text        = element_text(size = 13),
      legend.key.height  = unit(1.8, "lines"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.5),
      plot.margin        = margin(14, 30, 14, 14)
    )

  ggsave(out_file, p, width = 14, height = 11, dpi = 300)
  cat("Saved figure:", basename(out_file), "\n")
}

all_max <- ceiling(max(res_all_perm$nlog10p) * 10) / 10
make_perm_plot(res_all_perm, "All",
               file.path(OUT_DIR, "SI_all_vs_nonsig_permutation.png"),
               scale_max = all_max, scale_breaks = pretty(c(0, all_max), n = 4))

make_perm_plot(res_pos_perm, "Positive-Beta",
               file.path(OUT_DIR, "SI_pos_vs_nonsig_permutation.png"),
               scale_max = global_max_nlog, scale_breaks = global_breaks)

make_perm_plot(res_neg_perm, "Negative-Beta",
               file.path(OUT_DIR, "SI_neg_vs_nonsig_permutation.png"),
               scale_max = global_max_nlog, scale_breaks = global_breaks)


# --- positive against negative, on one panel ---------------------------------
# Ordered by the gap between the two odds ratios, so the categories where the
# two directions behave most differently sit at the top. Uses the permutation
# results already computed.

group_colors <- c("Positive-Beta cCpGs" = "#BB9457",
                  "Negative-Beta cCpGs" = "#344E41")

res_pos_div <- res_pos_perm %>% mutate(Group = "Positive-Beta cCpGs")
res_neg_div <- res_neg_perm %>% mutate(Group = "Negative-Beta cCpGs")

cat_order_div <- res_neg_div %>%
  left_join(res_pos_div %>% select(Category, OR_pos = OR), by = "Category") %>%
  mutate(abs_diff = abs(OR - OR_pos)) %>%
  arrange(abs_diff) %>%
  pull(Category)

res_both <- bind_rows(res_pos_div, res_neg_div) %>%
  mutate(
    Category = factor(Category, levels = cat_order_div),
    Group    = factor(Group, levels = c("Positive-Beta cCpGs", "Negative-Beta cCpGs"))
  )

prom_y <- which(cat_order_div == "promoters")

p_div <- ggplot(res_both, aes(x = Category, y = OR, fill = Group)) +
  annotate("rect", xmin = prom_y - 0.5, xmax = prom_y + 0.5,
           ymin = 0, ymax = Inf, fill = "#FFF3CD", alpha = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.72),
           width = 0.68, color = "black", linewidth = 0.35) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1.1) +
  geom_text(aes(label = sig, y = OR + 0.03),
            position = position_dodge(width = 0.72),
            hjust = -0.1, size = 6.5, fontface = "bold", color = "black") +
  coord_flip(ylim = c(0, max(res_both$OR) * 1.42)) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(breaks = seq(0, 2, 0.25)) +
  labs(
    title    = "Annotation Enrichment: Positive-Beta vs Negative-Beta cCpGs",
    subtitle = paste0("vs Non-significant Background (n=", n_bg_ns,
                      ")  |  Permutation test  |  ", N_PERM, " permutations"),
    x = NULL, y = "Odds Ratio (vs Background)", fill = NULL,
    caption = paste0(
      "Ordered by |ORneg - ORpos|, most dissociated at top. Yellow band: promoters.\n",
      "* nominal permutation p < 0.05 (enriched);  ** nominal p < 0.05 (depleted).\n",
      "Benjamini-Hochberg FDR values across the fifteen categories are reported in the accompanying tables."
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title         = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle      = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption       = element_text(size = 9.5, color = "grey35", hjust = 0, margin = margin(t = 8)),
    axis.text.y        = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x        = element_text(size = 13),
    axis.title.x       = element_text(size = 14, face = "bold", margin = margin(t = 8)),
    legend.position    = "bottom",
    legend.text        = element_text(size = 13, face = "bold"),
    legend.key.size    = unit(1.2, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.5),
    plot.margin        = margin(14, 30, 14, 14)
  )

ggsave(file.path(OUT_DIR, "SI_pos_neg_diverging_enrichment.png"), p_div,
       width = 14, height = 12, dpi = 300)
cat("Saved: SI_pos_neg_diverging_enrichment.png\n")

cat("\nWritten to ", OUT_DIR, "\n")
