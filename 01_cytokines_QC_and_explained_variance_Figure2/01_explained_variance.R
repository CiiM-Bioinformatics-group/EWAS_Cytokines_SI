# -----------------------------------------------------------------------------
# Variance in stimulation-induced cytokine responses explained by
# DNA methylation and by genetics
#
# For every cytokine-stimulus pair we fit three models on the principal
# components of each layer: methylation only, genetics only, and both together.
# The layers are also residualised on each other to separate the shared signal
# from the part that is specific to one layer.
#
# Input : cytokine matrix, methylation PCA, genotype PCA and covariates
# Output: adjusted R2 per cytokine-stimulus pair, the Figure 2 panels and the
#         Wilcoxon tests reported in the text
#
# Set work_dir below, then run top to bottom.
# -----------------------------------------------------------------------------

rm(list = ls())

work_dir <- "data/SI_cytokines"
setwd(work_dir)

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)
library(patchwork)
library(ggpubr)


# -----------------------------------------------------------------------------
# Cytokines: keep the conditions used in the analysis and drop the ones with
# too few measured samples
# -----------------------------------------------------------------------------

cyto2 <- readRDS("match_cyto2.rds")
colnames(cyto2) <- gsub("-", ".", colnames(cyto2))

# all.txt lists the cytokine_stimulus conditions to analyse. The unstimulated
# covid control is not a condition of interest.
valid_conditions <- readLines("all.txt")
valid_conditions <- gsub("-", ".", valid_conditions)
valid_conditions <- valid_conditions[!grepl("cov\\.ctrl", valid_conditions)]

cyto2_filtered <- cyto2[, colnames(cyto2) %in% valid_conditions]
cyto2 <- cyto2_filtered

sample_size_data <- data.frame(
  Column      = character(),
  Sample_Size = numeric(),
  stringsAsFactors = FALSE
)

for (colname in colnames(cyto2)) {
  sample_size <- sum(!is.na(cyto2[[colname]]))
  sample_size_data <- rbind(sample_size_data,
                            data.frame(Column = colname, Sample_Size = sample_size))
}

print(sample_size_data)

# A cytokine needs at least 200 measured samples to give a stable R2
valid_columns <- sample_size_data[sample_size_data$Sample_Size > 200, "Column"]
cyto_filtered <- cyto2[, colnames(cyto2) %in% valid_columns]
cyto2 <- cyto_filtered
dim(cyto_filtered)

png("samples_size_all_cyto.png")
ggplot(sample_size_data, aes(x = reorder(Column, Sample_Size), y = Sample_Size)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Distribution of Sample Sizes by Column",
       x = "Column",
       y = "Sample Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# -----------------------------------------------------------------------------
# Methylation principal components
# -----------------------------------------------------------------------------

PCA  <- readRDS("SI_PCAV2.rds")
cov2 <- readRDS("match_cov2.rds")

PCA2 <- as.data.frame(PCA$x)
PCA2 <- PCA2[1:140]          # first 140 methylation PCs

explained_variance            <- (PCA$sdev^2) / sum(PCA$sdev^2)
cumulative_explained_variance <- cumsum(explained_variance)

variance_data <- data.frame(
  Principal_Component           = 1:length(cumulative_explained_variance),
  Cumulative_Explained_Variance = cumulative_explained_variance
)
final_cumulative_variance_100 <- cumulative_explained_variance[140]

png("comulative_explained_150PCs.png")
ggplot(variance_data[1:140, ], aes(x = Principal_Component, y = Cumulative_Explained_Variance)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Cumulative Explained Variance by Principal Component",
       x = "Principal Component",
       y = "Cumulative Explained Variance") +
  theme_minimal()
dev.off()

# Put the methylation PCs in the same sample order as the cytokines
cyto2_rownames <- rownames(cyto2)
PCA2_subset    <- PCA2[rownames(PCA2) %in% cyto2_rownames, ]
PCA2_ordered   <- PCA2_subset[cyto2_rownames, ]
head(rownames(PCA2_ordered))

PCA2 <- PCA2_ordered

data_merged <- merge(PCA2, cyto2, by = "row.names", all = TRUE)
rownames(data_merged) <- data_merged$Row.names


# -----------------------------------------------------------------------------
# Genotype principal components
# -----------------------------------------------------------------------------

PCSG <- readRDS("PCA_SNPs_moreThan50.rds")
PCA  <- PCSG

PCA3 <- as.data.frame(PCA$x)
PCA3 <- PCA3[1:250]          # first 250 genotype PCs

explained_variance            <- (PCA$sdev^2) / sum(PCA$sdev^2)
cumulative_explained_variance <- cumsum(explained_variance)

variance_data <- data.frame(
  Principal_Component           = 1:length(cumulative_explained_variance),
  Cumulative_Explained_Variance = cumulative_explained_variance
)
final_cumulative_variance_100 <- cumulative_explained_variance[250]

rownames(data_merged) <- data_merged$Row.names


# -----------------------------------------------------------------------------
# Split the merged table back into PCs and cytokines, then keep the samples
# that are present in the genotype data as well
# -----------------------------------------------------------------------------

pcs_data      <- data_merged[, grep("PC", colnames(data_merged))]
cytokine_data <- data_merged[, !colnames(data_merged) %in% colnames(pcs_data)]
cytokine_data <- as.data.frame(lapply(cytokine_data, function(x) as.numeric(as.character(x))))
cytokine_data <- cytokine_data[, colnames(cytokine_data) != "Row.names"]

combined_data <- cbind(cytokine_data, pcs_data, cov2)

subset_PCA3 <- PCA3[rownames(combined_data), , drop = FALSE]

common_rows          <- intersect(rownames(subset_PCA3), rownames(combined_data))
subset_PCA3_common   <- subset_PCA3[common_rows, ]
combined_data_common <- combined_data[common_rows, ]

identical(rownames(subset_PCA3_common), rownames(combined_data_common))
dim(subset_PCA3_common)
dim(combined_data_common)

methylation_PCs <- combined_data_common[, grep("PC", colnames(combined_data_common))]
snps_PCs        <- subset_PCA3_common

identical(rownames(methylation_PCs), rownames(combined_data_common))
identical(rownames(snps_PCs), rownames(combined_data_common))


# All PCs of both layers are carried forward. Age and sex are not used to
# select PCs, they enter every cytokine model below as covariates instead.

dim(snps_PCs)
dim(methylation_PCs)


# -----------------------------------------------------------------------------
# Residualise each layer on the other. What is left of the genetic PCs after
# regressing out methylation is the genetic signal that methylation does not
# already carry, and the other way round.
# -----------------------------------------------------------------------------

residuals_genetics <- list()
for (i in 1:ncol(snps_PCs)) {
  model <- lm(snps_PCs[, i] ~ ., data = methylation_PCs)
  residuals_genetics[[i]] <- residuals(model)
}

residuals_genetics_df <- as.data.frame(do.call(cbind, residuals_genetics))
colnames(residuals_genetics_df) <- colnames(snps_PCs)

residuals_methylation <- list()
for (i in 1:ncol(methylation_PCs)) {
  model <- lm(methylation_PCs[, i] ~ ., data = snps_PCs)
  residuals_methylation[[i]] <- residuals(model)
}

residuals_methylation_df <- as.data.frame(do.call(cbind, residuals_methylation))
colnames(residuals_methylation_df) <- colnames(methylation_PCs)


# -----------------------------------------------------------------------------
# For each cytokine, keep the PCs that correlate with it and fit a linear model
# with age and sex as covariates. The adjusted R2 of that model is the variance
# explained by the layer.
#
# combined_data and cytokine_data are read from the environment.
# -----------------------------------------------------------------------------

cytokine_data <- combined_data %>%
  select(starts_with("pbmc"))

head(cytokine_data)

perform_analysis <- function(pcs_data, results_name) {

  matched_rows <- Reduce(intersect, list(row.names(combined_data),
                                         row.names(cytokine_data),
                                         row.names(pcs_data)))

  combined_data_subset <- combined_data[matched_rows, ]
  cytokine_data_subset <- cytokine_data[matched_rows, ]
  pcs_data_subset      <- pcs_data[matched_rows, ]

  combined_data_subset <- combined_data_subset[matched_rows, ]
  cytokine_data_subset <- cytokine_data_subset[matched_rows, ]
  pcs_data_subset      <- pcs_data_subset[matched_rows, ]

  print(paste("Missing values in cytokine data:", sum(is.na(cytokine_data_subset))))
  print(paste("Missing values in pcs data:", sum(is.na(pcs_data_subset))))
  print(paste("Missing values in age:", sum(is.na(combined_data_subset$age))))
  print(paste("Missing values in gender:", sum(is.na(combined_data_subset$gender))))

  results <- data.frame(cytokine = character(), R2 = numeric(), stringsAsFactors = FALSE)

  for (cytokine in colnames(cytokine_data_subset)) {

    significant_pcs <- c()

    for (pc in colnames(pcs_data_subset)) {
      spearman_cor <- cor.test(cytokine_data_subset[[cytokine]],
                               pcs_data_subset[[pc]],
                               method = "spearman")

      if (spearman_cor$p.value < 0.05) {
        significant_pcs <- c(significant_pcs, pc)
      }
    }

    if (length(significant_pcs) > 0) {
      formula_str <- paste(cytokine, "~", paste(significant_pcs, collapse = " + "),
                           "+ combined_data_subset$age + combined_data_subset$gender")

      model <- lm(as.formula(formula_str),
                  data = cbind(cytokine_data_subset, pcs_data_subset,
                               combined_data_subset[, c("age", "gender")]))

      r_squared <- summary(model)$adj.r.squared

      results <- rbind(results, data.frame(cytokine = cytokine, R2 = r_squared))
    }
  }

  sorted_results <- unique(results[order(-results$R2), ])

  assign(results_name, sorted_results, envir = .GlobalEnv)

  return(sorted_results)
}


# Residualised layers, then each layer on its own, then both together
gresults      <- perform_analysis(pcs_data = residuals_genetics_df,    results_name = "gresults")
mresults      <- perform_analysis(pcs_data = residuals_methylation_df, results_name = "mresults")
mresults_FULL <- perform_analysis(pcs_data = methylation_PCs,          results_name = "mresults_FULL")
gresults_FULL <- perform_analysis(pcs_data = snps_PCs,                 results_name = "gresults_FULL")

joint_PCs    <- cbind(methylation_PCs, snps_PCs)
jointresults <- perform_analysis(pcs_data = joint_PCs, results_name = "jointresults")


gresults_FULL <- gresults_FULL[order(gresults_FULL$cytokine), ]
gresults      <- gresults[order(gresults$cytokine), ]
mresults      <- mresults[order(mresults$cytokine), ]
mresults_FULL <- mresults_FULL[order(mresults_FULL$cytokine), ]
jointresults  <- jointresults[order(jointresults$cytokine), ]

# Collect the five models into one table
gresults$methylation_raw      <- mresults$R2
gresults$genetics_raw         <- gresults$R2
gresults$methylation_original <- mresults_FULL$R2
gresults$genetics_original    <- gresults_FULL$R2
gresults$both                 <- jointresults$R2

df <- gresults

df2 <- df %>%
  tidyr::separate(cytokine, into = c("tissue", "timepoint", "cytokine_name", "stimulation"), sep = "_")

print(df2)


# -----------------------------------------------------------------------------
# Figure 2c: distribution of explained variance per layer
# -----------------------------------------------------------------------------

source_colors <- c("Methylation" = "#1f77b4",
                   "Genetics"    = "#ff7f0e",
                   "Both"        = "#2ca02c")

df3 <- df2

df3_long <- df3 %>%
  pivot_longer(cols = c(methylation_original, genetics_original, both),
               names_to  = "Source",
               values_to = "VarianceExplained") %>%
  mutate(Source = recode(Source,
                         "methylation_original" = "Methylation",
                         "genetics_original"    = "Genetics",
                         "both"                 = "Both"))

df3_long$VarianceExplained <- df3_long$VarianceExplained * 100

# The unstimulated conditions are not part of the comparison
df3_filtered <- df3_long %>%
  filter(!stimulation %in% c("rpmi", "cov.ctrl"))

# Only cytokines that have a value for all three models
df3_complete <- df3_filtered %>%
  group_by(tissue, timepoint, cytokine_name, stimulation) %>%
  filter(!any(is.na(VarianceExplained))) %>%
  ungroup()

df3_long <- df3_complete

p_all <- ggplot(df3_long, aes(x = Source, y = VarianceExplained, fill = Source)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.9) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", color = "black", size = 0.3) +
  scale_fill_manual(values = source_colors) +
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 40), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance explained (%)") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(size = 16, color = "black", angle = 0),
    axis.text.y  = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16),
    axis.line    = element_line(color = "black"),
    plot.title   = element_blank(),
    plot.margin  = margin(5, 5, 5, 5),
    legend.position = "none"
  )

ggsave("explained_variance_panel_c.png", plot = p_all,
       width = 3.2, height = 4.5, dpi = 600)


# Same panel, one file per stimulation
unique_stimulations <- unique(df3_long$stimulation)

for (stim in unique_stimulations) {
  df_subset <- df3_long %>% filter(stimulation == stim)

  p <- ggplot(df_subset, aes(x = Source, y = VarianceExplained, fill = Source)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.9) +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", color = "black", size = 0.3) +
    scale_fill_manual(values = source_colors) +
    scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 40), expand = c(0, 0)) +
    labs(x = NULL, y = "Variance explained (%)") +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x  = element_text(size = 16, color = "black"),
      axis.text.y  = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 16),
      axis.line    = element_line(color = "black"),
      plot.title   = element_blank(),
      legend.position = "none"
    )

  ggsave(filename = paste0("explained_variance_", stim, ".png"),
         plot = p, width = 5.5, height = 5.5, dpi = 600)
}


# -----------------------------------------------------------------------------
# Methylation and genetics stacked per cytokine, one panel per stimulation
# -----------------------------------------------------------------------------

stack_colors <- c("methylation" = "#1f77b4",
                  "genetics"    = "#ff7f0e")

df3_long_stack <- df3 %>%
  filter(!stimulation %in% c("rpmi", "cov.ctrl", "cpg")) %>%
  filter(!is.na(methylation_original) & !is.na(genetics_original) & !is.na(both)) %>%
  mutate(cytokine = cytokine_name) %>%
  select(cytokine, stimulation,
         methylation = methylation_original,
         genetics    = genetics_original,
         both        = both)

# Order cytokines by their total explained variance under LPS and keep that
# order in every panel, so the panels can be read side by side
reference_stim <- "lps"
cytokine_order <- df3_long_stack %>%
  filter(stimulation == reference_stim) %>%
  mutate(total = methylation + genetics) %>%
  arrange(desc(total)) %>%
  pull(cytokine) %>%
  unique()

df3_stacked_long_explained <- df3_long_stack %>%
  pivot_longer(cols = c(methylation, genetics),
               names_to = "Source", values_to = "VarianceExplained") %>%
  mutate(cytokine = factor(cytokine, levels = cytokine_order))

# One shared y limit across panels
max_explained <- df3_stacked_long_explained %>%
  group_by(cytokine, stimulation) %>%
  summarise(total_explained = sum(VarianceExplained), .groups = "drop") %>%
  pull(total_explained) %>%
  max(na.rm = TRUE)

max_explained <- ceiling(max_explained * 10) / 10

stim_list <- unique(df3_stacked_long_explained$stimulation)

for (stim in stim_list) {
  df_plot <- df3_stacked_long_explained %>% filter(stimulation == stim)

  p <- ggplot(df_plot, aes(x = cytokine, y = VarianceExplained, fill = Source)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = stack_colors) +
    coord_flip() +
    scale_y_continuous(limits = c(0, max_explained)) +
    labs(title = paste("Explained Variance -", stim),
         x = NULL, y = "Variance explained (proportion)") +
    theme_minimal(base_size = 12) +
    theme(
      legend.title = element_blank(),
      axis.text.y  = element_text(size = 15),
      plot.title   = element_text(hjust = 0.5)
    )

  plot_height <- max(6, length(unique(df_plot$cytokine)) * 0.25)

  ggsave(paste0("explained_variance_global_scale_", stim, ".png"),
         plot   = p,
         width  = 6,
         height = plot_height,
         dpi    = 300)
}


# -----------------------------------------------------------------------------
# Genetic variance per cytokine, with one point per stimulation
# -----------------------------------------------------------------------------

stim_colors <- c(
  "LPS"      = "#e41a1c",
  "Pam3cys"  = "#377eb8",
  "Cov.N"    = "#4daf4a",
  "Cov.ctrl" = "#984ea3",
  "rpmi"     = "#bcbd22"
)

genetic_df <- df3_long_stack %>%
  filter(stimulation != "cpg") %>%
  mutate(
    stimulation = case_when(
      stimulation == "lps"      ~ "LPS",
      stimulation == "pam3cys"  ~ "Pam3cys",
      stimulation == "cov.ctrl" ~ "Cov.ctrl",
      stimulation == "cov.N"    ~ "Cov.N",
      stimulation == "rpmi"     ~ "rpmi",
      TRUE ~ stimulation
    ),
    cytokine = factor(cytokine, levels = cytokine_order)
  )

p <- ggplot(genetic_df, aes(x = cytokine, y = genetics)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA, width = 0.6) +
  geom_jitter(aes(color = stimulation), width = 0.15, size = 3, alpha = 0.8) +
  scale_color_manual(values = stim_colors) +
  scale_y_continuous(limits = c(0, 0.3)) +
  coord_flip() +
  labs(
    title = "Genetic Variance Explained per Cytokine Across Stimulations",
    y = "Variance explained (genetics)", x = NULL,
    color = "Stimulation"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

plot_height <- max(6, length(unique(genetic_df$cytokine)) * 0.25)

ggsave("G_variance_per_cytokine_clean_boxplot.png",
       plot   = p,
       width  = 6,
       height = plot_height,
       dpi    = 300)


saveRDS(df, "4modesl_SI_V2.rds")


# -----------------------------------------------------------------------------
# Annotate the results table: tissue, timepoint, cytokine and stimulation are
# encoded in the column name, and cytokines are grouped by cell origin
# -----------------------------------------------------------------------------

df <- readRDS("4modesl_SI_V2.rds")
head(df)

df_annotated <- df %>%
  separate(
    cytokine,
    into   = c("tissue", "time", "cytokine_name", "stimulation"),
    sep    = "_",
    remove = FALSE
  )

table(df_annotated$tissue)
table(df_annotated$time)
table(df_annotated$stimulation)

cytokine_groups <- list(
  chemokine  = c("ip10", "mig", "rantes", "sdf1a", "groa",
                 "mcp1", "mcp3", "mip1a", "mip1b", "eotaxin", "ctack"),
  monocyte   = c("il1a", "il1b", "il6", "tnfa", "tnfb", "il1ra",
                 "il18", "gcsf", "gmcsf", "mcsf", "mif",
                 "mip1a", "mip1b", "il8", "mcp1", "mcp3",
                 "lif", "scf", "eotaxin", "trail"),
  tcell      = c("ifng", "il17", "il4", "il5", "il13", "il2", "il7",
                 "il9", "il15", "il16", "il12.p40", "il12.p70", "il2ra"),
  regulatory = c("il10", "il1ra", "il13"),
  growth     = c("vegf", "pdgfbb", "fgf.basic", "hgf", "scgfb", "bngf")
)

# A cytokine can belong to more than one group
assign_group <- function(cyt) {
  hits <- names(Filter(function(x) cyt %in% x, cytokine_groups))
  if (length(hits) == 0) return(NA)
  paste(hits, collapse = ";")
}

df_annotated$cytokine_group <- sapply(df_annotated$cytokine_name, assign_group)

table(df_annotated$cytokine_group, useNA = "ifany")


# -----------------------------------------------------------------------------
# Does adding methylation to genetics buy anything beyond the better of the two
# single layers? Tested per stimulation, one sided.
# -----------------------------------------------------------------------------

df_simple <- df_annotated %>%
  mutate(
    best_single = pmax(genetics_original, methylation_original),
    delta_both  = both - best_single
  )

wilcox_by_stim <- df_simple %>%
  group_by(stimulation) %>%
  summarise(
    p_value       = wilcox.test(delta_both, mu = 0, alternative = "greater")$p.value,
    median_delta  = median(delta_both),
    prop_positive = mean(delta_both > 0),
    n             = n(),
    .groups = "drop"
  )

wilcox_by_stim

# Adjusted R2 penalises extra predictors, so the combined model can come out
# lower when the two layers carry the same information.


# -----------------------------------------------------------------------------
# Does residualising cost each layer any explained variance? Comparing the
# full model against the residualised one, per stimulation.
# -----------------------------------------------------------------------------

df_delta <- df %>%
  mutate(
    delta_meth = methylation_original - methylation_raw,
    delta_gen  = genetics_original   - genetics_raw
  )

wilcox_meth_by_stim <- df_delta %>%
  group_by(stimulation) %>%
  summarise(
    p_value       = wilcox.test(delta_meth, mu = 0, alternative = "greater")$p.value,
    median_delta  = median(delta_meth),
    prop_positive = mean(delta_meth > 0),
    n             = n(),
    .groups = "drop"
  )

wilcox_meth_by_stim

wilcox_gen_by_stim <- df_delta %>%
  group_by(stimulation) %>%
  summarise(
    p_value       = wilcox.test(delta_gen, mu = 0, alternative = "greater")$p.value,
    median_delta  = median(delta_gen),
    prop_positive = mean(delta_gen > 0),
    n             = n(),
    .groups = "drop"
  )

wilcox_gen_by_stim


# -----------------------------------------------------------------------------
# Numbers quoted in the manuscript text
# -----------------------------------------------------------------------------

gen_more_than_both_and_meth <- sum(df$genetics_original > df$methylation_original &
                                     df$genetics_original > df$both)
cat("Genetics > Methylation & Joint:", gen_more_than_both_and_meth, "\n")

meth_more_than_both_and_gen <- sum(df$methylation_original > df$genetics_original &
                                     df$methylation_original > df$both)
cat("Methylation > Genetics & Joint:", meth_more_than_both_and_gen, "\n")

both_more_than_meth_and_gen <- sum(df$both > df$genetics_original &
                                     df$both > df$methylation_original)
cat("Joint > Methylation & Genetics:", both_more_than_meth_and_gen, "\n")

meth_more_than_gen <- sum(df$methylation_original > df$genetics_original)
cat("Cytokines where methylation explains more variance than genetics:", meth_more_than_gen, "\n")

mean_methylation <- mean(df$methylation_original, na.rm = TRUE)
mean_genetics    <- mean(df$genetics_original, na.rm = TRUE)
mean_joint       <- mean(df$both, na.rm = TRUE)

max_methylation <- max(df$methylation_original, na.rm = TRUE)
max_genetics    <- max(df$genetics_original, na.rm = TRUE)
max_joint       <- max(df$both, na.rm = TRUE)

cat("Mean variance explained:\n")
cat("Methylation-only: ", round(mean_methylation * 100, 2), "%\n")
cat("Genetics-only:    ", round(mean_genetics * 100, 2), "%\n")
cat("Joint model:      ", round(mean_joint * 100, 2), "%\n\n")

cat("Maximum variance explained:\n")
cat("Methylation-only: ", round(max_methylation * 100, 2), "%\n")
cat("Genetics-only:    ", round(max_genetics * 100, 2), "%\n")
cat("Joint model:      ", round(max_joint * 100, 2), "%\n")


# Which cytokine-stimulus pair reaches the maximum in each model
get_max_info <- function(df, col_name) {
  max_val <- max(df[[col_name]], na.rm = TRUE)
  max_row <- df[df[[col_name]] == max_val, ]
  return(list(max = max_val, row = max_row))
}

max_meth_info     <- get_max_info(df, "methylation_original")
max_genetics_info <- get_max_info(df, "genetics_original")
max_joint_info    <- get_max_info(df, "both")

cat("Maximum variance explained:\n")
cat("Methylation-only:\n")
cat("  Value: ", round(max_meth_info$max * 100, 2), "%\n")
cat("  Cytokine-Stimulation:", max_meth_info$row$cytokine, "\n\n")

cat("Genetics-only:\n")
cat("  Value: ", round(max_genetics_info$max * 100, 2), "%\n")
cat("  Cytokine-Stimulation:", max_genetics_info$row$cytokine, "\n\n")

cat("Joint model:\n")
cat("  Value: ", round(max_joint_info$max * 100, 2), "%\n")
cat("  Cytokine-Stimulation:", max_joint_info$row$cytokine, "\n")

n_gen_over_meth <- sum(df$genetics_original > df$methylation_original, na.rm = TRUE)
cat("Cytokines where genetics explains more variance than methylation:",
    n_gen_over_meth,
    "(", round(100 * n_gen_over_meth / nrow(df), 1), "% )\n")


# -----------------------------------------------------------------------------
# Figure 2e: methylation and genetics side by side per cytokine,
# faceted by stimulation
# -----------------------------------------------------------------------------

df_dot <- df3 %>%
  select(cytokine_name, stimulation,
         methylation_original, genetics_original) %>%
  pivot_longer(cols = c(methylation_original, genetics_original),
               names_to = "source", values_to = "variance_explained") %>%
  mutate(source = recode(source,
                         "methylation_original" = "Methylation",
                         "genetics_original"    = "Genetics"))

p_dot <- ggplot(df_dot, aes(x = variance_explained,
                            y = reorder(cytokine_name, variance_explained))) +
  geom_point(aes(color = source), size = 3, position = position_dodge(width = 0.6)) +
  facet_wrap(~stimulation, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Methylation" = "#1f77b4", "Genetics" = "#ff7f0e")) +
  labs(x = "Explained Variance (%)", y = "Cytokine", color = "Source") +
  theme_bw() +
  theme(
    strip.text   = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18),
    axis.text.x  = element_text(size = 18)
  )

ggsave("panel_e_dotplot_horizontal.png",
       plot   = p_dot,
       width  = 30, height = 10, dpi = 500)


# -----------------------------------------------------------------------------
# Compact panels for the multi-panel figure layout
# -----------------------------------------------------------------------------

nature_colors <- c(
  "Genetics"    = "#E68613",
  "Methylation" = "#1F78B4"
)

df_long <- df_annotated %>%
  pivot_longer(
    cols      = c(genetics_original, methylation_original),
    names_to  = "source",
    values_to = "variance"
  ) %>%
  mutate(
    source = recode(source,
                    genetics_original    = "Genetics",
                    methylation_original = "Methylation"),
    source = factor(source, levels = c("Genetics", "Methylation"))
  )

# Two themes, the bar panels keep vertical guides, the boxplots have none
theme_nature_bars <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.line          = element_line(color = "black", linewidth = 0.3),
      axis.ticks         = element_line(color = "black", linewidth = 0.3),
      axis.title         = element_text(size = 11),
      axis.text          = element_text(size = 9),
      plot.title         = element_text(size = 12, face = "bold"),
      legend.title       = element_blank(),
      legend.position    = "top",
      legend.key.size    = unit(0.4, "cm"),
      plot.margin        = margin(5, 5, 5, 5)
    )
}

theme_nature_box <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.3),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      plot.title       = element_text(size = 13, face = "bold"),
      legend.title     = element_blank(),
      legend.position  = "top",
      plot.margin      = margin(5, 5, 5, 5)
    )
}

# Cytokines ordered by methylation R2, highest at the top
make_steak_plot <- function(stim) {

  df_stim <- df_long %>% filter(stimulation == stim)

  cytokine_order <- df_stim %>%
    filter(source == "Methylation") %>%
    arrange(desc(variance)) %>%
    pull(cytokine_name)

  df_stim <- df_stim %>%
    mutate(cytokine_name = factor(cytokine_name, levels = cytokine_order))

  ggplot(df_stim, aes(x = variance, y = cytokine_name, fill = source)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.6) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = nature_colors) +
    labs(
      title = paste("Explained Variance -", stim),
      x = "Variance explained (proportion)",
      y = ""
    ) +
    theme_nature_bars()
}

unique_stims <- unique(df_annotated$stimulation)

for (stim in unique_stims) {
  p <- make_steak_plot(stim)
  ggsave(
    filename = paste0("steakplot_", stim, "_compact_nature.png"),
    plot     = p,
    width    = 4.5,
    height   = 7,
    dpi      = 300
  )
}


# Genetics against methylation per stimulation, paired Wilcoxon on the panel
df_long2 <- df_annotated %>%
  select(cytokine_name, stimulation, genetics_original, methylation_original) %>%
  pivot_longer(
    cols      = c(genetics_original, methylation_original),
    names_to  = "source",
    values_to = "variance"
  ) %>%
  mutate(
    source = recode(source,
                    genetics_original    = "Genetics",
                    methylation_original = "Methylation"),
    source = factor(source, levels = c("Genetics", "Methylation"))
  )

make_boxplot <- function(stim) {

  df_stim <- df_long2 %>% filter(stimulation == stim)

  pval <- wilcox.test(
    df_stim$variance[df_stim$source == "Genetics"],
    df_stim$variance[df_stim$source == "Methylation"],
    paired = TRUE
  )$p.value

  p_label <- paste0("Wilcoxon p = ", signif(pval, 3))

  ggplot(df_stim, aes(x = source, y = variance, fill = source)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.1, size = 1.7, alpha = 0.7) +
    scale_fill_manual(values = nature_colors) +
    labs(
      title = paste("Genetics vs. Methylation -", stim),
      y = "Explained variance",
      x = ""
    ) +
    annotate("text",
             x = 1.5,
             y = max(df_stim$variance) * 1.05,
             label = p_label,
             size = 3.3) +
    theme_nature_box() +
    ylim(0, max(df_stim$variance) * 1.15)
}

for (stim in unique_stims) {
  p <- make_boxplot(stim)
  ggsave(
    filename = paste0("boxplot_", stim, "_nature.png"),
    plot     = p,
    width    = 4,
    height   = 4.5,
    dpi      = 300
  )
}


# -----------------------------------------------------------------------------
# Delta R2 of the combined model over the better single layer, per stimulation
# -----------------------------------------------------------------------------

df_delta <- df_annotated %>%
  mutate(
    best_single = pmax(genetics_original, methylation_original),
    delta_R2    = both - best_single
  ) %>%
  select(cytokine_name, stimulation, delta_R2)

wilcox_results <- df_delta %>%
  group_by(stimulation) %>%
  summarise(
    p = wilcox.test(delta_R2, mu = 0, alternative = "greater")$p.value
  ) %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

print(wilcox_results)

theme_nature_delta <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.3),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 10),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.margin      = margin(5, 5, 5, 5)
    )
}

df_plot <- df_delta %>%
  left_join(wilcox_results, by = "stimulation")

df_plot$stimulation <- factor(df_plot$stimulation)

p <- ggplot(df_plot, aes(x = stimulation, y = delta_R2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7, color = "grey30") +
  geom_text(
    aes(label = sig),
    y = max(df_plot$delta_R2) * 1.05,
    size = 5
  ) +
  labs(
    title = "Methylation adds stimulation-specific explained variance",
    x = "Stimulation",
    y = expression(Delta~R^2~"(both - best single layer)")
  ) +
  ylim(
    min(df_plot$delta_R2) * 1.2,
    max(df_plot$delta_R2) * 1.2
  ) +
  theme_nature_delta()

print(p)

ggsave(
  filename = "delta_R2_by_stimulation_nature.png",
  plot     = p,
  width    = 5.5,
  height   = 5.5,
  dpi      = 300
)
