# -----------------------------------------------------------------------------
# Quality control and normalisation of stimulation-induced (SI) cytokine levels
#
# Input : raw cytokine matrix (cytokine_stimulus x sample) and cohort phenotypes
# Output: filtered and rank-normalised cytokine matrices, QC plots,
#         24h / 7d subsets, correlation heatmaps and the age-sex distribution panel
#
# Set the three paths below before running. Everything else runs top to bottom.
# -----------------------------------------------------------------------------

data_dir <- "data/SI_cytokines"      # raw cytokine matrix + phenotype exports
meth_dir <- "data/methylation"       # methylation objects used for the age split
out_dir  <- file.path(data_dir, "aged_young")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(data_dir)

library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(Hmisc)
library(reshape)
library(pheatmap)
library(RColorBrewer)


# Raw matrix. Cytokines and stimulations are already recoded as cytokine_stimulus.
cytokines <- read.csv("cytokines_for_cqtl_researchset.tsv", sep = "\t")
dim(cytokines)

row.names(cytokines) <- cytokines$cytstim
cytokines <- as.data.frame(cytokines)
dat2 <- subset(cytokines, select = -cytstim)

# Samples in rows, cytokine_stimulus in columns
df <- t(dat2)


# -----------------------------------------------------------------------------
# Overview of the raw data: NAs, blanks and zero records per cytokine
# -----------------------------------------------------------------------------

zero_counts <- colSums(df == 0, na.rm = TRUE)
na_counts   <- colSums(is.na(df))
gap_counts  <- colSums(df == "", na.rm = TRUE)

result <- data.frame(
  Column     = colnames(df),
  NA_Count   = na_counts,
  Gap_Count  = gap_counts,
  Zero_Count = zero_counts
)
head(result)

saveRDS(result, "counts_na_zeros_gap_SI.rds")

# Zeros and blanks are below the detection limit of the assay, so treat them as missing
df_modified <- df
df_modified[df_modified == 0]  <- NA
df_modified[df_modified == ""] <- NA

na_counts <- colSums(is.na(df_modified))
result <- data.frame(
  Column   = colnames(df),
  NA_Count = na_counts
)
head(result)

df <- df_modified


# -----------------------------------------------------------------------------
# STEP 1: density plots per cytokine, 20 per page
# -----------------------------------------------------------------------------

create_density_plots <- function(df, cols_per_plot = 20) {
  num_cols <- ncol(df)

  for (i in seq(1, num_cols, by = cols_per_plot)) {
    col_indices <- i:min(i + cols_per_plot - 1, num_cols)

    plot_list <- lapply(col_indices, function(col_idx) {
      col_name <- colnames(df)[col_idx]
      ggplot(df, aes(x = .data[[col_name]])) +
        geom_density(na.rm = TRUE, fill = "blue", color = "blue") +
        labs(title = paste("cytokine: ", col_name), x = col_name) +
        theme_classic()
    })

    combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

    filename <- paste0("density_plots_", i, "_to_", min(i + cols_per_plot - 1, num_cols), ".png")
    ggsave(filename, plot = combined_plot, width = 14, height = 10, dpi = 300)
  }
}

create_density_plots(df)


# -----------------------------------------------------------------------------
# STEP 2: box plots per cytokine and Shapiro-Wilk test on the raw values
# -----------------------------------------------------------------------------

create_box_plots <- function(df, cols_per_plot = 20) {
  num_cols <- ncol(df)

  for (i in seq(1, num_cols, by = cols_per_plot)) {
    col_indices <- i:min(i + cols_per_plot - 1, num_cols)

    plot_list <- lapply(col_indices, function(col_idx) {
      col_name <- colnames(df)[col_idx]
      ggplot(df, aes(y = .data[[col_name]])) +
        geom_boxplot(fill = "blue", color = "blue") +
        labs(title = paste("Box Plot of", col_name), y = col_name) +
        theme_classic()
    })

    combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

    filename <- paste0("box_plots_", i, "_to_", min(i + cols_per_plot - 1, num_cols), ".png")
    ggsave(filename, plot = combined_plot, width = 14, height = 10, dpi = 300)
  }
}

create_box_plots(df)


data <- as.data.frame(df)

shapiro_results <- list()
for (col in colnames(data)) {
  test_result <- shapiro.test(data[[col]])
  shapiro_results[[col]] <- test_result
}

results_df <- data.frame(
  Variable = names(shapiro_results),
  W        = sapply(shapiro_results, function(x) x$statistic),
  p_value  = sapply(shapiro_results, function(x) x$p.value)
)
print(dim(results_df))

significant_columns <- which(as.numeric(results_df$p_value) < 0.05)
not_significnat     <- which(as.numeric(results_df$p_value) > 0.05)

write.csv(results_df, "shapiro_test_SI_raw.csv")

significant <- results_df[significant_columns, ]
not_sig     <- results_df[not_significnat, ]

print("signifincat SHAPIRO")
print(dim(significant))
print("NOT SIGNIFNCAT SHAPIRO")
print(dim(not_sig))


# -----------------------------------------------------------------------------
# STEP 3: drop cytokines missing in more than 50% of the samples
# -----------------------------------------------------------------------------

missing_percentage <- function(column) {
  mean(is.na(column)) * 100
}

missing_values_df <- data.frame(
  Variable           = colnames(data),
  Missing_Percentage = sapply(data, missing_percentage)
)

columns_above_50 <- missing_values_df %>%
  filter(Missing_Percentage > 50)

if (nrow(columns_above_50) > 0) {
  print("Columns with more than 50% missing values:")
  print(columns_above_50)
} else {
  print("No columns have more than 50% missing values.")
}

print("Columns >50% missing values:")
dim(columns_above_50)

columns_less_50 <- missing_values_df %>%
  filter(as.numeric(Missing_Percentage) < 50)

if (nrow(columns_less_50) > 0) {
  print("Columns with less than 50% missing values:")
  print(columns_less_50)
} else {
  print("No columns have less than 50% missing values.")
}

print("Columns <50% missing values:")
dim(columns_less_50)

print("Columns >50% missing values:")
dim(columns_above_50)

write.csv(columns_above_50, "Table1.Cytokines_missingMoreThan50.csv")
write.csv(columns_less_50,  "Table1.Cytokines_missinglessThan50.csv")

cleaned_data1 <- data %>%
  select(-all_of(columns_above_50$Variable))

print("Cleaned datasets from samples with NA more than 50%")
dim(cleaned_data1)

write.csv(cleaned_data1,    "Table3.cleand-missing.csv")
write.csv(t(cleaned_data1), "Table3.csv")


# -----------------------------------------------------------------------------
# STEP 4: drop cytokines where a single value is repeated in more than 50%
#         of the samples, these carry almost no information
# -----------------------------------------------------------------------------

redundancy_results <- data.frame(
  Column         = character(),
  Max_Redundancy = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in colnames(cleaned_data1)) {
  value_counts <- table(cleaned_data1[[col_name]][!is.na(cleaned_data1[[col_name]])])

  max_redundancy <- if (length(value_counts) > 0) max(value_counts) else 0

  redundancy_results <- rbind(redundancy_results,
                              data.frame(Column = col_name,
                                         Max_Redundancy = max_redundancy))
}

sample_size <- nrow(cleaned_data1)
redundancy_results$percentage <- (redundancy_results$Max_Redundancy / sample_size) * 100

removed <- redundancy_results[which(redundancy_results$percentage > 50), ]

if (nrow(removed) > 0) {
  cleaned_data2 <- cleaned_data1 %>%
    select(-all_of(removed$Column))
} else {
  cleaned_data2 <- cleaned_data1
}

print("datasets for Normalizations")
dim(cleaned_data2)

# Sample size and number of missing records per remaining cytokine
complete_counts <- colSums(!is.na(cleaned_data2))
na_counts       <- colSums(is.na(cleaned_data2))

complete_cases_df <- data.frame(
  Cytokine            = colnames(cleaned_data2),
  Non_Missing_Samples = complete_counts,
  Missing_Samples     = na_counts
)

print("Cytokines with number of samples (non-missing & missing):")
write.csv(complete_cases_df, "Complete_SI_cytokines_records_NA_stat.csv", row.names = TRUE)


# -----------------------------------------------------------------------------
# STEP 5: rank based inverse normal transformation
# -----------------------------------------------------------------------------

set.seed(123)

rank_normalize <- function(value) {
  qnorm((rank(value, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(value)))
}

normalized_data <- as.data.frame(lapply(cleaned_data2, rank_normalize))
row.names(normalized_data) <- row.names(cleaned_data2)

print("Normalized DataFrame:")
print(normalized_data)

write.csv(normalized_data, "normalized_data.csv", row.names = TRUE)


# -----------------------------------------------------------------------------
# STEP 6: same distribution checks on the normalised values
# -----------------------------------------------------------------------------

create_density_plots <- function(normalized_data, cols_per_plot = 20) {
  num_cols <- ncol(normalized_data)

  for (i in seq(1, num_cols, by = cols_per_plot)) {
    col_indices <- i:min(i + cols_per_plot - 1, num_cols)

    plot_list <- lapply(col_indices, function(col_idx) {
      col_name <- colnames(normalized_data)[col_idx]
      ggplot(normalized_data, aes(x = .data[[col_name]])) +
        geom_density(na.rm = TRUE, fill = "blue", color = "blue") +
        labs(title = col_name, x = col_name) +
        theme_classic()
    })

    combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

    filename <- paste0("density_plots_Normalized", i, "_to_", min(i + cols_per_plot - 1, num_cols), ".png")
    ggsave(filename, plot = combined_plot, width = 14, height = 10, dpi = 300)
  }
}

create_density_plots(normalized_data)


create_box_plots <- function(normalized_data, cols_per_plot = 20) {
  num_cols <- ncol(normalized_data)

  for (i in seq(1, num_cols, by = cols_per_plot)) {
    col_indices <- i:min(i + cols_per_plot - 1, num_cols)

    plot_list <- lapply(col_indices, function(col_idx) {
      col_name <- colnames(normalized_data)[col_idx]
      ggplot(normalized_data, aes(y = .data[[col_name]])) +
        geom_boxplot(fill = "blue", color = "blue") +
        labs(title = paste(".", col_name), y = col_name) +
        theme_classic()
    })

    combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

    filename <- paste0("Normalized_box_plots_", i, "_to_", min(i + cols_per_plot - 1, num_cols), ".png")
    ggsave(filename, plot = combined_plot, width = 14, height = 10, dpi = 300)
  }
}

create_box_plots(normalized_data)


rm(shapiro_results)
rm(significant)
rm(results_df)
rm(not_sig)

normalized_data <- as.data.frame(normalized_data)

shapiro_results <- list()
for (col in colnames(normalized_data)) {
  test_result <- shapiro.test(normalized_data[[col]])
  shapiro_results[[col]] <- test_result
}

results_df <- data.frame(
  Variable = names(shapiro_results),
  W        = sapply(shapiro_results, function(x) x$statistic),
  p_value  = sapply(shapiro_results, function(x) x$p.value)
)
print(dim(results_df))

significant_columns <- which(as.numeric(results_df$p_value) < 0.05)
not_significnat     <- which(as.numeric(results_df$p_value) > 0.05)

write.csv(results_df, "shapiro_test_SI_normalized.csv")

significant <- results_df[significant_columns, ]
not_sig     <- results_df[not_significnat, ]

print("signifincat SHAPIRO")
print(dim(significant))
print("NOT SIGNIFNCAT SHAPIRO")
print(dim(not_sig))


# make.names turned the dashes in the stimulation labels into dots, put them back
col1 <- colnames(normalized_data)
colnames(normalized_data) <- gsub("\\.(.+)", "-\\1", colnames(normalized_data))

saveRDS(normalized_data, "Normalized.Cytokines.rds")


# Raw distributions again, this time only for the cytokines that survived the filters
df <- as.data.frame(t(dat2))

columns_to_keep <- intersect(names(normalized_data), names(df))
df_subset <- df %>% select(all_of(columns_to_keep))

create_density_plots <- function(df_subset, cols_per_plot = 20) {
  num_cols <- ncol(df_subset)

  for (i in seq(1, num_cols, by = cols_per_plot)) {
    col_indices <- i:min(i + cols_per_plot - 1, num_cols)

    plot_list <- lapply(col_indices, function(col_idx) {
      col_name <- colnames(df_subset)[col_idx]
      ggplot(df_subset, aes(x = .data[[col_name]])) +
        geom_density(na.rm = TRUE, fill = "blue", color = "blue") +
        labs(title = paste("r.", col_name), x = col_name) +
        theme_classic()
    })

    combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

    filename <- paste0("density_plots_raw", i, "_to_", min(i + cols_per_plot - 1, num_cols), ".png")
    ggsave(filename, plot = combined_plot, width = 14, height = 10, dpi = 300)
  }
}

create_density_plots(df_subset)


# One histogram per normalised cytokine
c <- readRDS("Normalized.Cytokines.rds")
normalized_data <- c

for (col in names(c)) {
  if (is.numeric(c[[col]])) {
    filename <- paste0("histogram_", col, ".png")
    png(filename)
    hist(c[[col]],
         main   = paste("Histogram of", col),
         xlab   = col,
         col    = "lightblue",
         border = "black")
    dev.off()
  }
}


# -----------------------------------------------------------------------------
# Split by stimulation duration and report the missingness that is left
# -----------------------------------------------------------------------------

df <- normalized_data
dim(df)

data7 <- normalized_data[grepl("7d_", colnames(normalized_data))]
saveRDS(data7, "data7.rds")

data24 <- normalized_data[grepl("24", colnames(normalized_data))]
saveRDS(data24, "data24.rds")

write.csv(data7,     "data7.csv")
write.csv(data24,    "data24.csv")
write.csv(t(data7),  "tdata7.csv")

# 7 days
na_counts <- sapply(data7, function(x) sum(is.na(x)))

na_summary <- data.frame(
  Column   = names(na_counts),
  NA_Count = na_counts,
  stringsAsFactors = FALSE
)
print(na_summary)
dim(na_summary)

sample_size <- nrow(data7)
na_summary$NA_Count   <- as.numeric(na_summary$NA_Count)
na_summary$percentage <- (na_summary$NA_Count / sample_size) * 100

write.table(na_summary, "Precentage_of_missing_valuse_for7d_NormalizedCol.csv")

png("percentage_of_missing_values_for7d_normalised_cytokines.png",
    width = 10, height = 7, units = "in", res = 300)

ggplot(na_summary, aes(x = reorder(Column, percentage), y = percentage)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5, size = 1.5) +
  labs(
    title = "Percentage of Missing Values for Normalised Cytokines",
    x = "Column",
    y = "Percentage of Missing Values"
  ) +
  theme(axis.text.x = element_text(angle = 80, hjust = 1))

dev.off()

rm(na_summary)


# 24 hours
data24 <- normalized_data[grepl("24h_", colnames(normalized_data))]

na_counts <- sapply(data24, function(x) sum(is.na(x)))

na_summary <- data.frame(
  Column   = names(na_counts),
  NA_Count = na_counts,
  stringsAsFactors = FALSE
)
print(na_summary)
dim(na_summary)

sample_size <- nrow(data24)
na_summary$NA_Count   <- as.numeric(na_summary$NA_Count)
na_summary$percentage <- (na_summary$NA_Count / sample_size) * 100

write.table(na_summary, "Precentage_of_missing_valuse_for24d_NormalizedCol.csv")

png("percentage_of_missing_values_for24h_normalised_cytokines.png",
    width = 10, height = 7, units = "in", res = 300)

ggplot(na_summary, aes(x = reorder(Column, percentage), y = percentage)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5, size = 1.5) +
  labs(
    title = "Percentage of Missing Values for Normalised Cytokines",
    x = "Column",
    y = "Percentage of Missing Values"
  ) +
  theme(axis.text.x = element_text(angle = 80, hjust = 1))

dev.off()

rm(na_summary)


# ---------------------- normalisation and plots are done ---------------------


# -----------------------------------------------------------------------------
# Phenotypes: age, sex, smoking and BMI
# -----------------------------------------------------------------------------

pheno <- read.csv(file.path("phenotypes", "cohort_phenotype_export.csv"), sep = ";")
pheno <- as.data.frame(pheno)

pheno1 <- data.frame(as.numeric(pheno$modys_age), pheno$modys_gender,
                     pheno$blood_smoke, pheno$anthro_bmi)
colnames(pheno1) <- c("modys_age", "modys_gender", "blood_smoke", "anthro_bmi")

png("age.dist.png", width = 10, height = 7, units = "in", res = 300)
ggplot(pheno1, aes(x = modys_age)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Count")
dev.off()

png("Gender.dist.png", width = 10, height = 7, units = "in", res = 300)
ggplot(pheno1, aes(x = modys_gender)) +
  geom_bar(fill = "lightgreen", color = "black") +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")
dev.off()

png("Smoke.dist.png", width = 10, height = 7, units = "in", res = 300)
ggplot(pheno1, aes(x = blood_smoke)) +
  geom_bar(fill = "lightgreen", color = "black") +
  labs(title = "Smoke", x = "Smoking status", y = "Count")
dev.off()

png("bmi.dist.png", width = 10, height = 7, units = "in", res = 300)
ggplot(pheno1, aes(x = as.numeric(anthro_bmi))) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
  labs(title = "anthro_bmi", x = "anthro_bmi", y = "Count")
dev.off()


# -----------------------------------------------------------------------------
# Spearman correlation between cytokines, 24h and 7d separately
# -----------------------------------------------------------------------------

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row    = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor    = (cormat)[ut],
    p      = pmat[ut]
  )
}

# The labels look like tissue_time_cytokine_stimulation, so field 3 is the
# cytokine and everything after it is the stimulus
extract_info <- function(label) {
  parts <- strsplit(label, "_")[[1]]
  cytokine    <- parts[3]
  stimulation <- paste(parts[4:length(parts)], collapse = "_")
  return(c(cytokine = cytokine, stimulation = stimulation))
}


# 24 hours
cyto_norm <- readRDS("Normalized.Cytokines.rds")
data24    <- normalized_data[grepl("24h_", colnames(cyto_norm))]
cyto_norm <- data24
head(cyto_norm)

res2 <- rcorr(as.matrix(cyto_norm), type = "spearman")

heatmap_data <- melt(res2$r, varnames = c("Row", "Column"), value.name = "Correlation")
head(heatmap_data)

row_info <- t(apply(data.frame(heatmap_data$Row), 1, extract_info))
col_info <- t(apply(data.frame(heatmap_data$Column), 1, extract_info))

heatmap_data <- heatmap_data %>%
  mutate(
    Row_cytokine       = row_info[, "cytokine"],
    Row_stimulation    = row_info[, "stimulation"],
    Column_cytokine    = col_info[, "cytokine"],
    Column_stimulation = col_info[, "stimulation"]
  )

row_annotation <- unique(data.frame(heatmap_data$Row,
                                    heatmap_data$Row_cytokine,
                                    heatmap_data$Row_stimulation))

row_annotation2 <- data.frame(row_annotation$heatmap_data.Row_cytokine,
                              row_annotation$heatmap_data.Row_stimulation)
row_annotation2$row_annotation.heatmap_data.Row_stimulation <-
  gsub("\\.", "-", row_annotation2$row_annotation.heatmap_data.Row_stimulation)

row.names(row_annotation2) <- row_annotation$heatmap_data.Row
colnames(row_annotation2)  <- c("Cytokines", "Stimulation")

data_matrix <- as.matrix(res2$r)

png("heatmap24h.png", width = 10, height = 7, units = "in", res = 300)
pheatmap(
  data_matrix,
  annotation_row = row_annotation2,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = FALSE,
  show_colnames  = TRUE,
  main = "Normalised_SI_cytokines_Time_24_Tissue_PBMC",
  height = 100, width = 100, fontsize = 6
)
dev.off()

rm(cyto_norm)
rm(data_matrix)
rm(row_annotation2)


# 7 days
cyto_norm <- readRDS("Normalized.Cytokines.rds")
data7     <- normalized_data[grepl("7d_", colnames(cyto_norm))]
cyto_norm <- data7
head(cyto_norm)

res2 <- rcorr(as.matrix(cyto_norm), type = "spearman")

heatmap_data <- melt(res2$r, varnames = c("Row", "Column"), value.name = "Correlation")
head(heatmap_data)

row_info <- t(apply(data.frame(heatmap_data$Row), 1, extract_info))
col_info <- t(apply(data.frame(heatmap_data$Column), 1, extract_info))

heatmap_data <- heatmap_data %>%
  mutate(
    Row_cytokine       = row_info[, "cytokine"],
    Row_stimulation    = row_info[, "stimulation"],
    Column_cytokine    = col_info[, "cytokine"],
    Column_stimulation = col_info[, "stimulation"]
  )

row_annotation <- unique(data.frame(heatmap_data$Row,
                                    heatmap_data$Row_cytokine,
                                    heatmap_data$Row_stimulation))

row_annotation2 <- data.frame(row_annotation$heatmap_data.Row_cytokine,
                              row_annotation$heatmap_data.Row_stimulation)
row_annotation2$row_annotation.heatmap_data.Row_stimulation <-
  gsub("\\.", "-", row_annotation2$row_annotation.heatmap_data.Row_stimulation)

row.names(row_annotation2) <- row_annotation$heatmap_data.Row
colnames(row_annotation2)  <- c("Cytokines", "Stimulation")

data_matrix <- as.matrix(res2$r)

png("heatmap7d.png", width = 10, height = 7, units = "in", res = 300)
pheatmap(
  data_matrix,
  annotation_row = row_annotation2,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = FALSE,
  show_colnames  = TRUE,
  main = "Normalised_SI_cytokines_Time_7d_Tissue_PBMC",
  height = 100, width = 100, fontsize = 6
)
dev.off()


# -----------------------------------------------------------------------------
# Match the phenotype identifiers to the research identifiers used in the
# cytokine matrix, then keep only the samples present in both
# -----------------------------------------------------------------------------

match <- read.csv(file.path("phenotypes", "id_mapping.tsv"), sep = "\t")

colnames(pheno1)  <- c("Age", "Gender", "Smoke", "BMI")
row.names(pheno1) <- pheno$user_id

matches  <- pheno1[rownames(pheno1) %in% as.character(match$RESIST.ID), ]
matches2 <- match[as.character(match$RESIST.ID) %in% rownames(pheno1), ]
row.names(matches2) <- matches2$RESIST.ID

matches2 <- matches2[order(rownames(matches2)), ]
matches  <- matches[order(rownames(matches)), ]
data24   <- data24[order(rownames(data24)), ]

matches$ID <- matches2$ResearchID
row.names(matches) <- matches$ID

Data24 <- data24[rownames(data24) %in% rownames(matches), ]
Pheno2 <- matches[rownames(matches) %in% rownames(Data24), ]

Pheno2 <- Pheno2[order(rownames(Pheno2)), ]
Data24 <- Data24[order(rownames(Data24)), ]

identical(rownames(Pheno2), rownames(Data24))

pheno3 <- data.frame(Pheno2$Age, Pheno2$Gender)
rownames(pheno3) <- rownames(Pheno2)
colnames(pheno3) <- c("Age", "Gender")

dim(Data24)
dim(pheno3)


# -----------------------------------------------------------------------------
# Aged / young stratification at 50 years, used by the downstream analyses
# -----------------------------------------------------------------------------

setwd(out_dir)

aged_df  <- pheno3 %>% filter(Age >= 50)
young_df <- pheno3 %>% filter(Age < 50)

dim(aged_df)
dim(young_df)

aged_samples  <- rownames(aged_df)
young_samples <- rownames(young_df)

aged_normalized  <- normalized_data[rownames(normalized_data) %in% aged_samples, ]
young_normalized <- normalized_data[rownames(normalized_data) %in% young_samples, ]

dim(aged_normalized)
dim(young_normalized)

saveRDS(aged_normalized,  "aged_cyto_normalized.rds")
saveRDS(young_normalized, "young_cyto_normalized.rds")
saveRDS(aged_df,          "aged_pheno_normalized.rds")
saveRDS(young_df,         "young_pheno_normalized.rds")

# Same split on the methylation matrices
setwd(meth_dir)

MYoung <- readRDS("match_mythelation3.rds")
MAged  <- readRDS("match_mythelation.rds")

aged_Methylation  <- MAged[rownames(MAged) %in% aged_samples, ]
young_Methylation <- MYoung[rownames(MYoung) %in% young_samples, ]

setwd(out_dir)

saveRDS(aged_Methylation,  "aged_methylation_SI.rds")
saveRDS(young_Methylation, "young_methylation_SI.rds")

# Ses (cell composition estimates) and cov2_df (array plate covariates) are built
# in the covariate script. Load them first, then run the block below to get the
# matching aged and young subsets.
#
# aged_Ses   <- Ses[rownames(Ses) %in% aged_samples, ]
# young_Ses  <- Ses[rownames(Ses) %in% young_samples, ]
# saveRDS(aged_Ses,  "aged_ses_SI.rds")
# saveRDS(young_Ses, "young_ses_SI.rds")
#
# row.names(cov2_df) <- cov2_df$sample_id
# aged_plate  <- cov2_df[rownames(cov2_df) %in% aged_samples, ]
# young_plate <- cov2_df[rownames(cov2_df) %in% young_samples, ]
# saveRDS(aged_plate,  "aged_plate_SI.rds")
# saveRDS(young_plate, "young_plate_SI.rds")


# -----------------------------------------------------------------------------
# Age and sex distribution of the cohort
# -----------------------------------------------------------------------------

pheno3$Age <- as.factor(pheno3$Age)
pheno3     <- pheno3 %>% filter(as.numeric(as.character(Age)) > 50)

# The export stores sex in German
pheno3$Gender <- recode(pheno3$Gender, "Männlich" = "Male", "Weiblich" = "Female")

age_gender_count <- pheno3 %>%
  group_by(Age, Gender) %>%
  summarise(Count = n()) %>%
  ungroup()

colors <- c("Male" = "#1f78b4", "Female" = "#e31a1c")

p <- ggplot(age_gender_count, aes(x = Age, y = Count, fill = Gender)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Age and Gender Distribution",
       x = "Age",
       y = "Count",
       fill = "Gender") +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    legend.text     = element_text(size = 14),
    legend.title    = element_text(size = 16),
    plot.title      = element_text(size = 18, hjust = 0.5),
    legend.position = "top",
    plot.margin     = margin(10, 10, 10, 10)
  )

# A4 landscape
ggsave("age_gender_distribution_A4_highres.png", plot = p,
       width = 11.7, height = 8.3, dpi = 300)
