# -----------------------------------------------------------------------------
# Epigenome-wide association study of stimulation-induced cytokine responses
#
# One cytokine per run. For every CpG a robust linear model is fitted with the
# cytokine as outcome and methylation as the predictor of interest:
#
#   cytokine ~ methylation + age + sex + cell composition + plate + SV1..SVn
#
# Usage:
#   Rscript EWAS_run.R <cytokine_name>
#
# The cytokine name must match a column of the cytokine matrix. Hyphens are
# converted to dots so that the name is a valid term in the model formula.
#
# Expected inputs in data_dir, all indexed by sample identifier:
#   M_values_matrix.rds     methylation M-values, samples in rows, CpGs in columns
#   cytokine_data.rds       normalised cytokine responses, samples in rows
#   covariates.rds          age, gender, bmi, smoke, SAMPLE_PLATE
#   methylation_pca.rds     prcomp object of the methylation matrix
#   cell_counts.rds         estimated cell proportions
#   seasonality.tsv         season_sin, season_cos, numDaysFromJan2013, ResearchID
#   surrogate_variables.rds surrogate variables estimated on the methylation data
# -----------------------------------------------------------------------------

rm(list = ls())

library(MASS)
library(lmtest)
library(sandwich)
library(foreach)
library(doParallel)


# --- settings ----------------------------------------------------------------

data_dir <- "data"
out_dir  <- "results"

num_cores  <- 10    # cores for the parallel loop over CpGs
batch_size <- 500   # CpGs per batch

methylation_file <- file.path(data_dir, "M_values_matrix.rds")
cytokine_file    <- file.path(data_dir, "cytokine_data.rds")
covariate_file   <- file.path(data_dir, "covariates.rds")
pca_file         <- file.path(data_dir, "methylation_pca.rds")
cellcount_file   <- file.path(data_dir, "cell_counts.rds")
seasonality_file <- file.path(data_dir, "seasonality.tsv")
sva_file         <- file.path(data_dir, "surrogate_variables.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- inputs ------------------------------------------------------------------

cat("Loading datasets...\n")

Mval  <- t(readRDS(methylation_file))
cyto2 <- readRDS(cytokine_file)
cov2  <- readRDS(covariate_file)
pcs   <- readRDS(pca_file)
db    <- readRDS(cellcount_file)
ses   <- read.csv(seasonality_file, sep = "\t")

colnames(cyto2)  <- gsub("-", ".", colnames(cyto2))
row.names(ses)   <- ses$ResearchID
pcs              <- pcs$x


# --- target cytokine ---------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No target cytokine provided.")

target <- gsub("-", ".", args[1])
cat("Analyzing target:", target, "\n")

if (!target %in% colnames(cyto2)) {
  stop("Target cytokine not found in cytokine matrix: ", target)
}


# --- keep the samples present in every table ---------------------------------

shared_row_names <- Reduce(intersect, list(
  rownames(cyto2),
  rownames(cov2),
  rownames(pcs),
  rownames(ses),
  rownames(db)
))

shared_row_names <- intersect(colnames(Mval), shared_row_names)

Mval  <- Mval[, shared_row_names, drop = FALSE]
cyto2 <- cyto2[shared_row_names, , drop = FALSE]
cov2  <- cov2[shared_row_names, , drop = FALSE]
pcs   <- pcs[shared_row_names, , drop = FALSE]
ses   <- ses[shared_row_names, , drop = FALSE]
db    <- db[shared_row_names, , drop = FALSE]

cat("Checking row compatibility...\n")
datasets <- list(cyto2, cov2, pcs, ses, db)
names(datasets) <- c("cyto2", "cov2", "pcs", "ses", "db")
for (name in names(datasets)) {
  cat(name, "rows:", nrow(datasets[[name]]), "\n")
}


# --- covariate table ---------------------------------------------------------
# BMI, smoking, the methylation PCs and the seasonality terms are assembled here
# but are not in the model formula below. They are kept because they take part
# in the sample intersection above, so dropping them would change the set of
# samples analysed, and because they are the obvious candidates for a
# sensitivity analysis.

data <- data.frame(
  cyto2[[target]],
  cov2$age,
  cov2$gender,
  cov2$bmi,
  cov2$smoke,
  as.numeric(as.factor(cov2$SAMPLE_PLATE)),
  pcs[, 1:6, drop = FALSE],
  db,
  ses[, c("season_sin", "season_cos", "numDaysFromJan2013")]
)

colnames(data) <- c(
  target, "age", "sex", "bmi", "smoke", "plate",
  paste0("PC", 1:6), colnames(db), "sin", "cos", "num"
)
row.names(data) <- shared_row_names

cat("Filtering data for target:", target, "\n")
data <- data[!is.na(data[[target]]), , drop = FALSE]


# --- methylation matrix ------------------------------------------------------
# Samples must end up in the rows of Mval, whichever way round the input was
# stored.

ensure_correct_orientation <- function(Mval, sample_names) {
  if (any(rownames(Mval) %in% sample_names)) {
    cat("Mval row names contain sample names. No transposition needed.\n")
    return(Mval)
  }

  cat("Transposing Mval to check for sample names in row names...\n")
  Mval <- t(Mval)
  if (any(rownames(Mval) %in% sample_names)) {
    cat("Mval row names contain sample names after transposition.\n")
    return(Mval)
  }

  stop("Row names of Mval do not match sample names, even after transposition.")
}

Mval <- ensure_correct_orientation(Mval, shared_row_names)

Mval <- Mval[shared_row_names, , drop = FALSE]
Mval <- Mval[rownames(data), , drop = FALSE]

if (!identical(rownames(Mval), rownames(data))) {
  stop("Row names of Mval and data do not match.")
}


# --- surrogate variables -----------------------------------------------------
# Estimated beforehand on the methylation data and reused here so that every
# cytokine is adjusted for the same set.

cat("Adding surrogate variables (SVA)...\n")

SVA <- as.data.frame(readRDS(sva_file))
SVA <- SVA[shared_row_names, , drop = FALSE]
SVA <- SVA[rownames(data), , drop = FALSE]

data <- cbind(data, SVA)

head(data)
dim(data)


# --- model -------------------------------------------------------------------
# Granulocytes are left out of the cell composition terms because the
# proportions sum to one and would otherwise be collinear.

surrogate_vars <- paste(colnames(SVA), collapse = " + ")

formula <- as.formula(paste(target,
                            "~ methylation + age + sex + CD8T + Mono + Bcell + CD4T + NK + plate +",
                            surrogate_vars))
print(formula)


# --- fit one model per CpG ---------------------------------------------------

registerDoParallel(cores = num_cores)
cat("Using", num_cores, "cores for parallel processing.\n")

num_batches <- ceiling(length(colnames(Mval)) / batch_size)
cat("Processing in", num_batches, "batches...\n")

results <- foreach(batch_idx = seq_len(num_batches),
                   .combine  = rbind,
                   .packages = c("MASS", "lmtest", "sandwich")) %dopar% {

  start_idx  <- (batch_idx - 1) * batch_size + 1
  end_idx    <- min(batch_idx * batch_size, length(colnames(Mval)))
  batch_cols <- colnames(Mval)[start_idx:end_idx]

  batch_results <- lapply(batch_cols, function(cpG) {
    model_data <- data.frame(methylation = Mval[, cpG], data)

    # A CpG that fails to converge is reported as NA rather than dropped
    result <- c(CpGsite = cpG, Estimate = NA, `Std. Error` = NA,
                `z value` = NA, `Pr(>|z|)` = NA)

    tryCatch({
      rlm_model <- rlm(formula, data = model_data, maxit = 200)

      # Heteroskedasticity-consistent standard errors for the methylation term
      cf <- coeftest(rlm_model, vcov = vcovHC(rlm_model, type = "HC0"))
      result <- c(CpGsite = cpG, cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
    }, error = function(e) {
      cat("Error processing CpG site:", cpG, "- Message:", conditionMessage(e), "\n")
    })

    return(result)
  })

  return(do.call(rbind, batch_results))
}

stopImplicitCluster()


# --- multiple testing and genomic inflation ----------------------------------
# FDR and BH are the same correction, both are written so that either column
# name can be used downstream.

results <- as.data.frame(results, stringsAsFactors = FALSE)

pvalues <- as.numeric(as.character(results$`Pr(>|z|)`))

results$FDR <- p.adjust(pvalues, method = "fdr")
results$BH  <- p.adjust(pvalues, method = "BH")
results$BY  <- p.adjust(pvalues, method = "BY")

results$lambda <- median(qchisq(pvalues, df = 1, lower.tail = FALSE), na.rm = TRUE) /
  qchisq(0.5, 1)


# --- write -------------------------------------------------------------------

output_file <- file.path(out_dir, paste0(target, "_results_model2.cleanSVA.txt"))
write.table(results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Results saved to:", output_file, "\n")
