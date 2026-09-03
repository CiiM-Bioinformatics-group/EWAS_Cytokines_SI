#!/usr/bin/env Rscript
# =============================================================================
# Longitudinal methylation remodeling at cytokine-associated CpGs
#
# The cCpGs tested here were identified in the SI cohort, which is not a
# vaccination cohort and was sampled at a single time point, so it cannot show
# whether these sites change over time. They are therefore carried into the
# independent 300BCG cohort, which received BCG vaccination and was sampled
# repeatedly, to ask whether cCpGs undergo more methylation remodeling than
# comparable CpGs elsewhere in the genome.
#
# In 300BCG, blood samples were collected before vaccination (day 0), 14 days
# after (day 14) and 3 months after (day 90). Only individuals with complete
# data at all three time points are used.
#
# For every CpG and every interval a per-person change score is formed
#
#     dM = M(later) - M(earlier)
#
# and modelled as a function of the change in estimated cell-type proportions
# and sample plate, with age, sex, BMI and smoking status as subject-level
# terms. The intercept of that model is the covariate-adjusted magnitude of
# methylation change for the CpG, and its absolute value carries into the
# enrichment step.
#
# Changes in cell proportion are mean-centred, so zero on those columns is the
# average observed change rather than "no change". The plate change indicators
# are left uncentred, so zero there does mean the sample stayed on the same
# plate.
#
# Enrichment: the mean absolute change across a CpG set is compared with the
# same quantity in matched background sets of equal size, drawn from the set's
# own pool. Backgrounds are matched on baseline methylation level and on the
# variability of methylation change within the interval. The ratio of the two
# is the enrichment; a ratio above one indicates increased remodeling. P values
# are empirical two-sided probabilities from the matched background
# distribution, adjusted across CpG sets within each interval with the
# Benjamini-Hochberg procedure.
#
# The smallest achievable empirical p is 1/(nboot+1).
#
# Usage
#   Rscript delta_longitudinal_methylation_change.R \
#     --mval=M_values_matrix.rds --pheno=phenotypes.rds \
#     --dir=data/cpg_sets --outdir=results/longitudinal_change \
#     --nboot=2000 --nbin=8 --dropplate=none
#
# Inputs
#   --mval      methylation M-values, CpGs in rows, samples in columns
#   --pheno     one row per sample, with columns: id, time (0/14/90), the cell
#               proportions, a plate column, age, sex, BMI, smoking
#   --dir       directory holding, for each CpG set <nm>:
#                 target_cpgs_<nm>.txt        the cytokine-associated CpGs
#                 iteration_sets_<nm>_LONG.tsv  the background pool, column cpg
# =============================================================================

argv <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NA_character_) {
  hit <- grep(paste0("^--", name, "="), argv, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[1])
}

MVAL     <- get_arg("mval")
PHENO    <- get_arg("pheno")
DIR      <- get_arg("dir", ".")
OUTDIR   <- get_arg("outdir", "longitudinal_change")
NBOOT    <- as.integer(get_arg("nboot", "2000"))
NBIN     <- as.integer(get_arg("nbin", "8"))
DROPPL   <- get_arg("dropplate", "none")
SEED     <- as.integer(get_arg("seed", "42"))

CELLTYPES <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")

set.seed(SEED)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# --- logging -----------------------------------------------------------------

log_lines <- character(0)
say <- function(...) {
  s <- paste0(...)
  log_lines <<- c(log_lines, s)
  cat(s, "\n", sep = "")
}
rule <- function(...) say(paste0("\n", strrep("=", 92), "\n", paste0(...), "\n", strrep("=", 92)))


# --- helpers -----------------------------------------------------------------

load_any <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    e <- new.env(); n <- load(path, envir = e); e[[n[1]]]
  }
}

mval_to_beta <- function(m) { x <- 2^m; x / (1 + x) }

# Row variances, tolerating missing values
row_var <- function(x) {
  n  <- rowSums(!is.na(x))
  mu <- rowMeans(x, na.rm = TRUE)
  v  <- rowSums((x - mu)^2, na.rm = TRUE) / (n - 1)
  v[n < 2] <- NA
  v
}

# Two-sided empirical p: how often the background deviates from its own centre
# at least as far as the observed value does
empirical_p <- function(observed, null_values) {
  null_values <- null_values[is.finite(null_values)]
  if (!length(null_values)) return(NA_real_)
  centre <- mean(null_values)
  (1 + sum(abs(null_values - centre) >= abs(observed - centre))) / (1 + length(null_values))
}


# =============================================================================
# 1. DATA
# =============================================================================
rule("1. DATA")

pd <- as.data.frame(load_any(PHENO))
pd$sample <- rownames(pd)
M <- as.matrix(load_any(MVAL))

PLATE_COL <- if ("Sample_Plate" %in% colnames(pd)) "Sample_Plate" else
  grep("plate", colnames(pd), ignore.case = TRUE, value = TRUE)[1]

say("plate column : ", PLATE_COL)
say("plate levels : ", paste(names(table(as.character(pd[[PLATE_COL]]))), collapse = ", "))

if (DROPPL != "none") {
  keep <- as.character(pd[[PLATE_COL]]) != DROPPL
  if (!sum(!keep)) say("  !! --dropplate=", DROPPL, " matched nothing, check the spelling")
  say("  dropping ", sum(!keep), " samples on ", DROPPL)
  pd <- pd[keep, ]
}

# CpG sets: the cytokine-associated targets and the pool they are matched against
read_targets <- function(nm) {
  f <- file.path(DIR, paste0("target_cpgs_", nm, ".txt"))
  if (!file.exists(f)) character(0) else unique(trimws(readLines(f, warn = FALSE)))
}
read_pool <- function(nm) {
  f <- file.path(DIR, paste0("iteration_sets_", nm, "_LONG.tsv"))
  if (!file.exists(f)) character(0) else unique(read.delim(f, stringsAsFactors = FALSE)$cpg)
}

SETS <- list(
  ALL   = list(t = read_targets("ALL_all"),   p = read_pool("ALL_all")),
  lps   = list(t = read_targets("lps_all"),   p = read_pool("lps_all")),
  pam3  = list(t = read_targets("pam3_all"),  p = read_pool("pam3_all")),
  covid = list(t = read_targets("covid_all"), p = read_pool("covid_all"))
)
SETS$innate <- list(t = union(SETS$lps$t, SETS$pam3$t),
                    p = union(SETS$lps$p, SETS$pam3$p))

say("")
say(sprintf("  %-8s %7s %9s", "set", "target", "pool"))
for (nm in names(SETS)) {
  say(sprintf("  %-8s %7d %9d", nm, length(SETS[[nm]]$t), length(SETS[[nm]]$p)))
}
say("  innate = union(lps, pam3), duplicates removed, own pooled background.")

# Keep only the CpGs any set refers to, and only samples present in both tables
need <- intersect(unique(unlist(lapply(SETS, function(s) c(s$t, s$p)), use.names = FALSE)),
                  rownames(M))
M  <- M[need, intersect(colnames(M), pd$sample), drop = FALSE]
invisible(gc())
pd <- pd[pd$sample %in% colnames(M), ]

# Complete trios only
pd <- pd[pd$time %in% c(0, 14, 90), ]
subject_counts <- table(pd$id)
complete <- names(subject_counts)[subject_counts == 3]

say("")
say("subjects with all three timepoints : ", length(complete), " of ", length(subject_counts))

pd <- pd[pd$id %in% complete, ]
pd <- pd[order(pd$id, pd$time), ]
M  <- M[, pd$sample, drop = FALSE]
invisible(gc())

p0  <- pd[pd$time == 0,  ]
p14 <- pd[pd$time == 14, ]
p90 <- pd[pd$time == 90, ]
stopifnot(identical(as.character(p0$id), as.character(p14$id)),
          identical(as.character(p0$id), as.character(p90$id)))

s0  <- p0$sample
s14 <- p14$sample

say("  timepoint blocks person-aligned : OK (", length(s0), " subjects each)")
say("  CpGs loaded : ", length(need))


# =============================================================================
# 2. CHANGE-SCORE MODEL
# =============================================================================

BIOVARS <- grep("^age$|^sex$|^gender$|smok|^bmi$", colnames(pd),
                ignore.case = TRUE, value = TRUE)

rule("2. FIT   dM ~ dCells + dPlate + age + sex + BMI + smoking")
say("  intercept = adjusted mean methylation change over the interval.")
say("  subject-level covariates : ", paste(BIOVARS, collapse = ", "))

fit_window <- function(t1, t2) {

  a <- pd[pd$time == t1, ]
  b <- pd[pd$time == t2, ]
  stopifnot(identical(as.character(a$id), as.character(b$id)))

  # Change in cell proportions, centred. The most abundant cell type at the
  # earlier visit is dropped as reference because the proportions are compositional.
  d_cells   <- as.matrix(b[, CELLTYPES]) - as.matrix(a[, CELLTYPES])
  reference <- CELLTYPES[which.max(colMeans(as.matrix(a[, CELLTYPES]), na.rm = TRUE))]
  X <- cbind(`(Intercept)` = 1,
             scale(d_cells[, setdiff(CELLTYPES, reference), drop = FALSE], TRUE, FALSE))

  # Change of plate, as uncentred indicator differences
  levels_all <- sort(unique(c(as.character(a[[PLATE_COL]]), as.character(b[[PLATE_COL]]))))
  indicator <- function(z) {
    Z <- matrix(0, length(z), length(levels_all), dimnames = list(NULL, levels_all))
    Z[cbind(seq_along(z), match(as.character(z), levels_all))] <- 1
    Z
  }
  d_plate <- (indicator(b[[PLATE_COL]]) - indicator(a[[PLATE_COL]]))[, -1, drop = FALSE]
  if (ncol(d_plate)) {
    d_plate <- d_plate[, apply(d_plate, 2, function(z) stats::var(z) > 0), drop = FALSE]
  }
  n_changed_plate <- sum(as.character(a[[PLATE_COL]]) != as.character(b[[PLATE_COL]]))
  if (ncol(d_plate)) X <- cbind(X, d_plate)

  # Subject-constant covariates enter as covariate x time, i.e. whether these
  # participants remodel at a different rate
  for (v in BIOVARS) {
    x <- a[[v]]
    if (is.numeric(x)) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      X <- cbind(X, matrix(scale(x, TRUE, FALSE), ncol = 1, dimnames = list(NULL, v)))
    } else {
      f <- factor(as.character(x))
      if (nlevels(droplevels(f)) < 2) next
      mm <- stats::model.matrix(~ droplevels(f))[, -1, drop = FALSE]
      colnames(mm) <- paste0(v, "_", seq_len(ncol(mm)))
      X <- cbind(X, scale(mm, TRUE, FALSE))
    }
  }

  # Drop any columns that are rank deficient
  q <- qr(X)
  X <- X[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  XtX_inv <- solve(crossprod(X))

  # One least-squares fit for every CpG at once
  D  <- M[, b$sample, drop = FALSE] - M[, a$sample, drop = FALSE]
  ok <- rowSums(is.na(D)) == 0
  Y  <- t(D[ok, , drop = FALSE])
  B  <- XtX_inv %*% crossprod(X, Y)

  resid_var <- colSums((Y - X %*% B)^2) / (nrow(X) - ncol(X))

  estimate <- rep(NA_real_, nrow(D))
  std_err  <- estimate
  estimate[ok] <- as.numeric(B[1, ])
  std_err[ok]  <- sqrt(resid_var * XtX_inv[1, 1])
  names(estimate) <- rownames(D)
  names(std_err)  <- rownames(D)

  # Variability of the change itself, one of the two matching axes
  sd_change <- sqrt(row_var(D))

  rm(D, Y, B); invisible(gc())

  list(est = estimate, se = std_err, n = nrow(a),
       k = ncol(X) - 1, nchg = n_changed_plate, sdd = sd_change)
}

WINDOWS <- list(`d0->d14` = c(0, 14), `d14->d90` = c(14, 90), `d0->d90` = c(0, 90))

FIT <- list()
for (w in names(WINDOWS)) {
  FIT[[w]] <- fit_window(WINDOWS[[w]][1], WINDOWS[[w]][2])
  say(sprintf("  %-10s pairs %d   covariates %d   resid df %d   changed plate %d",
              w, FIT[[w]]$n, FIT[[w]]$k, FIT[[w]]$n - FIT[[w]]$k - 1, FIT[[w]]$nchg))
}


# =============================================================================
# 3. ENRICHMENT AGAINST MATCHED BACKGROUND
# =============================================================================
rule("3. ENRICHMENT   mean |estimate| in target vs matched background")

# Matching axis 1: baseline methylation as beta, at the earlier visit
baseline_d0  <- rowMeans(mval_to_beta(M[, s0,  drop = FALSE]), na.rm = TRUE)
baseline_d14 <- rowMeans(mval_to_beta(M[, s14, drop = FALSE]), na.rm = TRUE)

# Matching axis 2: variability of methylation change within the interval
AXES <- list(
  `d0->d14`  = list(baseline_d0,  FIT[["d0->d14"]]$sdd),
  `d14->d90` = list(baseline_d14, FIT[["d14->d90"]]$sdd),
  `d0->d90`  = list(baseline_d0,  FIT[["d0->d90"]]$sdd)
)

rm(M); invisible(gc())

INDEX <- setNames(seq_along(need), need)

# Draw nboot background sets that reproduce the target's joint distribution
# across the matching axes. Each axis is cut into NBIN quantile bins defined on
# the pool; a target CpG is replaced by pool CpGs from the same bin combination,
# or from the nearest occupied one when that combination is empty.
make_background <- function(target_ix, pool_ix, axes, nboot) {

  edges <- function(q) { q <- unique(q); q[1] <- -Inf; q[length(q)] <- Inf; q }

  bin_key <- function(ix) {
    do.call(paste, c(lapply(axes, function(v) {
      cut(v[ix],
          edges(stats::quantile(v[pool_ix], seq(0, 1, length.out = NBIN + 1), na.rm = TRUE)),
          include.lowest = TRUE, labels = FALSE)
    }), sep = ":"))
  }

  needed_per_bin <- table(bin_key(target_ix))
  pool_by_bin    <- split(pool_ix, bin_key(pool_ix))
  occupied       <- names(pool_by_bin)
  if (!length(occupied)) return(NULL)

  occupied_coords <- do.call(rbind, lapply(occupied, function(k) {
    as.integer(strsplit(k, ":")[[1]])
  }))

  source_bin <- setNames(vector("list", length(needed_per_bin)), names(needed_per_bin))
  for (k in names(needed_per_bin)) {
    source_bin[[k]] <- if (!is.null(pool_by_bin[[k]]) && length(pool_by_bin[[k]])) {
      pool_by_bin[[k]]
    } else {
      nearest <- which.min(rowSums(abs(sweep(occupied_coords, 2,
                                             as.integer(strsplit(k, ":")[[1]])))))
      pool_by_bin[[occupied[nearest]]]
    }
  }

  replicate(nboot,
            unlist(lapply(names(needed_per_bin), function(k) {
              p <- source_bin[[k]]
              if (!length(p)) return(integer(0))
              sample(p, needed_per_bin[[k]], replace = length(p) < needed_per_bin[[k]])
            }), use.names = FALSE),
            simplify = FALSE)
}

say("")
say(sprintf("  %-8s %-10s %6s %10s %8s %11s %10s %7s",
            "set", "window", "n", "obs", "ratio", "P", "base_chk", "sd_chk"))
say("  ", strrep("-", 82))

RES <- NULL
for (nm in names(SETS)) {

  target_ix <- unname(INDEX[intersect(SETS[[nm]]$t, need)])
  pool_ix   <- unname(INDEX[setdiff(intersect(SETS[[nm]]$p, need), SETS[[nm]]$t)])

  for (w in names(AXES)) {

    est <- FIT[[w]]$est
    tk  <- target_ix[is.finite(est[target_ix])]
    pk  <- pool_ix[is.finite(est[pool_ix])]
    if (length(tk) < 30 || length(pk) < 100) next

    draws <- make_background(tk, pk, AXES[[w]], NBOOT)
    if (is.null(draws)) next

    observed <- mean(abs(est[tk]))
    null_vals <- vapply(draws, function(s) mean(abs(est[s]), na.rm = TRUE), numeric(1))

    # Confirmation that the matching worked: both should be close to 1
    base_chk <- mean(AXES[[w]][[1]][tk]) /
      mean(vapply(draws, function(s) mean(AXES[[w]][[1]][s], na.rm = TRUE), numeric(1)))
    sd_chk <- mean(AXES[[w]][[2]][tk], na.rm = TRUE) /
      mean(vapply(draws, function(s) mean(AXES[[w]][[2]][s], na.rm = TRUE), numeric(1)), na.rm = TRUE)

    pval <- empirical_p(observed, null_vals)

    say(sprintf("  %-8s %-10s %6d %10.5f %8.3f %11.5f %10.3f %7.3f",
                nm, w, length(tk), observed, observed / mean(null_vals), pval, base_chk, sd_chk))

    RES <- rbind(RES, data.frame(
      set = nm, window = w, n = length(tk),
      observed = observed, null = mean(null_vals),
      ratio = observed / mean(null_vals), p_value = pval,
      base_chk = base_chk, sd_chk = sd_chk,
      stringsAsFactors = FALSE))

    rm(draws); invisible(gc())
  }
  say("")
}

# Sets are corrected within an interval. Intervals are not pooled, they answer
# three different questions.
RES$fdr <- ave(RES$p_value, RES$window, FUN = function(x) stats::p.adjust(x, "BH"))

write.csv(RES, file.path(OUTDIR, "MAIN_TABLE.csv"), row.names = FALSE)

for (w in names(FIT)) stopifnot(identical(names(FIT[[w]]$est), need))
saveRDS(c(lapply(FIT, function(f) list(est = f$est, se = f$se)), list(.cpg = need)),
        file.path(OUTDIR, "delta_estimates.rds"))


# =============================================================================
# 4. RESULT TABLE
# =============================================================================
rule("4. RESULT TABLE")

labels <- c(ALL = "ALL", lps = "LPS", pam3 = "Pam3", covid = "COVID",
            innate = "LPS+Pam3 innate")

stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else
    if (p < 0.05) "*" else ""
}

say(sprintf("  %-16s %6s | %8s %9s %4s | %8s %9s %4s | %8s %9s %4s",
            "Set", "n", "0->14", "P", "", "14->90", "P", "", "0->90", "P", ""))
say("  ", strrep("-", 88))

for (nm in names(labels)) {
  r <- RES[RES$set == nm, ]
  g <- function(w, f) { z <- r[r$window == w, ]; if (!nrow(z)) NA_real_ else z[[f]][1] }
  say(sprintf("  %-16s %6s | %8s %9s %4s | %8s %9s %4s | %8s %9s %4s", labels[nm],
              ifelse(nrow(r), format(r$n[1], big.mark = ","), "-"),
              sprintf("%.3f", g("d0->d14",  "ratio")), sprintf("%.5f", g("d0->d14",  "p_value")), stars(g("d0->d14",  "p_value")),
              sprintf("%.3f", g("d14->d90", "ratio")), sprintf("%.5f", g("d14->d90", "p_value")), stars(g("d14->d90", "p_value")),
              sprintf("%.3f", g("d0->d90",  "ratio")), sprintf("%.5f", g("d0->d90",  "p_value")), stars(g("d0->d90",  "p_value"))))
}

say("")
say("  * p<0.05   ** p<0.01   *** p<0.001")
say(sprintf("  P floor = 1/(nboot+1) = %.1e", 1 / (NBOOT + 1)))
say("  BH-adjusted p across the sets within each interval : `fdr` in MAIN_TABLE.csv")
say("  base_chk, sd_chk : ratio of target to background on each matching axis")

writeLines(log_lines, file.path(OUTDIR, "log.txt"))
say("\nwritten to ", OUTDIR)
