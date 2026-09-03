#!/usr/bin/env Rscript
# =============================================================================
# Background CpG sets for the longitudinal methylation change analysis
#
# Builds, for each cCpG group, a pool of background CpGs that reproduces the
# group's own chromosome x genomic annotation composition.
#
# Background pool
#   The full array annotation, after removing every cCpG and every CpG within
#   30 kb of a cCpG. Neighbouring CpGs are correlated, so leaving them in would
#   blunt the contrast between target and background.
#
# Matching
#   Draws are made stratum by stratum, a stratum being a chromosome x
#   annotation type combination, so each group's background carries the same
#   genomic composition as the group itself. Every group gets its own strata
#   table and its own draws: a background matched to all cCpGs is not a valid
#   background for the smaller single-stimulus subsets, which have a different
#   composition.
#
#   A CpG can belong to more than one annotation type, so the per-stratum draws
#   sum to more than the target size. The surplus is removed at random rather
#   than by taking the first n, which would cut in chromosome order and
#   under-represent the last chromosomes.
#
#   Matching on baseline methylation level and on the variability of
#   methylation change is applied downstream, in
#   delta_longitudinal_methylation_change.R, which draws from these pools.
#
# Groups
#   ALL_all     every cCpG
#   lps_all     cCpGs from 24h LPS stimulated cytokines
#   pam3_all    cCpGs from 24h Pam3Cys stimulated cytokines
#   covid_all   cCpGs from 7d SARS-CoV-2 stimulated cytokines
#
# Usage
#   Rscript build_iteration_sets.R --niter=100 --outdir=data/cpg_sets
#
# Outputs, per group
#   target_cpgs_<group>.txt          the cCpGs in the group
#   iteration_sets_<group>_LONG.tsv  the matched background draws
# and overall
#   target_group_composition.tsv     annotation composition of each group
#   iteration_sets_SUMMARY.tsv       draw sizes and strata coverage
# =============================================================================

argv <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), argv, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[1])
}

N_ITER <- as.integer(get_arg("niter", "100"))
FLANK  <- as.numeric(get_arg("flank", "30000"))
MINSET <- as.integer(get_arg("minset", "50"))
SEED   <- as.integer(get_arg("seed", "42"))
OUTDIR <- get_arg("outdir", "data/cpg_sets")

CCPG_LIST <- get_arg("cpgs",    "data/cpgs_all_SI.txt")
CC_ANNOT  <- get_arg("ccannot", "data/ccpg_annotation.rds")
BG_ANNOT  <- get_arg("bgannot", "data/array_annotation.rds")
SIGDIR    <- get_arg("sigdir",  "results/EWAS_discovery_SI")

set.seed(SEED)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

say  <- function(...) cat(paste0(...), "\n", sep = "")
rule <- function(...) cat("\n", strrep("=", 78), "\n", paste0(...), "\n", strrep("=", 78), "\n", sep = "")

# Annotation names carry a trailing array-design suffix that the CpG lists do not
strip <- function(x) sub("_[A-Z]{2}[0-9]{2}$", "", x)
KEY   <- function(a, b) paste(a, b, sep = "\r")


# =============================================================================
rule("1. cCpGs AND THEIR ANNOTATION")

ccpgs <- unique(strip(read.table(CCPG_LIST, header = FALSE, stringsAsFactors = FALSE)$V1))
say("cCpGs in list : ", length(ccpgs))

cc <- as.data.frame(readRDS(CC_ANNOT))
cc$name_clean <- strip(cc$name)
cc <- cc[cc$name_clean %in% ccpgs, ]

n_ann <- length(unique(cc$name_clean))
say("annotation rows for cCpGs : ", nrow(cc), "   unique cCpGs annotated : ", n_ann)
if (n_ann < length(ccpgs)) {
  say("  !! ", length(ccpgs) - n_ann, " cCpGs have no annotation row and cannot be",
      " matched; they are excluded from every group.")
}


# =============================================================================
rule("2. BACKGROUND POOL  (cCpGs and +/-", FLANK / 1000, " kb around them excluded)")

bgA <- as.data.frame(readRDS(BG_ANNOT))
bgA$name_clean <- strip(bgA$name)
say("background annotation rows : ", nrow(bgA),
    "   unique CpGs : ", length(unique(bgA$name_clean)))

poscol <- intersect(c("start", "pos", "position"), colnames(cc))[1]
if (is.na(poscol)) stop("no position column in the cCpG annotation - cannot apply the flank exclusion")
say("position column used : ", poscol)

cc_pos <- unique(cc[, c("seqnames", poscol, "name_clean")]);  names(cc_pos)[2] <- "pos"
bg_pos <- unique(bgA[, c("seqnames", poscol, "name_clean")]); names(bg_pos)[2] <- "pos"

# Distance to the nearest cCpG, per chromosome, by sorted search
near <- unlist(lapply(split(seq_len(nrow(bg_pos)), as.character(bg_pos$seqnames)), function(ix) {
  B <- bg_pos[ix, ]
  a <- sort(cc_pos$pos[as.character(cc_pos$seqnames) == as.character(B$seqnames[1])])
  if (!length(a)) return(character(0))
  i    <- findInterval(B$pos, a)
  d_lo <- ifelse(i >= 1,         abs(B$pos - a[pmax(i, 1)]),          Inf)
  d_hi <- ifelse(i < length(a),  abs(a[pmin(i + 1, length(a))] - B$pos), Inf)
  B$name_clean[pmin(d_lo, d_hi) <= FLANK]
}), use.names = FALSE)
near <- unique(near)
say("background CpGs within ", FLANK / 1000, " kb of a cCpG : ", length(near))

bg <- unique(bgA[, c("seqnames", "annot.type", "name_clean")])
bg <- bg[!bg$name_clean %in% ccpgs & !bg$name_clean %in% near, ]
say("background pool : ", length(unique(bg$name_clean)), " CpGs across ", nrow(bg), " strata rows")

bg_by <- split(bg$name_clean, KEY(bg$seqnames, bg$annot.type))


# =============================================================================
rule("3. STIMULUS GROUPS")

sig_files <- list.files(SIGDIR, pattern = "\\.sig\\.FDR\\.txt$", full.names = TRUE)

memb <- do.call(rbind, lapply(sig_files, function(f) {
  b <- sub("\\.sig\\.FDR\\.txt$", "", sub("^pbmc_", "", basename(f)))
  p <- strsplit(b, "_", fixed = TRUE)[[1]]
  if (length(p) < 3) return(NULL)
  x <- tryCatch(read.table(f, header = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(x)) return(NULL)
  data.frame(stim = paste0(p[1], "_", paste(p[-(1:2)], collapse = "_")),
             cpg  = strip(trimws(as.character(x[[1]]))),
             stringsAsFactors = FALSE)
}))
memb <- unique(memb)
say("stimulus labels found : ", paste(sort(unique(memb$stim)), collapse = ", "))

pick <- function(pat) unique(memb$cpg[grepl(pat, memb$stim, ignore.case = TRUE)])

STIM <- list(ALL   = ccpgs,
             lps   = intersect(pick("lps"),  ccpgs),
             pam3  = intersect(pick("pam3"), ccpgs),
             covid = intersect(pick("cov"),  ccpgs))

for (nm in names(STIM)) say(sprintf("  %-6s %5d cCpGs", nm, length(STIM[[nm]])))


# =============================================================================
rule("4. TARGET GROUPS")

TARGETS <- list()
for (s in names(STIM)) {
  TARGETS[[paste0(s, "_all")]] <- intersect(STIM[[s]], cc$name_clean)
}

say(sprintf("  %-16s %7s", "group", "n"))
for (nm in names(TARGETS)) {
  say(sprintf("  %-16s %7d%s", nm, length(TARGETS[[nm]]),
              if (length(TARGETS[[nm]]) < MINSET) "   <- below --minset, skipped" else ""))
}
TARGETS <- TARGETS[vapply(TARGETS, length, integer(1)) >= MINSET]

for (nm in names(TARGETS)) {
  writeLines(TARGETS[[nm]], file.path(OUTDIR, paste0("target_cpgs_", nm, ".txt")))
}
say("target CpG lists written : target_cpgs_<group>.txt")

# Annotation composition of each group, the table that justifies giving every
# group its own background
compo <- do.call(rbind, lapply(names(TARGETS), function(nm) {
  d <- unique(cc[cc$name_clean %in% TARGETS[[nm]], c("annot.type", "name_clean")])
  n <- length(unique(d$name_clean))
  p <- table(d$annot.type)
  data.frame(group = nm, n_cpg = n, annot.type = names(p),
             pct = round(100 * as.integer(p) / n, 1), stringsAsFactors = FALSE)
}))
write.table(compo, file.path(OUTDIR, "target_group_composition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("target composition written : target_group_composition.tsv")


# =============================================================================
rule("5. GENERATING ", N_ITER, " MATCHED SETS PER GROUP")

strata_of <- function(cpgs) {
  d  <- unique(cc[cc$name_clean %in% cpgs, c("seqnames", "annot.type", "name_clean")])
  tb <- table(KEY(d$seqnames, d$annot.type))
  list(key = names(tb), n = as.integer(tb))
}

SUM <- list()
for (nm in names(TARGETS)) {

  tg       <- TARGETS[[nm]]
  st       <- strata_of(tg)
  target_n <- length(tg)

  # Strata where the pool holds fewer CpGs than the target needs stay under-filled
  avail <- vapply(st$key, function(k) length(bg_by[[k]]), integer(1))
  short <- sum(avail < st$n)

  out <- vector("list", N_ITER)
  for (it in seq_len(N_ITER)) {

    drawn <- unlist(lapply(seq_along(st$key), function(i) {
      pool <- bg_by[[st$key[i]]]
      if (is.null(pool) || !length(pool)) return(character(0))
      sample(pool, min(st$n[i], length(pool)), replace = FALSE)
    }), use.names = FALSE)

    drawn <- unique(drawn)
    if (length(drawn) > target_n) drawn <- sample(drawn, target_n)

    out[[it]] <- data.frame(Iteration = it, slot = seq_along(drawn), cpg = drawn,
                            stringsAsFactors = FALSE)
  }

  L <- do.call(rbind, out)
  f <- file.path(OUTDIR, paste0("iteration_sets_", nm, "_LONG.tsv"))
  write.table(L, f, sep = "\t", quote = FALSE, row.names = FALSE)

  sz <- as.integer(table(L$Iteration))
  say(sprintf("  %-16s target %5d | drawn %5d/%5d/%5d | strata %3d (%d short) -> %s",
              nm, target_n, min(sz), as.integer(stats::median(sz)), max(sz),
              length(st$key), short, basename(f)))

  SUM[[nm]] <- data.frame(group = nm, n_target = target_n,
                          n_drawn_min = min(sz), n_drawn_median = stats::median(sz),
                          n_strata = length(st$key), n_strata_short = short,
                          file = basename(f), stringsAsFactors = FALSE)
}

S <- do.call(rbind, SUM)
write.table(S, file.path(OUTDIR, "iteration_sets_SUMMARY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

say("")
say("written to ", OUTDIR)
