# =============================================================================
# Immune-trait enrichment of cCpGs in the EWAS Catalog and the GWAS Catalog
#
# The same procedure is applied to both catalogs. For every trait, the number
# of cCpGs associated with that trait is compared with the number obtained in
# the 100 matched background sets. The mean across those sets is the expected
# value under the null, and enrichment is tested with Fisher's exact test.
# Traits with FDR < 0.05 and an odds ratio above one are taken as enriched.
#
# The two catalogs differ only in how a CpG is linked to a trait:
#   EWAS Catalog   direct CpG-trait associations at p < 1e-5, so the same
#                  reference serves as both the target and background lookup
#   GWAS Catalog   a CpG is linked to a trait when its surrounding region
#                  carries a cluster of trait-associated SNPs. The two lookup
#                  files are prepared beforehand and already carry that filter
#
# A trait is only tested when at least three cCpGs and at least three
# background CpGs are linked to it, and only reported when both the observed
# and the expected count reach three. Below that the odds ratio is too unstable
# to interpret.
#
# Input
#   iteration_cpg_sets_LONG.tsv                  matched background sets
#   cpgs_all_SI.txt                              the cCpGs
#   EWAS_catalogue_BLOOD_IMMUNE_p1e5_FINAL.tsv   EWAS Catalog reference
#   sig_cpg_trait_min3snps.tsv                   GWAS Catalog, cCpG regions
#   background_cpg_trait_min3snps.tsv            GWAS Catalog, background regions
#
# Output
#   EWAS_enrichment_results.tsv, GWAS_enrichment_results.tsv
#   EWAS_panel_counts.png, EWAS_panel_enrichment.png
#   GWAS_panel_counts.png, GWAS_panel_enrichment.png
#   GLOBAL_category_legend.png / .pdf
# =============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

data_dir <- "data"
out_dir  <- "results/EWAS_GWAS_enrichment"
gwas_dir <- file.path(data_dir, "GWAS_catalog")

ITERATION_SETS_PATH      <- file.path(out_dir,  "iteration_cpg_sets_LONG.tsv")
CCPG_LIST_PATH           <- file.path(data_dir, "cpgs_all_SI.txt")
EWAS_REF_PATH            <- file.path(out_dir,  "EWAS_catalogue_BLOOD_IMMUNE_p1e5_FINAL.tsv")
GWAS_SIG_REF_PATH        <- file.path(gwas_dir, "sig_cpg_trait_min3snps.tsv")
GWAS_BACKGROUND_REF_PATH <- file.path(gwas_dir, "background_cpg_trait_min3snps.tsv")

MIN_OBSERVED_FOR_RESULT <- 3
MIN_EXPECTED_FOR_RESULT <- 3

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# One palette shared by both count panels and the legend, so a category keeps
# the same colour across the EWAS and GWAS figures
GLOBAL_CATEGORY_COLORS <- c(
  "Autoimmune / Inflammatory Disease" = "#2E6C7E",
  "Cancer"                            = "#D9705A",
  "Environmental Exposure"            = "#4A85AD",
  "Inflammation Marker"               = "#2B2E1E",
  "Infectious Disease"                = "#6E7F9E",
  "Chronic Inflammatory Disease"      = "#7FA6B3",
  "Cytokine / Protein Level"          = "#B8B992",
  "Immune Cell Count"                 = "#8C6BAE",
  "Immune-Related Disorder"           = "#C9A87C",
  "Other Immune-Related"              = "#A9A9A9"
)


# --- shared inputs -----------------------------------------------------------

ccpgs_raw   <- fread(CCPG_LIST_PATH, header = FALSE)$V1
ccpgs_clean <- unique(sub("_[A-Z]{2}[0-9]{2}$", "", ccpgs_raw))
n_ccpg      <- length(ccpgs_clean)
cat("cCpGs:", n_ccpg, "\n")

iterations_long <- fread(ITERATION_SETS_PATH)
n_iterations    <- length(unique(iterations_long$Iteration))
iteration_sets  <- split(iterations_long$cpg, iterations_long$Iteration)
cat("Matched background sets loaded:", n_iterations, "\n")


# --- enrichment --------------------------------------------------------------
# sig_lookup and background_lookup are named lists mapping a trait to its CpGs.
# For the EWAS Catalog both arguments are the same list.

run_enrichment <- function(sig_lookup, background_lookup, label) {

  all_traits <- unique(c(names(sig_lookup), names(background_lookup)))

  n_sig <- sapply(all_traits, function(tr) length(sig_lookup[[tr]]))
  n_bg  <- sapply(all_traits, function(tr) length(background_lookup[[tr]]))

  traits_to_test <- all_traits[n_sig >= MIN_OBSERVED_FOR_RESULT &
                                 n_bg >= MIN_OBSERVED_FOR_RESULT]

  cat("[", label, "] Traits before filtering:", length(all_traits), "\n")
  cat("[", label, "] Traits retained:", length(traits_to_test), "\n")

  null_overlap_matrix <- matrix(
    NA_integer_, nrow = n_iterations, ncol = length(traits_to_test),
    dimnames = list(names(iteration_sets), traits_to_test)
  )

  cat("[", label, "] Computing null overlap across", n_iterations, "sets...\n")

  for (iter_name in names(iteration_sets)) {
    null_set <- iteration_sets[[iter_name]]
    for (tr in traits_to_test) {
      null_overlap_matrix[iter_name, tr] <- length(intersect(null_set, background_lookup[[tr]]))
    }
  }

  results <- rbindlist(lapply(traits_to_test, function(tr) {

    observed        <- length(intersect(ccpgs_clean, sig_lookup[[tr]]))
    null_dist       <- null_overlap_matrix[, tr]
    expected        <- mean(null_dist)
    fold_enrichment <- ifelse(expected > 0, observed / expected, NA)

    p_enrich    <- (sum(null_dist >= observed) + 1) / (length(null_dist) + 1)
    p_deplete   <- (sum(null_dist <= observed) + 1) / (length(null_dist) + 1)
    empirical_p <- min(1, 2 * min(p_enrich, p_deplete))

    direction <- ifelse(observed > expected, "enriched",
                        ifelse(observed < expected, "depleted", "neutral"))

    expected_rounded <- round(expected)
    fisher_mat <- matrix(c(observed,         n_ccpg - observed,
                           expected_rounded, n_ccpg - expected_rounded), nrow = 2)
    fisher_res <- tryCatch(fisher.test(fisher_mat), error = function(e) NULL)

    data.table(
      trait              = tr,
      observed           = observed,
      expected_mean_null = round(expected, 2),
      fold_enrichment    = round(fold_enrichment, 2),
      direction          = direction,
      empirical_p        = empirical_p,
      fisher_OR          = if (!is.null(fisher_res)) round(fisher_res$estimate, 2)    else NA,
      fisher_CI_low      = if (!is.null(fisher_res)) round(fisher_res$conf.int[1], 2) else NA,
      fisher_CI_high     = if (!is.null(fisher_res)) round(fisher_res$conf.int[2], 2) else NA,
      fisher_p           = if (!is.null(fisher_res)) fisher_res$p.value               else NA
    )
  }))

  results[, empirical_p_FDR := p.adjust(empirical_p, method = "BH")]
  results[, fisher_p_FDR    := p.adjust(fisher_p,    method = "BH")]
  results <- results[order(fisher_p)]

  fwrite(results, file.path(out_dir, paste0(label, "_enrichment_results.tsv")), sep = "\t")

  cat("[", label, "] Done.\n\n")
  results
}


# --- EWAS Catalog ------------------------------------------------------------

ewas_ref    <- fread(EWAS_REF_PATH)
ewas_traits <- unique(ewas_ref$trait_canonical)
ewas_lookup <- lapply(ewas_traits, function(tr) unique(ewas_ref[trait_canonical == tr]$cpg))
names(ewas_lookup) <- ewas_traits

ewas_results <- run_enrichment(ewas_lookup, ewas_lookup, "EWAS")


# --- GWAS Catalog ------------------------------------------------------------

gwas_sig_ref        <- fread(GWAS_SIG_REF_PATH)
gwas_background_ref <- fread(GWAS_BACKGROUND_REF_PATH)

gwas_sig_traits <- unique(gwas_sig_ref$trait)
gwas_sig_lookup <- lapply(gwas_sig_traits, function(tr) unique(gwas_sig_ref[trait == tr]$cpg))
names(gwas_sig_lookup) <- gwas_sig_traits

gwas_bg_traits        <- unique(gwas_background_ref$trait)
gwas_background_lookup <- lapply(gwas_bg_traits, function(tr) unique(gwas_background_ref[trait == tr]$cpg))
names(gwas_background_lookup) <- gwas_bg_traits

gwas_results <- run_enrichment(gwas_sig_lookup, gwas_background_lookup, "GWAS")


# --- trait display names and categories --------------------------------------

trait_display_map <- c(
  "Monocyte_count" = "Monocyte Count", "Lymphocyte_count" = "Lymphocyte Count",
  "Neutrophil_count" = "Neutrophil Count", "White_blood_cell_count" = "White Blood Cell Count",
  "Eosinophil_count" = "Eosinophil Count", "Prostate_cancer" = "Prostate Cancer",
  "Breast_cancer" = "Breast Cancer", "Colorectal_cancer" = "Colorectal Cancer",
  "Lung_cancer" = "Lung Cancer", "Inflammatory_bowel_disease" = "Inflammatory Bowel Disease",
  "Inflammatory_bowel_disease_(MTAG)" = "Inflammatory Bowel Disease",
  "Asthma" = "Asthma", "Asthma_(childhood_onset)" = "Asthma (Childhood Onset)",
  "Crohn's_disease" = "Crohn's Disease", "Ulcerative_colitis" = "Ulcerative Colitis",
  "Psoriasis" = "Psoriasis", "Hypothyroidism" = "Hypothyroidism",
  "Rheumatoid_arthritis" = "Rheumatoid Arthritis", "Multiple_sclerosis" = "Multiple Sclerosis",
  "Systemic_lupus_erythematosus" = "Systemic Lupus Erythematosus",
  "Systemic_lupus_erythematosus_(MTAG)" = "Systemic Lupus Erythematosus",
  "Atopic_dermatitis" = "Atopic Dermatitis",
  "C-reactive_protein_levels_(MTAG)" = "C-Reactive Protein (CRP) Levels",
  "C-reactive protein (CRP) levels" = "C-Reactive Protein (CRP) Levels",
  "COVID-19_(hospitalized_covid_vs_population)" = "COVID-19 (Hospitalized)",
  "COVID-19 Severity (Hospitalisation)" = "COVID-19 (Hospitalized)",
  "Chronic_obstructive_pulmonary_disease" = "Chronic Obstructive Pulmonary",
  "Incident COPD" = "Chronic Obstructive Pulmonary",
  "Hypothyroidism_or_rheumatoid_arthritis_(pleiotropy)" = "Hypothyroidism or Rheumatoid Arthritis",
  "Platelet-to-lymphocyte_ratio" = "Platelet-to-Lymphocyte Ratio",
  "Cancer_(pleiotropy)" = "Cancer",
  "Chronic_inflammatory_diseases_(ankylosing_spondylitis,_Crohn's_disease,_psoriasis,_primary_sclerosing_cholangitis,_ulcerative_colitis)_(pleiotropy)" = "Chronic Inflammatory Diseases",
  "Hay_fever_and/or_eczema" = "Hay Fever / Eczema",
  "Chronic_lymphocytic_leukemia" = "Chronic Lymphocytic Leukemia",
  "Atopic_asthma" = "Atopic Asthma", "Skin_cancer" = "Skin Cancer",
  "Childhood_ear_infection" = "Childhood Ear Infection",
  "Allergic_disease_(asthma,_hay_fever_or_eczema)" = "Allergic Disease (Asthma/Hay Fever/Eczema)",
  "Colorectal_cancer_or_advanced_adenoma" = "Colorectal Cancer / Advanced Adenoma",
  "HIV infection" = "HIV Infection", "Rheumatoid arthritis" = "Rheumatoid Arthritis",
  "Primary Sjogrens syndrome" = "Primary Sjogren's Syndrome",
  "Graves' disease" = "Graves' Disease",
  "Inflammatory bowel disease" = "Inflammatory Bowel Disease",
  "Ulcerative colitis" = "Ulcerative Colitis", "Breast cancer" = "Breast Cancer",
  "Incident Lung Cancer" = "Lung Cancer", "Smoking" = "Smoking"
)

clean_display_name <- function(x) {
  nm <- ifelse(x %in% names(trait_display_map), trait_display_map[x],
               tools::toTitleCase(gsub("_", " ", gsub("\\s*\\(.*\\)$", "", x))))
  nm <- gsub("\\s*\\(?[Pp]leiotropy\\)?", "", nm)
  nm <- gsub("\\s*\\(?MTAG\\)?", "", nm)
  trimws(nm)
}

# Categories are assigned by keyword rather than by an exact-name dictionary, so
# that "Prevalent", "Incident" and self-report variants of the same trait land
# in the same category as the base trait. First match wins.
get_category <- function(tr) {

  tr_lower <- tolower(tr)

  if (grepl("protein levels", tr_lower)) return("Cytokine / Protein Level")
  if (grepl("monocyte count|lymphocyte count|neutrophil count|white blood cell count|eosinophil count|basophil count|platelet-to-lymphocyte", tr_lower)) return("Immune Cell Count")
  if (grepl("c-reactive protein|c reactive protein|\\bcrp\\b", tr_lower)) return("Inflammation Marker")
  if (grepl("cancer|carcinoma|leukemia|leukaemia|lymphoma|adenoma|melanoma|malignan|neoplas", tr_lower)) return("Cancer")
  if (grepl("covid|\\bhiv\\b|infection|influenza|hepatitis|tuberculosis|sepsis", tr_lower)) return("Infectious Disease")
  if (grepl("copd|chronic obstructive pulmonary", tr_lower)) return("Chronic Inflammatory Disease")
  if (grepl("arthritis|lupus|sclerosis|psoriasis|crohn|colitis|inflammatory bowel|dermatitis|sjogren|hypothyroidism|graves|thyroiditis|vasculitis|ankylosing spondylitis", tr_lower)) return("Autoimmune / Inflammatory Disease")
  if (grepl("asthma|hay fever|eczema|allerg|atopic", tr_lower)) return("Immune-Related Disorder")
  if (grepl("^smoking|smoking initiation", tr_lower)) return("Environmental Exposure")

  "Other Immune-Related"
}
get_category <- Vectorize(get_category)


# --- summary of enriched traits ----------------------------------------------

summarize_enrichment <- function(results_data, label) {

  results_data <- copy(results_data)
  results_data[, category := get_category(trait)]

  eligible_immune <- results_data[
    observed >= MIN_OBSERVED_FOR_RESULT &
      expected_mean_null >= MIN_EXPECTED_FOR_RESULT &
      category != "Other Immune-Related"
  ]

  cat("\n============================================================\n")
  cat(label, "ENRICHMENT SUMMARY\n")
  cat("============================================================\n")
  cat("Traits tested (FDR denominator):", nrow(results_data), "\n")
  cat("Traits eligible (observed>=3, expected>=3, immune-related):",
      nrow(eligible_immune), "\n")
  cat("  nominal p<0.05 and OR>1:", nrow(eligible_immune[fisher_p < 0.05 & fisher_OR > 1]), "\n")
  cat("  FDR<0.05, any direction:", nrow(eligible_immune[fisher_p_FDR < 0.05]), "\n")
  cat("  FDR<0.05 and OR>1:",       nrow(eligible_immune[fisher_p_FDR < 0.05 & fisher_OR > 1]), "\n")

  invisible(eligible_immune)
}

summarize_enrichment(ewas_results, "EWAS")
summarize_enrichment(gwas_results, "GWAS")


# --- figure panels -----------------------------------------------------------

categories_actually_used <- character(0)

# The EWAS reference lists every CpG in the catalog associated with each trait,
# not only the cCpGs, so it has to be restricted to the cCpGs before counting.
# The GWAS reference is already restricted by construction.
make_counts_panel <- function(ref_data, cpg_col, trait_col, exclude_pattern_traits, title_letter) {

  dt <- copy(ref_data)
  setnames(dt, cpg_col, "cpg_id")
  setnames(dt, trait_col, "trait")

  dt <- dt[cpg_id %in% ccpgs_clean]

  dist <- dt[, .(n_cpg_hits = uniqueN(cpg_id)), by = trait]

  dist <- dist[!grepl("^Red_blood_cell", trait) &
                 !grepl("[Rr]eticulocyte", trait) &
                 !grepl("\\(UKB_data_field_[0-9]+\\)", trait) &
                 !grepl("pleiotropy", trait, ignore.case = TRUE) &
                 !grepl("osteoarthritis", trait, ignore.case = TRUE) &
                 !grepl("^smoking", trait, ignore.case = TRUE) &
                 !trait %in% exclude_pattern_traits]

  dist[, trait_display := clean_display_name(trait)]
  dist <- dist[order(-n_cpg_hits)][!duplicated(trait_display)]

  d <- head(dist, 20)
  d[, pct_cpgs := 100 * n_cpg_hits / n_ccpg]
  d <- d[order(-n_cpg_hits)]
  d[, trait_display := factor(trait_display, levels = rev(trait_display))]
  d[, category := get_category(trait)]

  categories_actually_used <<- unique(c(categories_actually_used, unique(d$category)))

  ggplot(d, aes(x = pct_cpgs, y = trait_display, fill = category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(n_cpg_hits, " CpGs")), hjust = -0.1, size = 2.8, color = "black") +
    scale_fill_manual(values = GLOBAL_CATEGORY_COLORS, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = "% of Significant CpGs", y = NULL, title = title_letter) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 6, color = "black"),
      axis.title.x       = element_text(face = "bold", size = 9),
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      legend.position    = "none"
    )
}

make_enrichment_panel <- function(results_data, exclude_pattern_traits, title_letter,
                                  n_top = 10, require_or_gt1 = FALSE) {

  d <- results_data[observed >= MIN_OBSERVED_FOR_RESULT &
                      expected_mean_null >= MIN_EXPECTED_FOR_RESULT]
  if (require_or_gt1) d <- d[fisher_OR > 1]

  d <- d[!grepl("^Red_blood_cell", trait) &
           !grepl("[Rr]eticulocyte", trait) &
           !grepl("\\(UKB_data_field_[0-9]+\\)", trait) &
           !grepl("osteoarthritis", trait, ignore.case = TRUE) &
           !grepl("^smoking", trait, ignore.case = TRUE) &
           !trait %in% exclude_pattern_traits]

  d[, trait_display := clean_display_name(trait)]
  d <- d[order(fisher_p)][!duplicated(trait_display)]
  d <- head(d, n_top)
  d[, neg_log10_p := -log10(fisher_p)]
  d[, sig_marker := ifelse(fisher_p_FDR < 0.05, "*", "")]

  # Selected by p-value above, displayed by odds ratio
  d <- d[order(-fisher_OR)]
  d[, trait_display := factor(trait_display, levels = rev(trait_display))]

  ggplot(d, aes(x = fisher_OR, y = trait_display, fill = neg_log10_p)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sig_marker), hjust = -0.3, size = 5, color = "black", fontface = "bold") +
    scale_fill_gradient(low = "#D9EAD3", high = "#1B4D1B",
                        name = expression(-log[10]*"(Fisher's "*italic(p)*")")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Odds Ratio (Enrichment)", y = NULL, title = title_letter) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 6, color = "black"),
      axis.title.x       = element_text(face = "bold", size = 9),
      plot.title         = element_text(face = "bold", size = 11, hjust = 0),
      legend.position    = "right",
      legend.title       = element_text(size = 7),
      legend.text        = element_text(size = 6.5),
      legend.key.size    = unit(0.3, "cm")
    )
}

panel_EWAS_counts <- make_counts_panel(
  ewas_ref, cpg_col = "cpg", trait_col = "trait_canonical",
  exclude_pattern_traits = c("C-reactive protein (CRP) levels", "Smoking"),
  title_letter = "a)")

panel_EWAS_enrichment <- make_enrichment_panel(
  ewas_results, exclude_pattern_traits = c("Smoking"), title_letter = "b)")

panel_GWAS_enrichment <- make_enrichment_panel(
  gwas_results, exclude_pattern_traits = c("Smoking"), title_letter = "c)",
  n_top = 15, require_or_gt1 = TRUE)

panel_GWAS_counts <- make_counts_panel(
  gwas_sig_ref, cpg_col = "cpg", trait_col = "trait",
  exclude_pattern_traits = c("Mean_platelet_thrombocyte_volume_(UKB_data_field_30100)",
                             "Platelet-to-lymphocyte_ratio", "C-reactive_protein_levels",
                             "Smoking"),
  title_letter = "d)")

ggsave(file.path(out_dir, "EWAS_panel_counts.png"),     panel_EWAS_counts,     width = 105, height = 99, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "EWAS_panel_enrichment.png"), panel_EWAS_enrichment, width = 105, height = 99, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "GWAS_panel_enrichment.png"), panel_GWAS_enrichment, width = 105, height = 99, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "GWAS_panel_counts.png"),     panel_GWAS_counts,     width = 105, height = 99, units = "mm", dpi = 300)


# --- shared category legend --------------------------------------------------
# Built as its own plot rather than extracted from one of the panels, so both
# count panels can share a single legend placed independently in the figure.

all_categories <- intersect(names(GLOBAL_CATEGORY_COLORS), categories_actually_used)
n_cat <- length(all_categories)
cat("\nCategories used across both count panels:", paste(all_categories, collapse = ", "), "\n")

legend_df <- data.table(category = all_categories, row = (n_cat - 1):0)

legend_plot <- ggplot(legend_df) +
  geom_tile(aes(x = 0, y = row, fill = category), width = 0.8, height = 0.7, show.legend = FALSE) +
  geom_text(aes(x = 0.55, y = row, label = category), hjust = 0, size = 3.5, color = "black") +
  scale_fill_manual(values = GLOBAL_CATEGORY_COLORS) +
  scale_x_continuous(limits = c(-0.5, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, n_cat), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out_dir, "GLOBAL_category_legend.png"), legend_plot, width = 120, height = 90, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "GLOBAL_category_legend.pdf"), legend_plot, width = 120, height = 90, units = "mm")

cat("\nWritten to ", out_dir, "\n")
