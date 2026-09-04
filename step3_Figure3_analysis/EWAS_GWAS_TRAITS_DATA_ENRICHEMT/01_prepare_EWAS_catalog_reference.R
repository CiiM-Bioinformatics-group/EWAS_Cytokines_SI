# =============================================================================
# EWAS Catalog reference file for the immune-trait enrichment analysis
#
# Builds the set of previously reported CpG-trait associations that cCpGs are
# tested against, restricted to
#   1. associations at p < 1e-5
#   2. blood-derived tissues, since the methylation data come from PBMCs
#   3. immune-related traits
#
# Traits are kept in tiers so they can be reported together or separately:
#   Tier 1   immune disease and phenotype traits
#   Tier 1b  cancer traits
#   Tier 1c  environmental exposure
#   Tier 2   cytokine and protein level traits, which are mechanistically
#            relevant here because the phenotype is a cytokine response
#
# Tissue filtering is applied to the tissue column alone rather than to the
# whole row, so that trait names such as "blood pressure" are not matched.
#
# Age is deliberately excluded from every tier. Genome-wide age related drift
# affects a large share of the array, so including it would swamp the
# trait-specific signal rather than add to it.
#
# Input
#   ewascatalog-results.txt.gz   CpG level summary statistics
#   ewascatalog-studies.txt.gz   study level metadata including tissue and trait
#
# Output
#   EWAS_catalogue_BLOOD_IMMUNE_p1e5_FINAL.tsv   reference used by the
#                                                enrichment script
#   plus per-tier files and the tissue and trait audit tables
# =============================================================================

library(data.table)

data_dir <- "data/EWAS_catalog"
out_dir  <- "results/EWAS_GWAS_enrichment"

PVAL_THRESHOLD <- 1e-5

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# --- load and merge ----------------------------------------------------------

cat("Loading raw files...\n")
results <- fread(file.path(data_dir, "ewascatalog-results.txt.gz"))
studies <- fread(file.path(data_dir, "ewascatalog-studies.txt.gz"))

# The results file carries the CpG identifier twice
stopifnot(sum(names(results) == "cpg") == 2)
dup_idx <- which(names(results) == "cpg")[2]
names(results)[dup_idx] <- "cpg_dup"
stopifnot(identical(results$cpg, results$cpg_dup))
results[, cpg_dup := NULL]

n_before <- nrow(results)
results <- results[!is.na(p) & p != "" & !is.na(as.numeric(p))]
cat(sprintf("Rows excluded for missing or invalid p-value: %d\n", n_before - nrow(results)))
results[, p := as.numeric(p)]

cat("Merging results with study metadata on study_id...\n")
merged <- merge(results,
                studies[, .(study_id, trait, tissue, n, author, pmid)],
                by = "study_id", all.x = TRUE)

cat("Total merged rows:", nrow(merged), "\n")
cat("Rows with missing tissue after join:", sum(is.na(merged$tissue)), "\n")
cat("Rows with missing trait after join:",  sum(is.na(merged$trait)),  "\n")

fwrite(merged, file.path(out_dir, "EWAS_catalogue_prepared.tsv"), sep = "\t")

tissue_counts <- merged[, .N, by = tissue][order(-N)]
fwrite(tissue_counts, file.path(out_dir, "tissue_value_counts.tsv"), sep = "\t")


# --- blood tissues, p < 1e-5 -------------------------------------------------

blood_pattern <- paste(
  "whole blood", "peripheral blood", "cord blood", "blood spot",
  "dried blood spot", "white blood cell", "leukocyte", "leucocyte",
  "buffy coat", "PBMC", "peripheral blood mononuclear cell",
  "monocyte", "granulocyte", "neutrophil", "lymphocyte", "lymphoblast",
  "\\bT cell", "\\bT-cell", "\\bB cell", "\\bB-cell",
  "natural killer", "\\bNK cell", "^blood$",
  sep = "|"
)

blood_sig <- merged[p < PVAL_THRESHOLD & grepl(blood_pattern, tissue, ignore.case = TRUE)]

cat("Significant (p<1e-5) blood-tissue rows:", nrow(blood_sig), "\n")
cat("Unique CpGs (blood, p<1e-5):", length(unique(blood_sig$cpg)), "\n")

# Any tissue label containing "blood" that the pattern did not match
missed_blood <- tissue_counts[grepl("blood", tissue, ignore.case = TRUE) &
                                !grepl(blood_pattern, tissue, ignore.case = TRUE)]
if (nrow(missed_blood) > 0) {
  warning("Blood-tissue pattern missed some labels containing 'blood' - see missed_blood_tissues.tsv")
  fwrite(missed_blood, file.path(out_dir, "missed_blood_tissues.tsv"), sep = "\t")
} else {
  cat("No blood-tissue labels were missed by the pattern.\n")
}

fwrite(blood_sig[, .(cpg, trait, tissue, p, beta, se, study_id, gene, chr, pos)],
       file.path(out_dir, "EWAS_catalogue_BLOOD_TISSUE_p1e5.tsv"), sep = "\t")


# --- trait vocabulary --------------------------------------------------------

trait_counts <- blood_sig[, .N, by = trait][order(-N)]
fwrite(trait_counts, file.path(out_dir, "trait_counts_raw.tsv"), sep = "\t")

# Drop UniProt accessions used as trait placeholders, and covariate labels
uniprot_pattern <- "^[A-Z][0-9][A-Z0-9]{4,}$"
trait_counts_filtered <- trait_counts[
  !grepl(uniprot_pattern, trait) &
    !(trait %in% c("Sex", "sex", "Smoking", "Tissue"))
]
fwrite(trait_counts_filtered, file.path(out_dir, "trait_counts_filtered.tsv"), sep = "\t")

immune_pattern <- paste(
  "immune", "inflamm", "autoimmune", "arthritis", "lupus", "asthma",
  "allerg", "psoriasis", "eczema", "dermatitis", "crohn", "colitis",
  "IBD", "celiac", "coeliac", "multiple sclerosis", "type 1 diabetes",
  "thyroiditis", "graves", "hashimoto", "infection", "sepsis", "\\bHIV\\b",
  "hepatitis", "tuberculosis", "malaria", "covid", "coronavirus",
  "sars-cov-2", "sars cov", "cytokine", "interleukin", "\\bIL-?[0-9]",
  "TNF", "interferon", "CRP", "C-reactive", "rheumatoid", "vasculitis",
  "sarcoidosis", "Sjogren", "scleroderma", "myositis", "Guillain",
  "Kawasaki", "ankylosing spondylitis", "COPD", "periodontitis",
  "leukemia", "lymphoma", "immunodeficiency", "complement", "antibody",
  "autoantibod",
  sep = "|"
)

immune_traits_all <- trait_counts_filtered[grepl(immune_pattern, trait, ignore.case = TRUE)]
fwrite(immune_traits_all, file.path(out_dir, "immune_related_traits.tsv"), sep = "\t")

# Cancer and exposure traits are named explicitly rather than matched by
# keyword, which would pick up unrelated traits mentioning the word
cancer_traits_allowlist <- c("Breast cancer", "breast cancer",
                             "Incident Lung Cancer", "Incident Ovarian Cancer")
cancer_traits_all <- trait_counts_filtered[trait %in% cancer_traits_allowlist]
fwrite(cancer_traits_all, file.path(out_dir, "cancer_related_traits.tsv"), sep = "\t")
cat("Cancer traits found:", nrow(cancer_traits_all), "\n")

# Smoking modulates circulating inflammatory markers and cell counts, but it is
# an exposure rather than a disease, so it is kept as its own tier
exposure_traits_allowlist <- c("Smoking", "smoking")
exposure_traits_all <- trait_counts_filtered[trait %in% exposure_traits_allowlist]
fwrite(exposure_traits_all, file.path(out_dir, "exposure_related_traits.tsv"), sep = "\t")
cat("Environmental exposure traits found:", nrow(exposure_traits_all), "\n")


# --- tiers -------------------------------------------------------------------

is_protein_level <- grepl("protein levels", immune_traits_all$trait, ignore.case = TRUE)
is_gene_expr     <- grepl("Gene expression of Affymetrix", immune_traits_all$trait, ignore.case = TRUE)

tier1_disease_phenotype <- immune_traits_all[!is_protein_level & !is_gene_expr]
tier2_protein_levels    <- immune_traits_all[is_protein_level]
tier3_gene_expression   <- immune_traits_all[is_gene_expr]

fwrite(tier1_disease_phenotype, file.path(out_dir, "immune_traits_disease_phenotype.tsv"), sep = "\t")
fwrite(tier2_protein_levels,    file.path(out_dir, "immune_traits_protein_levels.tsv"),    sep = "\t")
fwrite(tier3_gene_expression,   file.path(out_dir, "immune_traits_gene_expression.tsv"),   sep = "\t")

cat("Tier 1 (disease/phenotype) traits:", nrow(tier1_disease_phenotype), "\n")
cat("Tier 2 (protein levels) traits:",    nrow(tier2_protein_levels),    "\n")
cat("Tier 3 (gene expression) traits:",   nrow(tier3_gene_expression),   "\n")

# Gene expression traits are an eQTL type association rather than a disease or
# protein level one, and are not carried into the final file


# --- canonical trait names ---------------------------------------------------
# The same underlying trait appears under different study phrasings, and the
# same protein appears under several aptamer identifiers. Collapsing them stops
# one trait being counted several times as if independent.

trait_map <- c(
  "Prevalent Rheumatoid Arthritis (Self-report)" = "Rheumatoid arthritis",
  "Incident Rheumatoid Arthritis"                = "Rheumatoid arthritis",
  "C-reactive protein"                           = "C-reactive protein (CRP) levels",
  "C-Reactive Protein"                           = "C-reactive protein (CRP) levels",
  "Prevalent COPD (Self-report)"                 = "Incident COPD",
  "Human immunodeficiency virus"                 = "HIV infection",
  "Prevalent Osteoarthritis (Self-report)"       = "Incident Osteoarthritis"
)

cancer_trait_map   <- c("breast cancer" = "Breast cancer")
exposure_trait_map <- c("smoking" = "Smoking")

tier2_protein_levels[, trait_canonical_protein := trimws(sub("\\s*\\(SeqId\\s*=.*\\)$", "", trait))]
cat(sprintf("Tier 2: %d SeqId-level labels collapse to %d canonical proteins\n",
            uniqueN(tier2_protein_levels$trait),
            uniqueN(tier2_protein_levels$trait_canonical_protein)))

tier1b_cancer   <- cancer_traits_all
tier1c_exposure <- exposure_traits_all


# --- final reference file ----------------------------------------------------

tier1_traits  <- tier1_disease_phenotype$trait
tier1b_traits <- tier1b_cancer$trait
tier1c_traits <- tier1c_exposure$trait
tier2_traits  <- tier2_protein_levels$trait

final_data <- blood_sig[trait %in% c(tier1_traits, tier1b_traits, tier1c_traits, tier2_traits)]

final_data[, trait_tier := fcase(
  trait %in% tier1_traits,  "Tier1_immune_disease_phenotype",
  trait %in% tier1b_traits, "Tier1b_cancer",
  trait %in% tier1c_traits, "Tier1c_environmental_exposure",
  trait %in% tier2_traits,  "Tier2_cytokine_protein_level"
)]

final_data[, trait_canonical := fcase(
  trait_tier == "Tier1_immune_disease_phenotype",
    ifelse(trait %in% names(trait_map), trait_map[trait], trait),
  trait_tier == "Tier1b_cancer",
    ifelse(trait %in% names(cancer_trait_map), cancer_trait_map[trait], trait),
  trait_tier == "Tier1c_environmental_exposure",
    ifelse(trait %in% names(exposure_trait_map), exposure_trait_map[trait], trait),
  trait_tier == "Tier2_cytokine_protein_level",
    trimws(sub("\\s*\\(SeqId\\s*=.*\\)$", "", trait))
)]

cat("\n=== FINAL REFERENCE FILE ===\n")
cat("Rows:", nrow(final_data), "\n")
cat("Unique CpGs:", length(unique(final_data$cpg)), "\n")
cat("Unique canonical traits:", length(unique(final_data$trait_canonical)), "\n")
print(final_data[, .(n_rows = .N, n_cpgs = uniqueN(cpg)), by = trait_tier])

keep_cols <- c("cpg", "trait", "trait_canonical", "trait_tier", "tissue", "p",
               "beta", "se", "study_id", "gene", "chr", "pos", "n", "author", "pmid")

fwrite(final_data[, ..keep_cols],
       file.path(out_dir, "EWAS_catalogue_BLOOD_IMMUNE_p1e5_FINAL.tsv"), sep = "\t")

# Per-tier files, for analyses that keep the tiers separate
tier_files <- c(
  Tier1_immune_disease_phenotype = "EWAS_catalogue_BLOOD_IMMUNE_TIER1_ONLY_p1e5.tsv",
  Tier1b_cancer                  = "EWAS_catalogue_BLOOD_CANCER_TIER1b_ONLY_p1e5.tsv",
  Tier1c_environmental_exposure  = "EWAS_catalogue_BLOOD_EXPOSURE_TIER1c_ONLY_p1e5.tsv",
  Tier2_cytokine_protein_level   = "EWAS_catalogue_BLOOD_IMMUNE_TIER2_ONLY_p1e5.tsv"
)

for (tier in names(tier_files)) {
  fwrite(final_data[trait_tier == tier, ..keep_cols],
         file.path(out_dir, tier_files[[tier]]), sep = "\t")
}

cat("\nWritten to ", out_dir, "\n")
