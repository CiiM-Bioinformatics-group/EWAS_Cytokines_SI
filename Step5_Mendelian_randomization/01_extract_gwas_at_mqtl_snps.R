#!/usr/bin/env Rscript
# Pull the disease GWAS effect estimates at the SNPs used as mQTL instruments.
#
# The GWAS summary statistics come as bgzipped VCFs covering the whole genome,
# but only the rows at mQTL SNP positions are needed. The file is streamed in
# chunks and filtered as it goes, so it never has to be held in memory.
#
# Usage:
#   Rscript 01_extract_gwas_at_mqtl_snps.R <gwas.vcf.gz> <mqtl_file> <trait> <out.txt>

args      <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 01_extract_gwas_at_mqtl_snps.R <gwas.vcf.gz> <mqtl_file> <trait> <out.txt>")
}

vcf_file  <- args[1]
mqtl_file <- args[2]
trait     <- args[3]
out_file  <- args[4]

CHUNK_SIZE <- 500000


# --- positions to keep --------------------------------------------------------

cat("Reading mQTL positions from", mqtl_file, "\n")

mqtl_snp_ids <- read.table(mqtl_file, header = FALSE, sep = " ",
                           stringsAsFactors = FALSE, quote = "\"")[, 1]

pos_keys <- unique(vapply(strsplit(mqtl_snp_ids, "_"),
                          function(x) paste(x[1], x[2], sep = "_"),
                          character(1)))

cat("  Unique mQTL SNP positions:", length(pos_keys), "\n\n")


# --- stream the VCF and write the matching rows -------------------------------

cat("Streaming", vcf_file, "\n")

con <- gzfile(vcf_file, "rt")

repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) stop("Reached end of file without finding the #CHROM header line.")
  if (startsWith(line, "#CHROM")) break
}

out_con <- file(out_file, "w")

total_scanned       <- 0
total_written       <- 0
total_skipped_parse <- 0

repeat {
  lines <- readLines(con, n = CHUNK_SIZE)
  if (length(lines) == 0) break
  total_scanned <- total_scanned + length(lines)

  # Split and match the whole chunk at once rather than line by line
  split_all <- strsplit(lines, "\t")
  chrom <- vapply(split_all, function(x) x[1], character(1))
  pos   <- vapply(split_all, function(x) x[2], character(1))
  key   <- paste(chrom, pos, sep = "_")

  match_idx <- which(key %in% pos_keys)

  for (i in match_idx) {

    parts  <- split_all[[i]]
    ref    <- parts[4]
    alt    <- parts[5]
    fmt    <- parts[9]
    sample <- parts[10]

    # Effect size, standard error and -log10 p, wherever they sit in FORMAT
    fmt_fields <- strsplit(fmt, ":")[[1]]
    es_idx <- match("ES", fmt_fields)
    se_idx <- match("SE", fmt_fields)
    lp_idx <- match("LP", fmt_fields)
    if (any(is.na(c(es_idx, se_idx, lp_idx)))) {
      total_skipped_parse <- total_skipped_parse + 1
      next
    }

    sv     <- strsplit(sample, ":")[[1]]
    es_val <- suppressWarnings(as.numeric(sv[es_idx]))
    se_val <- suppressWarnings(as.numeric(sv[se_idx]))
    lp_val <- suppressWarnings(as.numeric(sv[lp_idx]))
    if (any(is.na(c(es_val, se_val, lp_val))) || se_val == 0) {
      total_skipped_parse <- total_skipped_parse + 1
      next
    }

    # CHROM already carries the chr prefix
    snp_id <- paste0(chrom[i], "_", pos[i], "_", ref, "_", alt)
    pval   <- 10^(-lp_val)
    tstat  <- es_val / se_val

    writeLines(paste(snp_id, trait, es_val, tstat, pval, sep = "\t"), out_con)
    total_written <- total_written + 1
  }

  cat("  Scanned:", total_scanned, " | written so far:", total_written, "\n")
}

close(con)
close(out_con)

cat("\nVCF rows scanned            :", total_scanned, "\n")
cat("Rows at mQTL positions      :", total_written + total_skipped_parse, "\n")
cat("Rows written                :", total_written, "\n")
cat("Rows skipped, unparseable   :", total_skipped_parse, "\n")
cat("Output:", out_file, "\n")
