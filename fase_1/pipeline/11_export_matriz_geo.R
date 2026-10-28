# 11 — Processed data export for the GEO deposit
#
# Writes the harmonized beta matrices (one file per source series), the phenotype table and a README
# for the third-party reanalysis deposit GSE349814 (six cohorts, P2 probe mask of the primary analysis).
# Paper section: Data and Code Availability
# Input:  fase_1/data_harmonized_P2/betas_harmonized.rds, fase_1/data_combat_P2/phenotype_final.csv,
#         fase_1/data_pheno/GSE147221_phenodata_crudo.csv
# Output: fase_1/geo_deposito_P2_6coh/ (betas_procesadas_<GSE>.csv.gz, fenotipo_6cohortes.csv, README.txt)
# DECIMALES: decimal places kept for each beta value (default 6)

# ---- 1. Packages, paths and data -----------------------------------------------------------------
suppressMessages({ library(data.table); library(readr) })         # fast reading and writing

# Paths hang from the repository root, found by walking up until the fase_1 folder appears, so the
# script runs from fase_1/pipeline or from the root alike.
raiz <- normalizePath(getwd())                                    # starting directory
while (!dir.exists(file.path(raiz, "fase_1")) &&                  # fase_1 not visible yet
       dirname(raiz) != raiz) {                                   # and not at the filesystem root
  raiz <- dirname(raiz)                                           # go up one level
}
stopifnot(dir.exists(file.path(raiz, "fase_1")))                  # stop if the root was not found

DIR_HARM   <- file.path(raiz, "fase_1/data_harmonized_P2")        # harmonized betas, six cohorts
DIR_COMBAT <- file.path(raiz, "fase_1/data_combat_P2")            # final phenotype table (script 05)
DIR_OUT    <- file.path(raiz, "fase_1/geo_deposito_P2_6coh")      # deposit folder
DECIMALES  <- 6                                                   # decimal places kept per beta
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
cat("Insumo:", DIR_HARM, "-> salida:", DIR_OUT, "\n")

harm  <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))     # list(betas, cohorte) from script 04
betas <- harm$betas; coh <- harm$cohorte; rm(harm); gc()          # probes x samples, source series per column
pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id; pheno <- pheno[colnames(betas), ]   # phenotype in matrix column order
stopifnot(identical(rownames(pheno), colnames(betas)))        # samples aligned one to one
cat("Matriz:", nrow(betas), "sondas ×", ncol(betas), "muestras | cohortes:",
    paste(sort(unique(coh)), collapse = ", "), "\n")
print(table(cohorte = pheno$gse, grupo = pheno$group))            # cases and controls per series

# ---- 2. Beta matrices, one file per source series ------------------------------------------------
# A third-party reanalysis in GEO declares one "Reanalysis of: GSExxxxxx" relation per source series,
# so each file maps to one of those relations (same layout as GSE207605). Probes in rows, first column
# ID_REF, one column per sample named by its GSM. Betas are rounded to six decimals, well below the
# technical resolution of the array.
series <- sort(unique(coh))                                            # the six source series
resumen <- data.frame()                                                # per-file summary
for (g in series) {                                                    # one file per series
  cols <- which(coh == g)                                              # columns of that series
  M    <- round(betas[, cols, drop = FALSE], DECIMALES)                # round to DECIMALES places
  X    <- data.frame(ID_REF = rownames(betas), M, check.names = FALSE) # ID_REF = GEO convention
  ruta <- file.path(DIR_OUT, paste0("betas_procesadas_", g, ".csv.gz"))# file named after the series
  fwrite(X, ruta, compress = "gzip")                                   # gzip-compressed CSV
  mb <- file.info(ruta)$size / 1048576                                 # size in MB for the summary
  cat(sprintf("  %-10s %5d muestras  %8.1f MB\n", g, length(cols), mb))
  resumen <- rbind(resumen, data.frame(serie = g, muestras = length(cols), MB = round(mb, 1)))
  rm(M, X); gc()                                                       # free memory before the next series
}
cat("Total:", sum(resumen$muestras), "muestras en", nrow(resumen), "archivos,",
    round(sum(resumen$MB)), "MB\n")

# ---- 3. Phenotype table and README ---------------------------------------------------------------
# Sample type from source_name_ch1. In GSE147221 cases are whole blood (EDTA tube) and controls are
# buffy coat, so diagnosis and sample type are fully confounded; the other series are whole blood.
cru <- read_csv(file.path(raiz, "fase_1/data_pheno/GSE147221_phenodata_crudo.csv"), show_col_types = FALSE)
# GSE147221 samples labelled buffy coat versus whole blood, keyed by GSM
st  <- setNames(ifelse(grepl("buffy", cru$source_name_ch1, ignore.case = TRUE),
                       "buffy coat", "whole blood (EDTA tube)"), cru$gsm_id)
pheno$sample_type <- ifelse(pheno$gse == "GSE147221", unname(st[pheno$gsm_id]), "whole blood")
stopifnot(!anyNA(pheno$sample_type))                              # every sample has a sample type
print(table(pheno$gse, pheno$sample_type))                        # sample type by series

out_pheno <- file.path(DIR_OUT, "fenotipo_6cohortes.csv")
write_csv(pheno, out_pheno)                                       # one row per sample, key gsm_id
cat("Guardado:", basename(out_pheno), "—", nrow(pheno), "muestras ×", ncol(pheno), "variables\n")

# README of the deposit, in English for GEO curators and data users; matrix and table dimensions are
# read from the data so the text stays in step with the files.
readme <- c(
  "PROCESSED DATA - Blood DNA methylation in schizophrenia, six cohorts reprocessed uniformly",
  "==========================================================================================",
  "",
  "Third-party reanalysis of six publicly available Illumina HumanMethylation450 (GPL13534)",
  "blood series. No new samples were generated. The six series were re-processed from",
  "the least-processed data each of them deposited (methylated and unmethylated signal",
  "intensities; none deposited IDAT files) under a single quality-control, normalisation and",
  "harmonisation pipeline. This deposit contains the resulting processed beta-value matrices",
  "for all six cohorts, and the phenotype table.",
  "",
  "Pre-registration: OSF Registries 2u56r (https://osf.io/2u56r).",
  "",
  "CONTENTS",
  paste0("  betas_procesadas_<GSE>.csv.gz   one file per source series (", nrow(resumen),
         " files). Processed beta values,"),
  paste0("                                  ", format(nrow(betas), big.mark = ","), " probes x ",
         format(ncol(betas), big.mark = ","), " samples in total, rounded to ", DECIMALES, " decimals."),
  "                                  First column ID_REF = 450K probe identifier; one column per",
  "                                  sample, named with its GSM accession in the source series.",
  paste0("  fenotipo_6cohortes.csv          ", format(nrow(pheno), big.mark = ","),
         " rows (one per sample, key gsm_id) x ", ncol(pheno), " variables:"),
  "                                  source series, diagnosis (SCZ/CNT), sex, chronological age",
  "                                  where deposited, estimated leukocyte fractions, epigenetic",
  "                                  smoking class and score, epigenetic ages (Horvath, Hannum),",
  "                                  ancestry principal components, imputed age, and the 80/20",
  "                                  analytical partition (column 'set'), and the sample type",
  "                                  recorded by the source series (column 'sample_type').",
  "",
  "  Files are split by source series because this deposit is a third-party reanalysis: it",
  "  contributes no samples of its own and instead declares one \"Reanalysis of: GSExxxxxx\"",
  "  relation per source series. One file per series makes it explicit which data came from",
  "  where.",
  "",
  "PROCESSING APPLIED (identical for every cohort, starting from methylated/unmethylated signal)",
  "  1. Quality control: detection p-value threshold corrected for multiplicity (Bonferroni over",
  "     the analysed probes); probe- and sample-level call-rate filters; sex verification against",
  "     methylation at the sex chromosomes.",
  "  2. Probe masking: cross-reactive probes (Chen et al. 2013; Benton et al. 2015), the general",
  "     mask of Zhou et al. 2017, probes with a SNP at the CpG or at the single-base-extension",
  "     site (MAF >= 0.05), and the sex chromosomes.",
  "  3. Normalisation: dasen (joint between-array quantile normalisation, wateRmelon) followed",
  "     by BMIQ (within-array type I/II probe correction).",
  "  4. Cross-cohort harmonisation: regRCPqn (regional quantile normalisation by genomic region",
  "     and probe type, on M-values; GSE84727 as reference).",
  "  5. Covariates estimated from the methylation data: leukocyte composition (EpiDISH, RPC),",
  "     smoking (EpiSmokEr), epigenetic age (methylclock: Horvath, Hannum), ancestry principal",
  "     components (EPISTRUCTURE).",
  "",
  "SOURCE SERIES (GEO; Illumina HumanMethylation450 BeadChip, GPL13534)",
  paste0("  ", paste(sort(unique(coh)), collapse = ", ")),
  "",
  "NOTE ON COHORT COMPOSITION",
  "  All six reanalysed series are included. Users should be aware of a design feature of",
  "  GSE147221: according to the sample annotations deposited with that series, all cases were",
  "  sampled as whole blood in EDTA tubes and all controls as buffy coat, so diagnosis is fully",
  "  confounded with sample type in that cohort. In our analyses a subset of case-control",
  "  differences was present in GSE147221 only and absent from the other five cohorts. The",
  "  column 'sample_type' of the phenotype table records this information for every sample.",
  "  An earlier version of this deposit was built with the case/control labels of GSE80417 and",
  "  GSE84727 reversed in our phenotype-extraction code; that error is corrected here.",
  "",
  "PARTITION VARIABLE",
  "  Column 'set' gives the 80/20 partition used in the analysis (train / test). All inference",
  "  in the associated manuscript was restricted to the training partition; the test partition",
  "  is reserved for a subsequent classification study.",
  "",
  "CONTACT",
  "  Hernan Guillermo Hernandez-Hincapie, Universidad Santo Tomas, Bucaramanga, Colombia."
)
writeLines(readme, file.path(DIR_OUT, "README.txt"))             # plain-text README of the deposit
cat("Guardado: README.txt\n")

cat("\n=== ENTREGABLE PARA DEPÓSITO (GEO) ===\n")
# final listing of the deposit folder with file sizes
inf <- file.info(list.files(DIR_OUT, full.names = TRUE))
print(data.frame(archivo = basename(rownames(inf)), MB = round(inf$size / 1e6, 1)), row.names = FALSE)
