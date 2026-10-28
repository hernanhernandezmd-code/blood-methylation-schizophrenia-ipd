# 02 — Download of raw signal intensities and quality control
#
# Downloads the methylated/unmethylated signal (and detection p-values where deposited) of the six
# Illumina 450K cohorts from GEO, maps columns to GSM accessions, applies detection call-rate filters,
# masks problematic probes, checks declared sex and writes one curated M/U object per cohort. None of
# the cohorts deposited IDAT files, so QC starts from the least processed signal available.
# Paper section: Quality Control and Probe Filtering
# Input:  ../data_pheno/phenodata_combinado_imputar_edad.csv; GEO supplementary signal files;
#         ../refs/zhou2017_mask_general_450k.txt, ../refs/naeem2014_problematic_450k.txt
# Output: ../data_curated/<GSE>_curated[_P2|_P3].rds (list of M and U, probes x samples, no NA)
# VARIANTE_MASKING: "P1" | "P2" | "P3" (default "P2"); probe-masking stringency, see section 5b
# Labels in console output: E1/E2 = probe/sample call-rate filters; E8 = probe mask.
# Memory: the six cohorts are held in RAM between sections (about 35 GB at the peak).

# ---- 1. Utility functions and packages ----------------------------------------------------------------

# Installs missing Bioconductor/CRAN packages and loads them all
InsLoad.pks <- function(pqks, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))  # BiocManager installs from both repositories
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  ip <- rownames(installed.packages())                 # packages already installed
  new.pkg <- pqks[!(pqks %in% ip)]                     # packages still missing
  if (length(new.pkg) > 0)
    BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)   # handles CRAN too
  out <- sapply(pqks, require, character.only = TRUE, quietly = TRUE)          # TRUE = loaded
  print(out)                                            # load status of each package
  if (length(new.pkg) > 0)
    message("Instalados ", length(new.pkg), " paquetes nuevos.")
  else
    message("Todos los paquetes requeridos ya estaban instalados.")
  invisible(out)
}

# Wraps each named function so that it resets the seed before every call
sealseed <- function(func_names, seed = 7777777) {
  for (func_name in func_names) {
    original_func <- get(func_name, envir = .GlobalEnv)                 # function to wrap
    new_func <- function(...) {                         # wrapper with the same arguments
      set.seed(as.integer(get("seed", envir = environment())))          # same seed on every call
      get("original_func", envir = environment())(...)                  # then the original call
    }
    # Private environment holding the original function and the seed
    environment(new_func) <- new.env(parent = environment(original_func))
    assign("original_func", original_func, envir = environment(new_func))
    assign("seed", seed, envir = environment(new_func))
    assign(func_name, new_func, envir = .GlobalEnv)                     # replace it globally
  }
}

# Bioconductor packages used across the pipeline (installed here once)
pqks_BC <- c(
  "GEOquery", "minfi", "minfiData", "illuminaio", "wateRmelon",
  "ENmix", "sva", "limma", "BiocParallel",
  "DMRcate", "missMethyl", "bacon",
  "FlowSorted.Blood.450k", "FlowSorted.Blood.EPIC"
)
# CRAN packages
pqks_CRAN <- c(
  "data.table", "dplyr", "readr", "stringr",
  "tibble", "ggplot2", "matrixStats", "jsonlite",
  "yaml", "httr2", "purrr", "glue",
  "gtsummary", "officer", "caret", "renv", "remotes"
)

InsLoad.pks(c(pqks_BC, pqks_CRAN))                      # install what is missing, then load everything

# regRCPqn is distributed only on GitHub
if (!requireNamespace("regRCPqn", quietly = TRUE))
  remotes::install_github("regRCPqn/regRCPqn")
library(regRCPqn)                                       # between-cohort harmonization, used in 04

# ---- 2. Paths and parameters --------------------------------------------------------------------------
set.seed(7777777)                                       # fixed seed for the session

# Download and QC are deterministic; sealseed is applied in the downstream scripts.

# Run from fase_1/pipeline/ whether started from the repository root or from that folder
if (!file.exists("../data_pheno") && dir.exists("fase_1/pipeline")) {
  setwd("fase_1/pipeline")
}

# The six cohorts, in processing order
GSE_IDS <- c("GSE84727", "GSE80417", "GSE147221",
             "GSE116378", "GSE116379", "GSE152027")

DIR_INPUT   <- "../data_input"     # downloaded signal files and series matrices
DIR_CURATED <- "../data_curated"   # curated M/U per cohort (input of 03)
DIR_PHENO   <- "../data_pheno"     # phenotype tables written by 01

dir.create(DIR_INPUT,   showWarnings = FALSE, recursive = TRUE)  # created only if absent
dir.create(DIR_CURATED, showWarnings = FALSE, recursive = TRUE)  # created only if absent

options(timeout = 3600)   # large GEO files exceed the 60 s default

# Call-rate filters. Failed measurements are not set to NA; only whole probes or samples are
# dropped, as in wateRmelon pfilter. The probe threshold (3%) is stricter than the 5% of
# Lehne 2015; the sample threshold stays at 10% to keep samples that will need imputed age.
PVAL_MAX_FRAC_MU  <- 0.03   # max. fraction of failed samples per probe (E1)
PVAL_MAX_FRAC_SM  <- 0.10   # max. fraction of failed probes per sample (E2)

# Detection uses the detection p-values deposited by the authors, with a Bonferroni threshold
# over the analysed probes (section 7b), not 0.01 and not 1e-16 (Heiss & Just 2019).
# GSE80417 has no detection p-values and gets a same-laboratory intensity threshold (7b).

# Laboratory that produced each detP (Mill = methylumi/Exeter; Boks = GenomeStudio/Utrecht);
# used only to label diagnostics
LAB <- c(GSE84727 = "Mill", GSE80417 = "Mill", GSE147221 = "Mill",
         GSE152027 = "Mill", GSE116378 = "Boks", GSE116379 = "Boks")

# ---- 3. Harmonized phenotype data ---------------------------------------------------------------------
# Inclusive table (diagnosis and sex required, age optional): samples without age are downloaded
# and quality-controlled so that their age can be imputed from epigenetic clocks downstream.
pheno <- readr::read_csv(
  file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
  show_col_types = FALSE
)

pheno$edad_imputada <- is.na(pheno$edad)   # TRUE = chronological age missing, to be imputed

# Sample counts, diagnostic groups and missing age per cohort
cat("Número de muestras en phenodata:", nrow(pheno), "\n")
cat("Con edad real:", sum(!pheno$edad_imputada),
    "| sin edad (a imputar):", sum(pheno$edad_imputada), "\n")
cat("\nTabla de grupos por cohortes:\n")
print(table(pheno$gse, pheno$group))                    # cases and controls per cohort
cat("\nEdad faltante por cohorte (TRUE = a imputar):\n")
print(table(pheno$gse, edad_imputada = pheno$edad_imputada))  # samples to impute, per cohort

# ---- 4. Selective download of signal files ------------------------------------------------------------
# getGEOSuppFiles() is not used: the RAW.tar of these series holds no IDAT files (only the 450K
# manifest) and the beta files are large. Only raw methylated/unmethylated signal and, where
# available, detection p-values are fetched from the NCBI FTP, with a size check against the
# server Content-Length. Layouts differ: separate M/U files (GSE84727, GSE80417) or one combined
# file with Detection/Methylated/Unmethylated columns (the other four). GSE80417 has no detP.
SUPP_FILES <- list(
  GSE84727  = c("GSE84727_methylatedIntensities.csv.gz",
                "GSE84727_unmethylatedIntensities.csv.gz",
                "GSE84727_detectionP.csv.gz"),
  GSE80417  = c("GSE80417_methylatedIntensities.csv.gz",
                "GSE80417_unmethylatedIntensities.csv.gz"),   # no detP deposited
  GSE147221 = c("GSE147221_Dublin_blood_raw_Signal.csv.gz"),   # raw signal, not the processed file
  GSE116378 = c("GSE116378_DUTCHSCZ_Matrix_signal_intensities.csv.gz"),
  GSE116379 = c("GSE116379_CHINFAM_Matrix_signal_intensities.csv.gz"),
  GSE152027 = c("GSE152027_IOP_raw_Signal.csv.gz")
)

options(timeout = 9200)   # the largest files are over 2 GB

for (gse in GSE_IDS) {
  dir_gse <- file.path(DIR_INPUT, gse)                  # one folder per cohort
  dir.create(dir_gse, showWarnings = FALSE, recursive = TRUE)
  message("\n== ", gse, " ==")

  pref <- sub("[0-9]{3}$", "nnn", gse)   # FTP folder: GSE84727 -> GSE84nnn

  for (archivo in SUPP_FILES[[gse]]) {                  # only the files listed above
    url     <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
                       pref, gse, archivo)
    destino <- file.path(dir_gse, archivo)              # local copy

    # Expected size from the Content-Length header
    h        <- tryCatch(curlGetHeaders(url, redirect = TRUE),
                         error = function(e) character(0))           # no header if the request fails
    cl       <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)  # size header, if any
    esperado <- if (length(cl)) as.numeric(trimws(sub(".*:", "", tail(cl, 1)))) else NA_real_   # bytes

    # Skip only files that are already complete
    if (file.exists(destino) && !is.na(esperado) && file.size(destino) == esperado) {
      message("  ya completo: ", archivo, " (", round(esperado / 1048576), " MB)")
      next
    }

    # Download (or re-download a partial file) in binary mode
    message("  descargando: ", archivo,
            if (!is.na(esperado)) paste0(" (", round(esperado / 1048576), " MB)") else "")
    download.file(url, destino, mode = "wb", quiet = FALSE, method = "libcurl")

    # A size mismatch means a truncated download; warn instead of stopping
    if (!is.na(esperado) && file.size(destino) != esperado)
      warning(archivo, ": tamaño descargado (", file.size(destino),
              ") != esperado (", esperado, "). Revisar.")
  }
}

# List downloaded files
for (gse in GSE_IDS) {
  archivos <- list.files(file.path(DIR_INPUT, gse), full.names = FALSE)  # files present on disk
  cat(gse, ":", length(archivos), "archivo(s)\n")
  cat("  ", paste(archivos, collapse = "\n   "), "\n\n")
}

# ---- 4b. Series matrix for GSE116378 and GSE116379 ----------------------------------------------------
# These cohorts label signal columns with sample names (DUTCHSCZ-01, CHINFAM-01) instead of
# Sentrix IDs; the series matrix provides the Sample_title -> GSM accession table.
dir_sm <- file.path(DIR_INPUT, "geo_series_matrix")     # same cache folder used by 01
dir.create(dir_sm, showWarnings = FALSE, recursive = TRUE)

for (gse in c("GSE116378", "GSE116379")) {
  pref    <- sub("[0-9]{3}$", "nnn", gse)               # FTP folder, as in section 4
  archivo <- paste0(gse, "_series_matrix.txt.gz")       # GEO series matrix file
  destino <- file.path(dir_sm, archivo)                 # local copy
  url     <- sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s",
    pref, gse, archivo
  )

  # Same completeness check as in section 4
  h        <- tryCatch(curlGetHeaders(url, redirect = TRUE),
                       error = function(e) character(0))
  cl       <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)  # size header, if any
  esperado <- if (length(cl)) as.numeric(trimws(sub(".*:", "", tail(cl, 1)))) else NA_real_

  if (file.exists(destino) && !is.na(esperado) && file.size(destino) == esperado) {
    message(gse, " series_matrix: ya completo (", round(esperado / 1024), " KB)")
    next
  }

  message(gse, " series_matrix: descargando...")
  download.file(url, destino, mode = "wb", quiet = TRUE, method = "libcurl")
  message("  listo: ", round(file.size(destino) / 1024), " KB")
}

# ---- 5. QC parameters and 450K annotation -------------------------------------------------------------
# Signal file layout per cohort
FORMATO <- c(GSE84727 = "separado", GSE80417 = "separado",
             GSE147221 = "combinado", GSE152027 = "combinado",
             GSE116378 = "combinado", GSE116379 = "combinado")

# ---- 5a. SNP and sex-chromosome probes ----------------------------------------------------------------
# Illumina 450K annotation (hg19): chromosome and SNP columns for every probe
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       update = FALSE, ask = FALSE)
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))

anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)   # 485,512 probes

# Two SNP definitions; the variant chosen in 5b decides which one applies
MAF_SNP <- 0.05   # minimum minor allele frequency for a "common" SNP
# P1: any SNP in the probe body
sondas_snp_probe  <- rownames(anno)[!is.na(anno$Probe_rs)]
# P2/P3: SNP at the target CpG or single-base extension site with MAF >= MAF_SNP
snp_cpg <- !is.na(anno$CpG_rs) & !is.na(anno$CpG_maf) & anno$CpG_maf >= MAF_SNP   # at the CpG
snp_sbe <- !is.na(anno$SBE_rs) & !is.na(anno$SBE_maf) & anno$SBE_maf >= MAF_SNP   # at the extension base
sondas_snp_cpgsbe <- rownames(anno)[snp_cpg | snp_sbe]

# Sex-chromosome probes are removed from the analysis (sex-dependent methylation)
sondas_x    <- rownames(anno)[anno$chr == "chrX"]       # X-chromosome probes
sondas_y    <- rownames(anno)[anno$chr == "chrY"]       # Y-chromosome probes
sondas_sexo <- c(sondas_x, sondas_y)   # X and Y are also used separately for the sex check (10)

cat("Filtro SNP — amplio (Probe_rs):", length(sondas_snp_probe),
    "| protocolo (CpG/SBE, MAF>=", MAF_SNP, "):", length(sondas_snp_cpgsbe), "\n")
cat("Sondas en cromosoma X:", length(sondas_x),
    "| Y:", length(sondas_y),
    "| total sexuales:", length(sondas_sexo), "\n")

# ---- 5b. Problematic probes: masking variant P1/P2/P3 -------------------------------------------------
# Sources: cross-reactive probes (Chen 2013 + Benton 2015, via maxprobes), design/mapping mask
# (Zhou 2017, MASK_general) and "discard" probes (Naeem 2014). Zhou and Naeem lists are read from
# cached files in ../refs/ (see the *.SOURCE.txt files there).
if (!requireNamespace("maxprobes", quietly = TRUE))
  remotes::install_github("markgene/maxprobes")
sondas_chen_benton <- maxprobes::xreactive_probes(array_type = "450K")   # cross-reactive probes

# Zhou 2017 general mask (mapping, design and SNP problems); empty with a warning if missing
zhou_file   <- file.path("..", "refs", "zhou2017_mask_general_450k.txt")
sondas_zhou <- if (file.exists(zhou_file)) readLines(zhou_file) else character(0)
if (length(sondas_zhou) == 0)
  warning("Lista Zhou 2017 no encontrada en ", zhou_file,
          " — la unión no la incluye. Ver refs/*.SOURCE.txt para regenerarla.")

# Naeem 2014 "discard" list (used only by P1); empty with a warning if missing
naeem_file   <- file.path("..", "refs", "naeem2014_problematic_450k.txt")
sondas_naeem <- if (file.exists(naeem_file)) readLines(naeem_file) else character(0)
if (length(sondas_naeem) == 0)
  warning("Lista Naeem 2014 no encontrada en ", naeem_file,
          " — la unión no la incluye. Ver refs/*.SOURCE.txt para regenerarla.")

# Masking stringency:
#   "P1" = Chen+Benton + Zhou + Naeem, SNP anywhere in the probe (Probe_rs)
#   "P2" = Chen+Benton + Zhou, SNP at CpG/SBE with MAF >= 0.05 (the registered specification)
#   "P3" = Chen+Benton only, SNP at CpG/SBE
VARIANTE_MASKING <- "P2"                                # P2 = variant used in the paper

# Union of problematic-probe lists for each variant, and the SNP set that goes with it
xreact_por_variante <- list(
  P1 = unique(c(sondas_chen_benton, sondas_zhou, sondas_naeem)),
  P2 = unique(c(sondas_chen_benton, sondas_zhou)),
  P3 = unique(sondas_chen_benton))
sondas_xreact <- xreact_por_variante[[VARIANTE_MASKING]]  # problematic probes of the chosen variant
sondas_snp    <- if (VARIANTE_MASKING == "P1") sondas_snp_probe else sondas_snp_cpgsbe

# Size of each source and of the lists actually applied
cat("VARIANTE DE MASKING:", VARIANTE_MASKING, "\n")
cat("Sondas problemáticas por fuente (450K):\n")
cat("  Chen 2013 + Benton 2015 (maxprobes):", length(sondas_chen_benton), "\n")
cat("  Zhou 2017 (MASK_general):           ", length(sondas_zhou),
    if (VARIANTE_MASKING == "P3") "(NO usada en P3)" else "", "\n")
cat("  Naeem 2014 (discard):               ", length(sondas_naeem),
    if (VARIANTE_MASKING != "P1") "(NO usada fuera de P1)" else "", "\n")
cat("  cross-react/mapeo aplicadas:        ", length(sondas_xreact), "\n")
cat("  SNP aplicadas:                      ", length(sondas_snp),
    if (VARIANTE_MASKING == "P1") "(Probe_rs amplio)" else "(CpG/SBE protocolo)", "\n")

# Master exclusion list (SNP + sex chromosomes + problematic). Sex-chromosome probes are copied
# aside before removal, for the sex check in section 10.
sondas_quitar <- unique(c(sondas_snp, sondas_sexo, sondas_xreact))
cat("Total sondas en lista de exclusión:    ", length(sondas_quitar),
    "(de", nrow(anno), "→ quedarían", nrow(anno) - length(sondas_quitar), ")\n")

# ---- 6. Read M / U / detP matrices per cohort ---------------------------------------------------------
# Output: datos_brutos[[gse]] with M, U and P (probes x samples, raw column IDs; P is NULL when
# no detP exists). In combined files "Unmethylated" contains "Methylated", hence the setdiff.
datos_brutos <- list()                                  # raw matrices per cohort

for (gse in GSE_IDS) {
  message("\n######## ", gse, " (", FORMATO[gse], ") ########")
  dir_gse  <- file.path(DIR_INPUT, gse)
  archivos <- file.path(dir_gse, SUPP_FILES[[gse]])     # local paths of this cohort's files

  if (FORMATO[gse] == "combinado") {                    # combined layout
    # Locate the signal file of this cohort
    f  <- grep("Signal|signal|intensit", archivos, ignore.case = TRUE, value = TRUE)[1]
    message("  [6] leyendo (combinado): ", basename(f))
    dt <- data.table::fread(f, data.table = FALSE)       # fread: fast on multi-GB files
    rownames(dt) <- dt[[1]]; dt[[1]] <- NULL             # first column = probe ID
    cols     <- colnames(dt)                            # three columns per sample, told apart by suffix
    col_pval <- grep("Detection|Pval",   cols, ignore.case = TRUE, value = TRUE)   # detection p-values
    col_unm  <- grep("Unmethylated",     cols, ignore.case = TRUE, value = TRUE)   # unmethylated signal
    col_met  <- setdiff(grep("Methylated", cols, ignore.case = TRUE, value = TRUE), col_unm)   # methylated

    M <- as.matrix(dt[, col_met,  drop = FALSE])        # methylated signal
    U <- as.matrix(dt[, col_unm,  drop = FALSE])        # unmethylated signal
    P <- as.matrix(dt[, col_pval, drop = FALSE])        # detection p-values
    # Strip the suffix so the three matrices share the same sample IDs
    colnames(M) <- trimws(sub("[ _]+methylated[ _]*signal$",   "", col_met,  ignore.case = TRUE))
    colnames(U) <- trimws(sub("[ _]+unmethylated[ _]*signal$", "", col_unm,  ignore.case = TRUE))
    colnames(P) <- trimws(sub("[ _]+detection[ _]*pval$",      "", col_pval, ignore.case = TRUE))
    rm(dt)                                              # free the data frame

  } else {  # separate files
    # Methylated file name also matches "unmethylated", so exclude it explicitly
    fm <- archivos[grepl("methylatedInt",  archivos, ignore.case = TRUE) &
                   !grepl("unmethyl",      archivos, ignore.case = TRUE)]
    fu <- archivos[grepl("unmethylatedInt", archivos, ignore.case = TRUE)]  # unmethylated file
    fp <- archivos[grepl("detectionP",     archivos, ignore.case = TRUE)]   # empty for GSE80417
    message("  [6] leyendo (separado): ", basename(fm), " + ", basename(fu))
    dm <- data.table::fread(fm, data.table = FALSE); rownames(dm) <- dm[[1]]; dm[[1]] <- NULL
    M  <- as.matrix(dm); rm(dm)                         # methylated signal
    du <- data.table::fread(fu, data.table = FALSE); rownames(du) <- du[[1]]; du[[1]] <- NULL
    U  <- as.matrix(du); rm(du)                         # unmethylated signal
    if (length(fp) == 1) {                              # detP file deposited
      dp <- data.table::fread(fp, data.table = FALSE); rownames(dp) <- dp[[1]]; dp[[1]] <- NULL
      P  <- as.matrix(dp); rm(dp)                       # detection p-values
    } else {
      P <- NULL                                         # detection falls back to intensity (7b)
      message("  [6] sin detP (", gse, "): detección solo por intensidad.")
    }
  }

  # Same sample order in M, U and P
  ids <- intersect(colnames(M), colnames(U))
  if (!is.null(P)) ids <- intersect(ids, colnames(P))   # and in P, when deposited
  M <- M[, ids, drop = FALSE]                           # common column order
  U <- U[, ids, drop = FALSE]
  if (!is.null(P)) P <- P[, ids, drop = FALSE]
  message("  [6] sondas leídas: ", nrow(M), " | columnas (IDs brutos): ", length(ids))

  datos_brutos[[gse]] <- list(M = M, U = U, P = P)      # columns still carry raw IDs
  rm(M, U, P, ids); gc()                                # release the copies before the next cohort
}

# ---- 7. Column IDs to GSM accessions ------------------------------------------------------------------
# Mill/Exeter cohorts use Sentrix IDs (matched to pheno$sentrix_position); Boks/Utrecht cohorts use
# sample names, matched through the series matrix from 4b. Initial counts are kept for the summary.
datos_gsm <- list()                                     # matrices with GSM column names

for (gse in GSE_IDS) {
  M   <- datos_brutos[[gse]]$M                          # methylated signal
  U   <- datos_brutos[[gse]]$U                          # unmethylated signal
  P   <- datos_brutos[[gse]]$P                          # detection p-values (NULL for GSE80417)
  ids <- colnames(M)                                    # raw column IDs
  datos_brutos[[gse]] <- NULL   # free memory for this cohort
  gc()

  # Sentrix IDs if most columns match the phenotype table; sample names otherwise
  if (mean(ids %in% pheno$sentrix_position) > 0.5) {
    gsm <- pheno$gsm_id[match(ids, pheno$sentrix_position)]   # NA for samples not in the phenotype
    message(gse, " [7] crosswalk: sentrix → gsm_id")
  } else {
    sm_file <- file.path(DIR_INPUT, "geo_series_matrix",
                         paste0(gse, "_series_matrix.txt.gz"))
    sm      <- readLines(gzfile(sm_file))               # whole series matrix (small file)
    # Parallel rows of sample titles and GSM accessions in the series matrix header
    tit     <- gsub('"', '', strsplit(grep("^!Sample_title", sm, value = TRUE)[1], "\t")[[1]][-1])
    geo     <- gsub('"', '', strsplit(grep("^!Sample_geo_accession", sm, value = TRUE)[1], "\t")[[1]][-1])
    tit_pf  <- trimws(sub(":.*$", "", tit))   # "DUTCHSCZ-01: description" -> "DUTCHSCZ-01"
    gsm     <- geo[match(ids, tit_pf)]
    message(gse, " [7] crosswalk: nombre → gsm_id (vía series_matrix)")
    rm(sm, tit, geo, tit_pf)
  }

  # Keep mapped samples only (unmapped columns are ineligible: other diagnoses or array controls)
  ok <- !is.na(gsm)
  M  <- M[, ok, drop = FALSE]; U <- U[, ok, drop = FALSE]  # same columns in M and U
  if (!is.null(P)) P <- P[, ok, drop = FALSE]
  colnames(M) <- gsm[ok]; colnames(U) <- gsm[ok]        # columns renamed to GSM accessions
  if (!is.null(P)) colnames(P) <- gsm[ok]
  message(gse, " [7] muestras mapeadas a gsm_id: ", sum(ok),
          " (sin match: ", sum(!ok), ")")

  # Starting counts of probes and samples, reported in the funnels of section 12
  datos_gsm[[gse]] <- list(
    M             = M,
    U             = U,
    P             = P,
    n_muestras_ini = ncol(M),                           # samples mapped to GSM
    n_sondas_ini   = nrow(M)                            # probes in the deposited file
  )
  rm(M, U, P); gc()
}
rm(datos_brutos); gc()                                  # raw list no longer needed

# ---- 7b. Detection threshold: per-probe Bonferroni ----------------------------------------------------
# DETP_STAR = 0.05 / number of probes analysed (FWER 5%). The family is the probes that survive the
# list-based masking, not the whole array, because masked probes are discarded regardless.
n_sondas_analizadas <- nrow(anno) - length(intersect(sondas_quitar, rownames(anno)))
DETP_STAR <- 0.05 / n_sondas_analizadas                 # per-probe Bonferroni threshold
message("[7b] sondas del array: ", nrow(anno),
        " | excluidas por masking: ", length(intersect(sondas_quitar, rownames(anno))),
        " | analizadas: ", n_sondas_analizadas)
message("[7b] DETP_STAR (Bonferroni 0.05/", n_sondas_analizadas, ") = ",
        format(DETP_STAR, scientific = TRUE, digits = 3))

# Diagnostic only: percentage of measurements above DETP_STAR per cohort with detP
for (gse in GSE_IDS) {
  P <- datos_gsm[[gse]]$P
  if (is.null(P)) {
    message(gse, " [7b] (", LAB[[gse]], ") sin detP → usa T_STAR_80417 (abajo).")
    next
  }
  message(gse, " [7b] (", LAB[[gse]], ") % mediciones enmascaradas a DETP_STAR: ",
          round(100 * mean(P > DETP_STAR, na.rm = TRUE), 3), "%")
  rm(P)
}
gc()

# GSE80417 (no detP): T_STAR_80417 is the quantile of total intensity T = M + U that reproduces
# the mean failure rate at DETP_STAR of the three Mill cohorts with detP.
mill_detp <- c("GSE84727", "GSE147221", "GSE152027")   # same laboratory as GSE80417
ff <- c()                                              # failure percentage of each Mill cohort
for (g in mill_detp) {
  Pg <- datos_gsm[[g]]$P                                # detP of one Mill cohort
  ff <- c(ff, 100 * mean(Pg > DETP_STAR, na.rm = TRUE))  # % of measurements failing
}
target_f <- mean(ff)                                   # target failure rate for GSE80417 (%)

T80          <- datos_gsm[["GSE80417"]]$M + datos_gsm[["GSE80417"]]$U   # total intensity
T_STAR_80417 <- as.numeric(quantile(T80, probs = target_f / 100, na.rm = TRUE))
rm(T80); gc()

message("[7b] GSE80417: objetivo de fallo = ", round(target_f, 3),
        "% (media Mill en DETP_STAR) → T_STAR_80417 = ", round(T_STAR_80417))

# ---- 8. Detection calls and call-rate filters (filter-only, no NA) ------------------------------------
# A measurement fails if its detP > DETP_STAR (GSE80417: T = M + U < T_STAR_80417). Failed cells
# keep their raw value; the failure matrix only feeds the per-probe (> 3%) and per-sample (> 10%)
# exclusion fractions. Output: datos_det[[gse]] with filtered M and U.
datos_det     <- list()                                 # call-rate-filtered matrices per cohort
diag_frac     <- list()                      # probes that would be dropped at each candidate threshold
UMBRALES_DIAG <- c(0.01, 0.03, 0.05, 0.10)   # candidate probe thresholds, diagnostic only

for (gse in GSE_IDS) {
  M <- datos_gsm[[gse]]$M
  U <- datos_gsm[[gse]]$U
  P <- datos_gsm[[gse]]$P
  datos_gsm[[gse]]$M <- NULL   # free raw matrices; keep metadata
  datos_gsm[[gse]]$U <- NULL
  datos_gsm[[gse]]$P <- NULL
  gc()

  Tt <- M + U                                           # total intensity per measurement
  if (!is.null(P)) {                                    # cohorts with deposited detP
    fallo <- P > DETP_STAR                              # TRUE = failed detection
    # Diagnostic: agreement between the detP mask and an intensity threshold calibrated within the
    # same cohort to mask the same fraction; supports the intensity fallback used for GSE80417.
    frac_self <- mean(fallo, na.rm = TRUE)              # fraction failing by detP
    T_self    <- as.numeric(quantile(Tt, probs = frac_self, na.rm = TRUE))  # cutoff failing the same fraction
    conc      <- mean((Tt < T_self) == fallo, na.rm = TRUE)  # cell-by-cell agreement
    message(gse, " [8] regla: detP > ", DETP_STAR,
            " | concordancia intensidad auto-calibrada: ",
            round(100 * conc, 1), "%")
  } else {
    fallo <- Tt < T_STAR_80417   # GSE80417: same-laboratory intensity threshold
    message(gse, " [8] regla: T < ", round(T_STAR_80417),
            " (sin detP; intensidad calibrada mismo-lab).")
  }
  pct_fallidas <- round(100 * mean(fallo, na.rm = TRUE), 2)   # failed measurements, whole cohort
  message(gse, " [8] % mediciones fallidas: ", pct_fallidas, "%")

  frac_por_muestra <- colMeans(fallo, na.rm = TRUE)     # failed fraction of probes, per sample
  frac_por_sonda   <- rowMeans(fallo, na.rm = TRUE)     # failed fraction of samples, per probe

  # Probes exceeding each candidate threshold (reported, not applied)
  n_sobre <- integer(length(UMBRALES_DIAG))
  for (k in seq_along(UMBRALES_DIAG)) {
    n_sobre[k] <- sum(frac_por_sonda > UMBRALES_DIAG[k])
  }
  names(n_sobre) <- paste0(">", 100 * UMBRALES_DIAG, "%")  # labels ">1%" to ">10%"
  diag_frac[[gse]] <- n_sobre
  message(gse, " [8] sondas con fallo en >1/3/5/10% de muestras: ",
          paste(n_sobre, collapse = " / "),
          "  (umbral aplicado E1: ", round(100 * PVAL_MAX_FRAC_MU), "%)")

  muestras_malas   <- colnames(M)[frac_por_muestra > PVAL_MAX_FRAC_SM]   # E2: > 10% failed probes
  sondas_malas     <- rownames(M)[frac_por_sonda   > PVAL_MAX_FRAC_MU]   # E1: > 3% failed samples

  # Drop failing probes and samples; U follows the dimensions of M
  M <- M[!rownames(M) %in% sondas_malas, !colnames(M) %in% muestras_malas, drop = FALSE]
  U <- U[rownames(M), colnames(M), drop = FALSE]        # keep U aligned with M
  message(gse, " [8] filtro E1/E2: -", length(muestras_malas), " muestras, -",
          length(sondas_malas), " sondas")

  # Filtered matrices plus what each step removed, for the summary
  datos_det[[gse]] <- list(
    M                = M,
    U                = U,
    n_muestras_ini   = datos_gsm[[gse]]$n_muestras_ini,
    n_sondas_ini     = datos_gsm[[gse]]$n_sondas_ini,
    muestras_excl_det = muestras_malas,                 # samples removed by E2
    sondas_excl_det  = sondas_malas,                    # probes removed by E1
    pct_fallidas     = pct_fallidas                     # failed measurements (%)
  )
  rm(M, U, P, Tt, fallo); gc()
}
rm(datos_gsm); gc()                                     # previous stage no longer needed

# Probes dropped per cohort at each candidate threshold (">3%" is the applied one)
cat("\n[8] Sondas descartadas por umbral candidato (filas = cohorte):\n")
print(do.call(rbind, diag_frac))

# ---- 9. Problematic-probe filter ----------------------------------------------------------------------
# Sex-chromosome probes are copied to M_sex/U_sex for section 10, then sondas_quitar is removed.
datos_filt <- list()                                    # mask-filtered matrices per cohort

for (gse in GSE_IDS) {
  M <- datos_det[[gse]]$M
  U <- datos_det[[gse]]$U
  datos_det[[gse]]$M <- NULL                            # free memory; counts stay
  datos_det[[gse]]$U <- NULL
  gc()

  # Sex-chromosome signal kept aside for the sex check
  sx    <- intersect(sondas_sexo, rownames(M))          # X/Y probes present in this cohort
  M_sex <- M[sx, , drop = FALSE]
  U_sex <- U[sx, , drop = FALSE]

  # Remove every probe of the master exclusion list present in this cohort
  quitar <- intersect(sondas_quitar, rownames(M))
  M <- M[!rownames(M) %in% quitar, , drop = FALSE]
  U <- U[rownames(M), , drop = FALSE]                   # keep U aligned with M

  # Probes removed by each component (they overlap)
  n_snp_excl    <- length(intersect(sondas_snp,    quitar))
  n_sexo_excl   <- length(intersect(sondas_sexo,   quitar))
  n_xreact_excl <- length(intersect(sondas_xreact, quitar))
  message(gse, " [9] sondas excluidas — SNP: ", n_snp_excl,
          " | sexo: ", n_sexo_excl,
          " | problemáticas (unión): ", n_xreact_excl,
          " | restantes: ", nrow(M))

  # Carry forward matrices, sex-probe signal and all counts so far
  datos_filt[[gse]] <- list(
    M                = M,
    U                = U,
    M_sex            = M_sex,
    U_sex            = U_sex,
    n_muestras_ini   = datos_det[[gse]]$n_muestras_ini,
    n_sondas_ini     = datos_det[[gse]]$n_sondas_ini,
    muestras_excl_det = datos_det[[gse]]$muestras_excl_det,
    sondas_excl_det  = datos_det[[gse]]$sondas_excl_det,
    pct_fallidas     = datos_det[[gse]]$pct_fallidas,
    n_mascara_excl   = length(quitar),    # probes removed by the union of all masks
    n_snp_excl       = n_snp_excl,        # mask components overlap and do not add up
    n_sexo_excl      = n_sexo_excl,
    n_xreact_excl    = n_xreact_excl
  )
  rm(M, U, M_sex, U_sex); gc()
}
rm(datos_det); gc()                                     # previous stage no longer needed

# ---- 10. Declared versus predicted sex ----------------------------------------------------------------
# Sex is predicted from the contrast of median log2 total intensity on Y versus X probes, split into
# two k-means clusters (fixed seed); the cluster with the higher Y - X is male. Clustering instead of
# a median cut keeps the prediction independent of the cohort sex ratio. Discordant samples are
# removed.
for (gse in GSE_IDS) {
  M     <- datos_filt[[gse]]$M
  U     <- datos_filt[[gse]]$U
  M_sex <- datos_filt[[gse]]$M_sex
  U_sex <- datos_filt[[gse]]$U_sex
  datos_filt[[gse]]$M_sex <- NULL                       # sex-probe signal is not needed after this
  datos_filt[[gse]]$U_sex <- NULL
  gc()

  pheno_gse <- pheno[match(colnames(M), pheno$gsm_id), ]   # phenotype rows in sample order
  discord   <- character(0)                                # stays empty if sex cannot be checked

  p_x <- rownames(M_sex) %in% sondas_x                  # X-chromosome rows
  p_y <- rownames(M_sex) %in% sondas_y                  # Y-chromosome rows

  # Only possible when both X and Y probes survived (the Boks/Utrecht files carry none)
  if (sum(p_x) > 0 && sum(p_y) > 0) {
    Tsex <- log2(M_sex + U_sex + 1)                     # log2 total intensity; +1 avoids log(0)
    medX <- matrixStats::colMedians(Tsex[p_x, , drop = FALSE], na.rm = TRUE)  # median X per sample
    medY <- matrixStats::colMedians(Tsex[p_y, , drop = FALSE], na.rm = TRUE)  # median Y per sample
    yx   <- medY - medX                                 # high in males, low in females
    names(yx) <- colnames(Tsex)

    ok_yx <- is.finite(yx)                              # samples with a usable contrast
    pred  <- setNames(rep(NA_character_, length(yx)), names(yx))  # NA until predicted
    if (sum(ok_yx) > 2) {                               # k-means needs more than two values
      set.seed(7777777)                                 # reproducible k-means start
      km      <- kmeans(yx[ok_yx], centers = 2)         # two clusters on Y - X
      grupo_M <- which.max(tapply(yx[ok_yx], km$cluster, mean))  # higher Y - X = male
      pred[ok_yx] <- ifelse(km$cluster == grupo_M, "M", "F")  # label each sample by its cluster
    }
    decl    <- toupper(substr(pheno_gse$sexo, 1, 1))    # declared sex as "M"/"F"
    discord <- colnames(M)[!is.na(pred) & !is.na(decl) & pred != decl]  # predicted and declared sex differ
    message(gse, " [10] sexo predicho M/F: ", sum(pred == "M", na.rm = TRUE),
            "/", sum(pred == "F", na.rm = TRUE),
            " | declarado M/F: ", sum(decl == "M", na.rm = TRUE),
            "/", sum(decl == "F", na.rm = TRUE),
            " | discordancias: ", length(discord))
    if (length(discord) > 0)
      message("      GSMs excluidos: ", paste(discord, collapse = ", "))
  } else {
    message(gse, " [10] sin sondas X e Y suficientes para chequeo de sexo.")   # declared sex is used
  }

  # Remove discordant samples
  M <- M[, !colnames(M) %in% discord, drop = FALSE]
  U <- U[, colnames(M), drop = FALSE]                   # keep U aligned with M

  datos_filt[[gse]]$M       <- M                        # sex-checked matrices
  datos_filt[[gse]]$U       <- U
  datos_filt[[gse]]$discord <- discord                  # kept for the sample funnel
  rm(M, U, M_sex, U_sex); gc()
}

# ---- 11. Save curated matrices and per-cohort QC counts -----------------------------------------------
# Writes <GSE>_curated[_suffix].rds = list(M, U), probes x samples, input of 03_normalizacion.R.
resumen_qc <- list()                                    # per-cohort QC counts
paneles    <- list()   # probe set of each curated cohort, for the intersection in section 12

for (gse in GSE_IDS) {
  M <- datos_filt[[gse]]$M
  U <- datos_filt[[gse]]$U

  # P1 has no file suffix; P2/P3 add "_P2"/"_P3"
  suf_var <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)
  rds_out <- file.path(DIR_CURATED, paste0(gse, "_curated", suf_var, ".rds"))  # e.g. GSE84727_curated_P2.rds
  saveRDS(list(M = M, U = U), rds_out)                  # raw signal, no NA, ready for normalization

  paneles[[gse]] <- rownames(M)                         # curated probe set of this cohort

  # Counts per QC stage, used for the funnels in section 12
  resumen_qc[[gse]] <- data.frame(
    gse                    = gse,
    regla_deteccion        = if (gse == "GSE80417") paste0("T<", round(T_STAR_80417))
                             else paste0("detP>", signif(DETP_STAR, 3)),
    pct_fallidas           = datos_filt[[gse]]$pct_fallidas,  # failed measurements (%)
    # probe funnel
    sondas_in              = datos_filt[[gse]]$n_sondas_ini,  # probes in the deposited file
    sondas_excl_callrate   = length(datos_filt[[gse]]$sondas_excl_det),  # probe call-rate (3%)
    sondas_excl_mascara    = datos_filt[[gse]]$n_mascara_excl,           # probe mask (union)
    sondas_cur             = nrow(M),                   # probes after QC
    # mask components (overlapping)
    n_snp                  = datos_filt[[gse]]$n_snp_excl,  # SNP probes removed
    n_sexo                 = datos_filt[[gse]]$n_sexo_excl,  # X/Y probes removed
    n_xreact               = datos_filt[[gse]]$n_xreact_excl,  # problematic-list probes removed
    # sample funnel
    muestras_in            = datos_filt[[gse]]$n_muestras_ini,  # samples mapped to GSM
    muestras_excl_callrate = length(datos_filt[[gse]]$muestras_excl_det), # sample call-rate (10%)
    muestras_excl_sexo     = length(datos_filt[[gse]]$discord),  # sex-discordant samples
    muestras_cur           = ncol(M),                   # samples after QC
    stringsAsFactors       = FALSE
  )

  message(gse, " [11] guardado: ", rds_out,
          " (", nrow(M), " sondas × ", ncol(M), " muestras)")
  datos_filt[[gse]] <- NULL; rm(M, U); gc()              # release the cohort once saved
}

# ---- 12. Overall QC summary ---------------------------------------------------------------------------
# Probe and sample funnels (input -> exclusions per stage -> curated), mask composition, detection
# rule per cohort, and the probe set common to the six cohorts that enters 03_normalizacion.R.
q <- do.call(rbind, resumen_qc)                         # one row per cohort
rownames(q) <- NULL

# Probe funnel: array -> call-rate -> mask -> curated
ts <- data.frame(
  cohorte     = q$gse,
  entrada     = q$sondas_in,                            # probes in the deposited file
  menos_E1_3pct = q$sondas_excl_callrate,               # removed by probe call-rate
  menos_E8_mask = q$sondas_excl_mascara,                # removed by the mask
  curadas     = q$sondas_cur,                           # probes after QC
  stringsAsFactors = FALSE
)
# Totals row and percentage of probes retained
ts <- rbind(ts, data.frame(cohorte = "TOTAL",
  entrada = sum(ts$entrada), menos_E1_3pct = sum(ts$menos_E1_3pct),
  menos_E8_mask = sum(ts$menos_E8_mask), curadas = sum(ts$curadas)))
ts$pct_retenido <- paste0(round(100 * ts$curadas / ts$entrada, 1), "%")  # share of probes retained
for (col in c("entrada", "menos_E1_3pct", "menos_E8_mask", "curadas"))   # thousands separator
  ts[[col]] <- format(ts[[col]], big.mark = ".", decimal.mark = ",", trim = TRUE)
cat("══════ Embudo de SONDAS ══════\n")
cat("entrada  −menos_E1_3pct (call-rate detP)  −menos_E8_mask (máscara unión)  = curadas\n\n")
print(ts, row.names = FALSE)

# Mask composition; components overlap, so they do not add up to the union removed
tc <- data.frame(
  cohorte        = q$gse,
  SNP            = q$n_snp,                             # SNP filter of the chosen variant
  sexo_XY        = q$n_sexo,                            # sex-chromosome probes
  xreact_Naeem   = q$n_xreact,            # cross-reactive/mapping lists of the chosen variant
  UNION_quitada  = q$sondas_excl_mascara,               # probes actually removed
  stringsAsFactors = FALSE
)
for (col in c("SNP", "sexo_XY", "xreact_Naeem", "UNION_quitada"))   # thousands separator
  tc[[col]] <- format(tc[[col]], big.mark = ".", decimal.mark = ",", trim = TRUE)
cat("\n══════ Composición de la máscara E8 (componentes con solapamiento) ══════\n")
cat("SNP + sexo_XY + xreact_Naeem NO suman la UNION_quitada (se cruzan entre listas).\n\n")
print(tc, row.names = FALSE)

# Sample funnel: input -> call-rate -> sex discordance -> curated
tm <- data.frame(
  cohorte          = q$gse,
  entrada          = q$muestras_in,                     # samples mapped to GSM
  menos_E2_10pct   = q$muestras_excl_callrate,          # removed by sample call-rate
  menos_sexo       = q$muestras_excl_sexo,              # removed for sex discordance
  curadas          = q$muestras_cur,                    # samples after QC
  stringsAsFactors = FALSE
)
# Totals row
tm <- rbind(tm, data.frame(cohorte = "TOTAL",
  entrada = sum(tm$entrada), menos_E2_10pct = sum(tm$menos_E2_10pct),
  menos_sexo = sum(tm$menos_sexo), curadas = sum(tm$curadas)))
cat("\n══════ Embudo de MUESTRAS ══════\n")
cat("entrada  −menos_E2_10pct (call-rate)  −menos_sexo (discordancia)  = curadas\n\n")
print(tm, row.names = FALSE)

# Detection rule and percentage of failed measurements per cohort
td <- data.frame(
  cohorte          = q$gse,
  regla            = q$regla_deteccion,                 # detP or intensity rule
  pct_med_falladas = paste0(q$pct_fallidas, "%"),       # failed measurements (%)
  stringsAsFactors = FALSE
)
cat("\n══════ Detección (regla y % de mediciones falladas, filter-only sin NA) ══════\n\n")
print(td, row.names = FALSE)

# Probe set common to the six cohorts (03 normalizes on this intersection)
sondas_inter <- NULL                                    # running intersection
for (gse in GSE_IDS) {
  sondas_inter <- if (is.null(sondas_inter)) paneles[[gse]]  # first cohort starts the set
                  else intersect(sondas_inter, paneles[[gse]])  # then intersect with each next one
}
mediana_cur <- median(q$sondas_cur)             # median curated probes per cohort
cat("\n══════ PANEL FINAL para el paso 03 ══════\n")
cat("Intersección común a las 6 cohortes : ",           # 397,194 probes under P2
    format(length(sondas_inter), big.mark = ".", decimal.mark = ","), " CpG\n", sep = "")
cat("Retención al cruzar (vs mediana)    : ",                 # how much the intersection costs
    round(100 * length(sondas_inter) / mediana_cur, 1), "%\n", sep = "")
cat("Muestras curadas totales            : ", sum(q$muestras_cur), "\n", sep = "")  # 2,931 of 2,945

resumen_qc <- q   # keep the full per-cohort table
