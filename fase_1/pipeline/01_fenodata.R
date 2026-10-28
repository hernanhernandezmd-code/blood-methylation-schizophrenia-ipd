# 01 — Phenotype tables
#
# Builds one phenotype table per cohort and two combined tables from the GEO series matrices. Every
# diagnostic mapping is checked against the case and control counts that GEO declares for the series,
# and the script stops if they differ.
# Paper section: Design, Registration and Data Sources
# Input:  GEO series matrices (downloaded with GEOquery into ../data_input/geo_series_matrix)
# Output: <DIR_OUT>/<GSE>_phenodata.csv, phenodata_combinado_sin_na.csv (group, sex and age complete),
#         phenodata_combinado_imputar_edad.csv (group and sex complete; age may be missing)
# DIR_OUT: environment variable PHENO_OUT, default ../data_pheno

suppressMessages({ library(GEOquery); library(Biobase); library(data.table) })
# Allow running from the repository root as well as from fase_1/pipeline
if (!dir.exists("../data_input") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

DIR_GEO <- "../data_input/geo_series_matrix"           # cache of downloaded series matrices
DIR_OUT <- Sys.getenv("PHENO_OUT", "../data_pheno")    # output folder, overridable from the environment
dir.create(DIR_GEO, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

GSE_IDS <- c("GSE80417", "GSE84727", "GSE116378", "GSE116379", "GSE147221", "GSE152027")
SENTRIX_EN_TITLE <- c("GSE147221", "GSE152027")     # chip code is the first token of the title

# Schizophrenia-spectrum diagnoses recorded in GSE147221 (consensus_diagnosis, lower case)
DIAGNOSTICOS_SCZ <- c("schizophrenia", "schizoaffective disorder", "schizophreniform disorder", "schizoid pd",
                      "schizoaffective disorder (galway sz)", "schizophrenia (galway cp)",
                      "schizophrenia (galway prcl sz)")

# Cases and controls declared in each GEO series ("overall design"). The mapping must reproduce them.
# GSE116379 declares more famine cases than it deposited; the deposited counts are used.
CONTEO_GEO <- list(GSE80417 = c(SCZ = 353, CNT = 322), GSE84727 = c(SCZ = 414, CNT = 433),
                   GSE116378 = c(SCZ = 15, CNT = 49), GSE116379 = c(SCZ = 74, CNT = 79),
                   GSE147221 = c(SCZ = 364, CNT = 349), GSE152027 = c(SCZ = 290, CNT = 203))

# Two impossible ages deposited in GSE80417 (a spurious trailing digit); the true value is kept.
CORRECCIONES_EDAD <- c(GSM2126406 = "90", GSM2126750 = "83")

VALORES_AUSENTES <- c("na", "n/a", "nan", "")       # strings treated as missing (compared in lower case)

# Splits "<chip>_R<row>C<col>" into chip, position and chip_position
extraer_sentrix <- function(x) {
  x <- trimws(as.character(x)); i <- regexpr("_R", x, fixed = TRUE)     # locate the "_R" separator
  ok <- !is.na(x) & i > 0                                                 # only well-formed codes are split
  sentrix <- ifelse(ok, substr(x, 1, i - 1), NA); position <- ifelse(ok, substr(x, i + 1, nchar(x)), NA)
  data.frame(sentrix = sentrix, position = position,
             sentrix_position = ifelse(ok, paste0(sentrix, "_", position), NA), stringsAsFactors = FALSE)
}

# ---- 1. One table per cohort --------------------------------------------------------------------------

resultados <- list()                                   # one data frame per cohort, combined in section 2
for (gse_id in GSE_IDS) {
  cat("\n", strrep("=", 60), "\nProcesando:", gse_id, "\n")
  # Sample annotation only; the platform annotation (GPL) is not needed
  pd <- pData(getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = DIR_GEO)[[1]])
  cat("Muestras cargadas:", nrow(pd), "\n")

  # Array methylation controls carry no phenotype
  if (gse_id %in% SENTRIX_EN_TITLE) {
    meth <- grepl("Meth Control", pd$title, ignore.case = TRUE)          # control wells flagged in the title
    if (any(meth)) { cat("Excluidas (Meth Control en title):", sum(meth), "\n"); pd <- pd[!meth, ] }
  }

  # Clinical characteristics: GEOquery stores each "key: value" pair as a "<key>:ch1" column
  ch1 <- grep(":ch1$", names(pd), value = TRUE)
  t <- data.frame(gsm_id = pd$geo_accession, stringsAsFactors = FALSE)   # one row per GSM accession
  for (k in ch1) t[[sub(":ch1$", "", k)]] <- as.character(pd[[k]])      # drop the ":ch1" suffix

  # Chip code and array position, wherever each series stored them
  if (gse_id %in% SENTRIX_EN_TITLE) {
    # Title holds the code, followed by ":" or a space
    ti <- trimws(as.character(pd$title))
    cod <- ifelse(grepl(":", ti, fixed = TRUE), sub(":.*$", "", ti), sub(" .*$", "", ti))
    t <- cbind(t[, setdiff(names(t), c("sentrix", "position")), drop = FALSE], extraer_sentrix(trimws(cod)))
  } else if (all(c("sentrix", "position") %in% names(t))) {
    t$sentrix_position <- paste0(trimws(t$sentrix), "_", trimws(t$position))   # chip and position apart
  } else if ("sentrixids" %in% names(t)) {
    t <- cbind(t, extraer_sentrix(t$sentrixids))                               # combined code field
  } else if ("description" %in% names(pd)) {
    t <- cbind(t, extraer_sentrix(pd$description))                             # code in the description
  }

  # Diagnostic group. In GSE80417 and GSE84727, code 1 = control and 2 = case
  # (checked against the case and control counts declared by GEO; see CONTEO_GEO).
  g <- rep("EXCLUIR", nrow(t))                          # anything not matched below is dropped
  for (i in seq_len(nrow(t))) {
    if (gse_id %in% c("GSE80417", "GSE84727")) {
      # PLINK-style coding; the column name differs between the two series
      v <- trimws(if ("disease status" %in% names(t)) t[i, "disease status"] else t[i, "disease_status"])
      g[i] <- if (identical(v, "1")) "CNT" else if (identical(v, "2")) "SCZ" else "EXCLUIR"
    } else if (gse_id == "GSE116378") {
      # SCZ / CTR labels
      v <- trimws(t$group[i]); g[i] <- if (identical(v, "SCZ")) "SCZ" else if (identical(v, "CTR")) "CNT" else "EXCLUIR"
    } else if (gse_id == "GSE116379") {
      # Labels carry a suffix, so only the prefix is matched
      v <- trimws(t$group[i]); g[i] <- if (is.na(v)) "EXCLUIR" else if (startsWith(v, "SCZ")) "SCZ" else if (startsWith(v, "CTR")) "CNT" else "EXCLUIR"
    } else if (gse_id == "GSE147221") {
      # Controls by status; cases only with a schizophrenia-spectrum consensus diagnosis
      st <- tolower(trimws(t$status[i])); dx <- tolower(trimws(t$consensus_diagnosis[i]))
      g[i] <- if (identical(st, "control")) "CNT" else if (!is.na(dx) && dx %in% DIAGNOSTICOS_SCZ) "SCZ" else "EXCLUIR"
    } else if (gse_id == "GSE152027") {
      # SCZ / CON labels
      v <- trimws(t$status[i]); g[i] <- if (identical(v, "SCZ")) "SCZ" else if (identical(v, "CON")) "CNT" else "EXCLUIR"
    }
  }
  t$group <- g
  # Drop other diagnoses and unlabelled samples
  if (any(g == "EXCLUIR")) { cat("Excluidas por group no clasificable:", sum(g == "EXCLUIR"), "\n"); t <- t[g != "EXCLUIR", ] }
  print(table(t$group))

  # Stop if the mapping does not reproduce the counts declared by GEO
  for (k in c("SCZ", "CNT")) {
    if (sum(t$group == k) != CONTEO_GEO[[gse_id]][[k]])
      stop(gse_id, ": ", k, " = ", sum(t$group == k), ", GEO declara ", CONTEO_GEO[[gse_id]][[k]])
  }
  cat("Contraste con GEO: SCZ =", sum(t$group == "SCZ"), "| CNT =", sum(t$group == "CNT"), "OK\n")

  # Age and sex under the names each series used
  col_edad <- intersect(c("age", "Age", "age (years)", "ageatbloodcollection"), names(t))[1]   # first match
  col_sexo <- intersect(c("Sex", "sex", "gender", "Gender"), names(t))[1]
  # Harmonized per-cohort table; age and sex stay as deposited text (cleaned downstream)
  out <- data.frame(gsm_id = t$gsm_id, sentrix_position = t$sentrix_position, gse = gse_id,
                    sentrix = t$sentrix, position = t$position,
                    edad = if (!is.na(col_edad)) t[[col_edad]] else NA_character_,
                    sexo = if (!is.na(col_sexo)) t[[col_sexo]] else NA_character_,
                    group = t$group, stringsAsFactors = FALSE)
  # Apply the two age corrections, with a console trace of each
  for (gsm in names(CORRECCIONES_EDAD)) {
    j <- which(out$gsm_id == gsm)
    if (length(j)) { cat("Corrección de edad:", gsm, out$edad[j], "->", CORRECCIONES_EDAD[[gsm]], "\n"); out$edad[j] <- CORRECCIONES_EDAD[[gsm]] }
  }
  fwrite(out, file.path(DIR_OUT, paste0(gse_id, "_phenodata.csv")), na = "")   # missing written as empty
  resultados[[gse_id]] <- out
}

# ---- 2. Combined tables -------------------------------------------------------------------------------

todos <- rbindlist(resultados)                         # all eligible samples of the six cohorts
valido <- function(x) !is.na(x) & !(tolower(trimws(x)) %in% VALORES_AUSENTES)   # non-missing value

estricta <- todos[valido(todos$sexo) & valido(todos$edad), ]      # group, sex and age complete
imputar  <- todos[valido(todos$sexo), ]                            # age may be missing (imputed later)
fwrite(estricta, file.path(DIR_OUT, "phenodata_combinado_sin_na.csv"), na = "")
fwrite(imputar,  file.path(DIR_OUT, "phenodata_combinado_imputar_edad.csv"), na = "")   # used by 02

# Report any diagnostic group that a filter removed entirely
for (nm in c("estricta", "imputar")) {
  # Cohort x group counts before (a) and after (b) the filter, on the same levels
  d <- get(nm); a <- table(todos$gse, todos$group); b <- table(factor(d$gse, rownames(a)), factor(d$group, colnames(a)))
  cat("\n[", nm, "] total =", nrow(d), "\n"); print(b)
  vacios <- which(a > 0 & b == 0, arr.ind = TRUE)      # cells present before and empty after
  if (nrow(vacios)) cat("AVISO: grupos eliminados por completo:", paste(rownames(a)[vacios[, 1]], colnames(a)[vacios[, 2]]), "\n")
}
# Samples without age in the inclusive table; their age is imputed in 05
cat("\nSin edad en la base inclusiva:", sum(!valido(imputar$edad)), "\n")
