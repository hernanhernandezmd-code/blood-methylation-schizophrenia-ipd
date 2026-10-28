# 03 — Normalization (joint between-array quantile normalization + BMIQ)
#
# Assembles the curated methylated/unmethylated signals of all cohorts on their common probes,
# applies dasen (wateRmelon) jointly to every sample, then BMIQ within each array to correct the
# type I/II design bias. One method per axis: dasen aligns samples to each other, BMIQ maps type II
# onto the type I scale. Between-cohort harmonization (regRCPqn) follows in 04.
# Paper section: Normalization and Harmonization
# Input:  data_curated/<GSE>_curated_P2.rds (M/U matrices from 02, no NA)
# Output: data_normalized_P2/betas_normalized.rds = list(betas, cohorte);
#         data_normalized_P2/betas_dasen.rds (dasen checkpoint)
# VARIANTE_MASKING: "P1" | "P2" | "P3" (default "P2"); must match 02.
# EXCLUIR_COHORTE: "" (default, all cohorts) or a GSE ID to drop; adds a suffix to folders.
# CORRER_DIAGNOSTICOS: FALSE (default) | TRUE; runs optional method comparisons (blocks 3b-3e).
# FORCE_DASEN: FALSE (default) | TRUE; recomputes dasen even if the checkpoint exists.

# ---- 1. Helper functions and packages ---------------------------------------

# Installs missing packages (Bioconductor or CRAN) and loads all of them.
InsLoad.pks <- function(pqks, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")   # installer for Bioc
  ip <- rownames(installed.packages())                  # packages already present
  new.pkg <- pqks[!(pqks %in% ip)]                      # those still missing
  if (length(new.pkg) > 0)
    BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)  # Bioc and CRAN alike
  out <- sapply(pqks, require, character.only = TRUE, quietly = TRUE)   # load; TRUE/FALSE per package
  print(out)                                            # shows any package that failed to load
  if (length(new.pkg) > 0)
    message("Instalados ", length(new.pkg), " paquetes nuevos.")
  else
    message("Todos los paquetes requeridos ya estaban instalados.")
  invisible(out)   # returned for inspection
}


# Wraps stochastic functions so that each call sets the same seed (applied to BMIQ in step 5).
# BMIQ draws random subsets to fit its mixture, so each sample would otherwise vary run to run.
sealseed <- function(func_names, seed = 7777777) {
  for (func_name in func_names) {
    original_func <- get(func_name, envir = .GlobalEnv)   # function to wrap
    new_func <- function(...) {
      set.seed(as.integer(get("seed", envir = environment())))   # reseed before every call
      get("original_func", envir = environment())(...)           # then run the original
    }
    environment(new_func) <- new.env(parent = environment(original_func))  # private closure
    assign("original_func", original_func, envir = environment(new_func))  # keep original inside
    assign("seed", seed, envir = environment(new_func))                    # and the seed
    assign(func_name, new_func, envir = .GlobalEnv)       # replace the global name with the wrapper
  }
}


pqks_BC <- c(
  "minfi", "wateRmelon", "ENmix",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",       # probe type I/II for BMIQ
  "BiocParallel"
)
pqks_CRAN <- c(
  "data.table", "matrixStats",
  "readr", "stringr", "ggplot2"
)

InsLoad.pks(c(pqks_BC, pqks_CRAN))                      # only what normalization needs


# ---- 2. Paths and parameters ------------------------------------------------
set.seed(7777777)                                       # project seed for this script

# Moves to the pipeline folder when started from the repository root.
if (!file.exists("../data_curated") && dir.exists("fase_1/pipeline")) {
  setwd("fase_1/pipeline")                              # all paths below are relative to it
}

GSE_IDS <- c("GSE84727", "GSE80417", "GSE147221",   # order of assembly
             "GSE116378", "GSE116379", "GSE152027")     # the six 450K blood cohorts

# Optional cohort exclusion, applied before assembly so it does not enter the joint normalization.
EXCLUIR_COHORTE <- ""                                   # "" keeps all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""  # suffix only when a cohort is dropped
GSE_IDS <- setdiff(GSE_IDS, EXCLUIR_COHORTE)            # no-op with ""
cat("Cohortes en esta corrida:", length(GSE_IDS), "|", paste(GSE_IDS, collapse=", "),
    if(nchar(EXCLUIR_COHORTE)) paste0(" (excluida ", EXCLUIR_COHORTE, ")") else "", "\n")

# Optional method comparisons (3b-3e); not needed for the normalized output and memory-heavy.
CORRER_DIAGNOSTICOS <- FALSE   # TRUE runs blocks 3b-3e

# P1 reads and writes unsuffixed files; P2/P3 use suffixed files and folders.
VARIANTE_MASKING <- "P2"                                # P2 = registered probe panel
SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)

DIR_CURATED    <- "../data_curated"                     # per-cohort M/U from 02
DIR_NORMALIZED <- paste0("../data_normalized", SUF_VAR, SUF_COH)   # output of this script
DIR_PHENO      <- "../data_pheno"                       # phenotype tables (diagnostics only)

dir.create(DIR_NORMALIZED, showWarnings = FALSE, recursive = TRUE)   # created if absent


# ---- 3. Load curated cohorts and assemble the joint matrix ------------------
# Joint quantile normalization needs the same probe set in every sample, so each cohort is
# reduced to the probes shared by all cohorts. Loop 1 finds the intersection; loop 2 stacks
# the trimmed M and U matrices column-wise. Reading twice keeps only one cohort in memory at a time.

sondas_comunes <- NULL                                  # running intersection of probe IDs
for (gse in GSE_IDS) {
  # Curated file of this masking variant.
  ruta <- file.path(DIR_CURATED, paste0(gse, "_curated", SUF_VAR, ".rds"))
  x    <- readRDS(ruta)                                  # list with $M and $U (probes x samples)
  sgse <- rownames(x$M)                                  # probes kept by 02 in this cohort
  sondas_comunes <- if (is.null(sondas_comunes)) sgse
                    else intersect(sondas_comunes, sgse)   # first cohort seeds the set
  message(gse, ": ", length(sgse), " sondas | intersección acumulada: ",   # progress per cohort
          length(sondas_comunes))
  rm(x); gc()                                            # free before the next cohort
}
message("Sondas comunes a las 6 cohortes: ", length(sondas_comunes))


M_list <- vector("list", length(GSE_IDS))                # methylated signal per cohort
U_list <- vector("list", length(GSE_IDS))                # unmethylated signal per cohort
names(M_list) <- GSE_IDS; names(U_list) <- GSE_IDS   # index lists by GSE ID
cohorte_de_muestra <- character(0)                       # cohort of each sample

for (gse in GSE_IDS) {
  # Re-read, now to keep the trimmed matrices.
  x <- readRDS(file.path(DIR_CURATED, paste0(gse, "_curated", SUF_VAR, ".rds")))
  M_list[[gse]] <- x$M[sondas_comunes, , drop = FALSE]   # subset and reorder to the common probes
  U_list[[gse]] <- x$U[sondas_comunes, , drop = FALSE]   # same rows, same order as M
  cohorte_de_muestra <- c(cohorte_de_muestra,
                          setNames(rep(gse, ncol(x$M)), colnames(x$M)))   # named by GSM ID
  message(gse, ": ", ncol(x$M), " muestras añadidas")   # progress per cohort
  rm(x); gc()
}

M <- do.call(cbind, M_list)                              # probes x all samples
U <- do.call(cbind, U_list)   # U stacked in the same column order
rm(M_list, U_list); gc()   # free the per-cohort copies


# Checks dimensions, M/U alignment and absence of NA (02 filters without masking to NA).
cat("Dimensiones M:", nrow(M), "sondas ×", ncol(M), "muestras\n")
cat("Dimensiones U:", nrow(U), "sondas ×", ncol(U), "muestras\n")
cat("¿M y U alineadas (sondas y muestras)?: ",
    identical(dimnames(M), dimnames(U)), "\n")           # must be TRUE

n_na <- sum(is.na(M))                                    # any NA would be a reading bug
cat("NA en M:", n_na, sprintf("(%.4f%% de las celdas) — esperado 0\n", 100 * n_na / length(M)))
cat("Sondas con ≥1 NA:", sum(matrixStats::rowAnyNAs(M)),
    "| Muestras con ≥1 NA:", sum(matrixStats::colAnyNAs(M)), "\n")

cat("\nMuestras por cohorte (orden de columnas):\n")
print(table(cohorte_de_muestra[colnames(M)]))   # column order follows GSE_IDS


# ---- 3b. Optional: raw beta with offsets 0, 80 and 100 ----------------------
# Method comparison only. Computes raw beta = M / (M + U + alpha) for three offsets and the raw
# SCZ - CNT delta beta for each. A second part compares the tail of raw |delta beta| (offset 100)
# between probes kept under P2 but excluded under P1 (Naeem list, Probe_rs) and the remaining
# core probes. Delta beta here uses all samples: a check of scale, not the inference of 06.
if (CORRER_DIAGNOSTICOS) {
GUARDAR_BETAS_CRUDAS <- FALSE                            # TRUE also writes each raw beta matrix
offsets <- c(sin_offset = 0, off80 = 80, off100 = 100)   # 100 = Illumina/minfi/dasen convention

pheno_d  <- read.csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                     stringsAsFactors = FALSE)           # inclusive phenotype table used by 02
grp      <- pheno_d$group[match(colnames(M), pheno_d$gsm_id)]   # SCZ/CNT aligned to columns
scz      <- which(grp == "SCZ"); cnt <- which(grp == "CNT")     # column indices per group
# Group sizes; samples without a group are left out of delta beta.
cat("Muestras con grupo para Δβ — SCZ:", length(scz), "| CNT:", length(cnt),
    "| sin grupo:", sum(is.na(grp)), "\n\n")

tot      <- M + U                                        # total intensity per cell
lo_idx   <- tot <= quantile(tot, 0.05, na.rm = TRUE)     # low-intensity cells (lowest 5%)
dbeta_formas <- vector("list", length(offsets))          # delta beta vector per offset
names(dbeta_formas) <- names(offsets)
diag_filas <- vector("list", length(offsets))            # one summary row per offset

for (k in seq_along(offsets)) {
  a    <- offsets[[k]]   # offset value
  nm   <- names(offsets)[k]   # label used in names and files
  beta <- M / (tot + a)                                  # raw beta with this offset

  dbeta <- rowMeans(beta[, scz, drop = FALSE]) -
           rowMeans(beta[, cnt, drop = FALSE])           # raw SCZ - CNT difference per probe
  dbeta_formas[[nm]] <- dbeta   # kept for the 80 vs 100 comparison

  # Distribution summary: spread, saturation near 0/1, behavior at low intensity.
  diag_filas[[nm]] <- data.frame(
    forma        = nm, offset = a,
    media        = mean(beta),
    sd           = sd(beta),
    min          = min(beta), max = max(beta),
    pct_lt_0.01  = 100 * mean(beta < 0.01),
    pct_gt_0.99  = 100 * mean(beta > 0.99),
    media_baja_int = mean(beta[lo_idx]),                 # the offset matters most here
    sd_dbeta     = sd(dbeta), max_abs_dbeta = max(abs(dbeta)),
    stringsAsFactors = FALSE)

  cat(sprintf("forma %-10s (α=%3d): sd(β)=%.5f | %%β<0.01=%.2f %%β>0.99=%.2f | sd(Δβ)=%.5f max|Δβ|=%.4f\n",
              nm, a, sd(beta), 100*mean(beta<0.01), 100*mean(beta>0.99), sd(dbeta), max(abs(dbeta))))

  if (GUARDAR_BETAS_CRUDAS)
    saveRDS(beta, file.path(DIR_NORMALIZED, paste0("betas_crudas_", nm, ".rds")))
  rm(beta); gc()                                         # one full beta matrix at a time
}

diag_tab <- do.call(rbind, diag_filas)   # one row per offset
cat("\nTabla comparativa de las 3 formas:\n"); print(diag_tab, row.names = FALSE)

# How far apart the offsets leave delta beta (80 vs 100 should be negligible).
cat(sprintf("\nDif. media |Δβ| entre offset 80 y 100: %.6f (debería ser diminuta)\n",
            mean(abs(dbeta_formas[["off80"]] - dbeta_formas[["off100"]]))))
cat(sprintf("Dif. media |Δβ| entre offset 0 y 100:  %.6f\n",
            mean(abs(dbeta_formas[["sin_offset"]] - dbeta_formas[["off100"]]))))

saveRDS(list(dbeta = dbeta_formas, diagnostico = diag_tab, offsets = offsets),
        file.path(DIR_NORMALIZED, "betas_crudas_diagnostico.rds"))   # reused in 3e
cat("\nGuardado: betas_crudas_diagnostico.rds (3 vectores Δβ + tabla)\n")

# Tail of raw |delta beta|: probes excluded under P1 but present in this panel vs core probes.
dbeta_off100 <- dbeta_formas[["off100"]]                 # standard offset
p1_excl <- unique(c(   # union of the P1-only exclusions
  tryCatch(readLines(file.path("..", "refs", "naeem2014_problematic_450k.txt")), warning = function(w) character(0)),
  { a <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19); rownames(a)[!is.na(a$Probe_rs)] }))
recuperadas <- intersect(rownames(M), p1_excl)           # probes P1 would have removed
core        <- setdiff(rownames(M), p1_excl)             # probes kept by both variants
cat(sprintf("\n=== Δβ MÁXIMA — recuperadas (P1 las botaba) vs core (%d vs %d sondas) ===\n",
            length(recuperadas), length(core)))
if (length(recuperadas) > 0) {                           # empty under P1
  # Tail statistics: maximum, upper quantiles and share of probes above two effect sizes.
  tail_stats <- function(v) c(max = max(abs(v), na.rm = TRUE),
                              q99.9 = quantile(abs(v), 0.999, na.rm = TRUE),
                              q99 = quantile(abs(v), 0.99, na.rm = TRUE),
                              pct_gt_0.035 = 100 * mean(abs(v) > 0.035, na.rm = TRUE),
                              pct_gt_0.02 = 100 * mean(abs(v) > 0.02, na.rm = TRUE))
  cmp <- rbind(recuperadas = tail_stats(dbeta_off100[recuperadas]),
               core        = tail_stats(dbeta_off100[core]))   # heavier tail = signal in removed probes
  print(round(cmp, 4))
  cat("Top 10 |Δβ| entre las RECUPERADAS:\n")
  # Largest effects among recovered probes.
  print(round(sort(abs(dbeta_off100[recuperadas]), decreasing = TRUE)[1:10], 4))
  saveRDS(list(recuperadas = recuperadas, core = core, cmp = cmp,
               dbeta_recuperadas = dbeta_off100[recuperadas]),
          file.path(DIR_NORMALIZED, paste0("diagnostico_masking_dbeta", SUF_VAR, ".rds")))
  cat("Guardado: diagnostico_masking_dbeta", SUF_VAR, ".rds\n", sep = "")
} else {
  cat("(Corriendo P1: no hay sondas recuperadas. Correr con VARIANTE_MASKING='P2' para este diagnóstico.)\n")
}
rm(tot, lo_idx); gc()   # free the intensity matrices
}

# ---- 3c. Optional: six wateRmelon beta methods ------------------------------
# Method comparison only. Without IDAT files, the methods that need red/green channels or control
# probes are not applicable; the wateRmelon functions that take M/U directly are compared instead
# (naten, nanet, nasen, danes, danet, betaqn), with the correlation of their delta beta vectors.
if (CORRER_DIAGNOSTICOS) {
GUARDAR_BETAS_METODOS <- FALSE                           # TRUE also writes each method's betas
metodos <- c("naten", "nanet", "nasen", "danes", "danet", "betaqn")   # M/U-native methods; no IDAT needed

anno2  <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
onetwo <- anno2[rownames(M), "Type"]; names(onetwo) <- rownames(M)   # probe design I/II
phe    <- read.csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"), stringsAsFactors = FALSE)
grp    <- phe$group[match(colnames(M), phe$gsm_id)]      # recomputed so the block stands alone
scz    <- which(grp == "SCZ"); cnt <- which(grp == "CNT")

# One delta beta vector per method.
dbeta_met  <- vector("list", length(metodos)); names(dbeta_met) <- metodos
diag_met   <- vector("list", length(metodos)); names(diag_met)  <- metodos   # one summary row per method
beta_naten <- NULL                                        # reused as input to betaqn

for (mt in metodos) {
  cat("método", mt, "... ")   # progress
  # naten: raw (offset only); nanet: QN of M and U; nasen: QN by probe type, no background;
  # danes/danet: dye-bias correction + QN; betaqn: QN on the beta scale.
  beta <- tryCatch(switch(mt,
           naten  = naten(M, U),
           nanet  = nanet(M, U),
           nasen  = nasen(M, U, onetwo),
           danes  = danes(M, U, onetwo),
           danet  = danet(M, U, onetwo),
           betaqn = betaqn(if (!is.null(beta_naten)) beta_naten else naten(M, U))),
          error = function(e) { cat("✗", conditionMessage(e), "\n"); NULL })   # skip a failing method
  if (is.null(beta)) next
  if (is(beta, "MethylSet")) beta <- getBeta(beta)       # some methods return a MethylSet
  if (mt == "naten") beta_naten <- beta                  # keep for betaqn

  dbeta <- rowMeans(beta[, scz, drop = FALSE]) - rowMeans(beta[, cnt, drop = FALSE])   # SCZ - CNT
  dbeta_met[[mt]] <- dbeta
  diag_met[[mt]]  <- data.frame(   # same summary columns as 3b
    metodo = mt, media = mean(beta, na.rm = TRUE), sd = sd(beta, na.rm = TRUE),
    min = min(beta, na.rm = TRUE), max = max(beta, na.rm = TRUE),
    pct_lt_0.01 = 100 * mean(beta < 0.01, na.rm = TRUE),
    pct_gt_0.99 = 100 * mean(beta > 0.99, na.rm = TRUE),
    sd_dbeta = sd(dbeta, na.rm = TRUE), max_abs_dbeta = max(abs(dbeta), na.rm = TRUE),
    stringsAsFactors = FALSE)
  cat(sprintf("sd(β)=%.5f | sd(Δβ)=%.5f | max|Δβ|=%.4f\n",
              sd(beta, na.rm = TRUE), sd(dbeta, na.rm = TRUE), max(abs(dbeta), na.rm = TRUE)))
  if (GUARDAR_BETAS_METODOS)
    saveRDS(beta, file.path(DIR_NORMALIZED, paste0("betas_metodo_", mt, ".rds")))
  if (mt != "naten") { rm(beta); gc() }                  # naten copy is still needed
}

diag_met_tab <- do.call(rbind, diag_met)   # one row per method
cat("\nComparativa de las 6 estrategias:\n"); print(diag_met_tab, row.names = FALSE)

ok <- names(dbeta_met)[!vapply(dbeta_met, is.null, logical(1))]   # methods that ran
if (length(ok) >= 2) {
  cm <- cor(do.call(cbind, dbeta_met[ok]), use = "complete.obs")  # 1 = identical delta beta
  cat("\nCorrelación de Δβ entre estrategias (1 = idénticas):\n"); print(round(cm, 3))
}

saveRDS(list(dbeta = dbeta_met, diagnostico = diag_met_tab, metodos = metodos),
        file.path(DIR_NORMALIZED, "betas_metodos_diagnostico.rds"))   # reused in 3e
cat("\nGuardado: betas_metodos_diagnostico.rds (6 vectores Δβ + tabla + correlaciones)\n")
rm(beta_naten); gc()
}

# ---- 3d. Optional: minfi on a MethylSet rebuilt from M/U --------------------
# Method comparison only. Compares getBeta(type = "Illumina") and preprocessQuantile; the
# declared sex is passed because the sex-chromosome probes were removed in 02. Uses phe, scz
# and cnt from 3c.
if (CORRER_DIAGNOSTICOS) {
GUARDAR_BETAS_MINFI <- FALSE                             # TRUE also writes the minfi betas
ms <- MethylSet(Meth = M, Unmeth = U,   # minfi container from M/U only
                annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))
cat("MethylSet reconstruido:", class(ms), "|", nrow(ms), "×", ncol(ms), "\n")

sexo_decl <- phe$sexo[match(colnames(M), phe$gsm_id)]    # declared sex, aligned to columns

dbeta_minfi <- list(); diag_minfi <- list()   # results per minfi method
for (mt in c("minfi_illumina", "minfi_quantile")) {
  cat("método", mt, "... ")
  # minfi_illumina: raw beta with offset 100 (cross-package check of naten);
  # minfi_quantile: stratified quantile normalization by probe type and region.
  beta <- tryCatch(switch(mt,
           minfi_illumina = getBeta(ms, type = "Illumina"),
           minfi_quantile = getBeta(preprocessQuantile(ms, sex = sexo_decl))),
          error = function(e) { cat("✗", conditionMessage(e), "\n"); NULL })
  if (is.null(beta)) next
  beta  <- beta[, colnames(M), drop = FALSE]            # preprocessQuantile may reorder samples
  dbeta <- rowMeans(beta[, scz, drop = FALSE]) - rowMeans(beta[, cnt, drop = FALSE])   # SCZ - CNT
  dbeta_minfi[[mt]] <- dbeta
  # Same summary columns as 3b/3c.
  diag_minfi[[mt]]  <- data.frame(metodo = mt, media = mean(beta, na.rm = TRUE), sd = sd(beta, na.rm = TRUE),
    min = min(beta, na.rm = TRUE), max = max(beta, na.rm = TRUE),
    pct_lt_0.01 = 100 * mean(beta < 0.01, na.rm = TRUE), pct_gt_0.99 = 100 * mean(beta > 0.99, na.rm = TRUE),
    sd_dbeta = sd(dbeta, na.rm = TRUE), max_abs_dbeta = max(abs(dbeta), na.rm = TRUE), stringsAsFactors = FALSE)
  cat(sprintf("sd(β)=%.5f | sd(Δβ)=%.5f | max|Δβ|=%.4f\n",
              sd(beta, na.rm = TRUE), sd(dbeta, na.rm = TRUE), max(abs(dbeta), na.rm = TRUE)))
  if (GUARDAR_BETAS_MINFI) saveRDS(beta, file.path(DIR_NORMALIZED, paste0("betas_", mt, ".rds")))
  rm(beta); gc()
}
diag_minfi_tab <- do.call(rbind, diag_minfi)   # one row per method
cat("\nComparativa minfi (MethylSet reconstruido):\n"); print(diag_minfi_tab, row.names = FALSE)

# Correlation of delta beta across wateRmelon (3c) and minfi methods.
todos <- c(dbeta_met, dbeta_minfi)   # wateRmelon and minfi together
todos <- todos[!vapply(todos, is.null, logical(1))]      # drop methods that failed
if (length(todos) >= 2) {
  cm2 <- cor(do.call(cbind, todos), use = "complete.obs")
  cat("\nCorrelación de Δβ — TODAS las estrategias (wateRmelon + minfi):\n"); print(round(cm2, 3))
}
saveRDS(list(dbeta = dbeta_minfi, diagnostico = diag_minfi_tab),
        file.path(DIR_NORMALIZED, "betas_minfi_diagnostico.rds"))   # reused in 3e
cat("\nGuardado: betas_minfi_diagnostico.rds\n")
rm(ms); gc()
}

# ---- 3e. Optional: lumi quantile normalization of M-values ------------------
# Method comparison only. normalizeMethylation.quantile normalizes the M-value distribution as a
# whole, without separating M/U or probe types. The block ends with a correlation table of delta
# beta across all methods from 3b-3e plus dasen and dasen + BMIQ (read from their saved files).
if (CORRER_DIAGNOSTICOS) {
GUARDAR_BETAS_LUMI <- FALSE                              # TRUE also writes the lumi betas
# Loaded here only: not needed elsewhere.
suppressMessages(suppressWarnings({ library(lumi); library(methylumi) }))

# Phenotype and group indices recomputed so the block stands alone.
phe <- read.csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"), stringsAsFactors = FALSE)
grp <- phe$group[match(colnames(M), phe$gsm_id)]; scz <- which(grp == "SCZ"); cnt <- which(grp == "CNT")

dbeta_lumi <- NULL; dbeta_lumi_crudo <- NULL; diag_lumi <- NULL   # stay NULL if lumi fails
res_lumi <- tryCatch({
  mval <- log2((M + 1) / (U + 1))                       # +1 avoids log(0)
  # Assay data for the lumi object.
  adL  <- as.environment(list(methylated = M, unmethylated = U, exprs = mval)); rm(mval); gc()
  esL  <- new("MethyLumiM", assayData = adL)            # lumi container built from M/U
  b_cr <- methylumi::betas(esL)[rownames(M), colnames(M)]      # without normalization
  qnL  <- normalizeMethylation.quantile(esL); rm(esL); gc()    # QN of the joint M-value distribution
  b_qn <- methylumi::betas(qnL)[rownames(M), colnames(M)]; rm(qnL); gc()   # after quantile normalization
  list(crudo = b_cr, qn = b_qn)   # raw and normalized betas
}, error = function(e) { cat("lumi ✗", conditionMessage(e), "\n"); NULL })

if (!is.null(res_lumi)) {
  diag_lumi <- list()
  for (nv in c("crudo", "qn")) {                         # raw and normalized lumi betas
    b  <- res_lumi[[nv]]; nm <- paste0("lumi_", nv)   # matrix and its label
    db <- rowMeans(b[, scz, drop = FALSE]) - rowMeans(b[, cnt, drop = FALSE])   # SCZ - CNT
    if (nv == "crudo") dbeta_lumi_crudo <- db else dbeta_lumi <- db   # store under the right name
    # Same summary columns as 3b-3d.
    diag_lumi[[nm]] <- data.frame(metodo = nm, media = mean(b, na.rm = TRUE), sd = sd(b, na.rm = TRUE),
      min = min(b, na.rm = TRUE), max = max(b, na.rm = TRUE),
      pct_lt_0.01 = 100 * mean(b < 0.01, na.rm = TRUE), pct_gt_0.99 = 100 * mean(b > 0.99, na.rm = TRUE),
      sd_dbeta = sd(db, na.rm = TRUE), max_abs_dbeta = max(abs(db), na.rm = TRUE), stringsAsFactors = FALSE)
    cat(sprintf("%-11s sd(β)=%.5f | sd(Δβ)=%.5f | max|Δβ|=%.4f\n", nm,
                sd(b, na.rm = TRUE), sd(db, na.rm = TRUE), max(abs(db), na.rm = TRUE)))
    if (GUARDAR_BETAS_LUMI) saveRDS(b, file.path(DIR_NORMALIZED, paste0("betas_", nm, ".rds")))
  }
  saveRDS(list(dbeta = list(lumi_crudo = dbeta_lumi_crudo, lumi_qn = dbeta_lumi),
               diagnostico = do.call(rbind, diag_lumi)),
          file.path(DIR_NORMALIZED, "betas_lumi_diagnostico.rds"))
  cat("Guardado: betas_lumi_diagnostico.rds\n")
}
rm(res_lumi); gc()

# Combined comparison of all methods, aligned by CpG.
d1 <- readRDS(file.path(DIR_NORMALIZED, "betas_crudas_diagnostico.rds"))$dbeta[c("off80", "off100")]  # 3b
d2 <- readRDS(file.path(DIR_NORMALIZED, "betas_metodos_diagnostico.rds"))$dbeta   # 3c
d3 <- readRDS(file.path(DIR_NORMALIZED, "betas_minfi_diagnostico.rds"))$dbeta     # 3d

dbeta_desde <- function(ruta) {                          # delta beta from a saved beta matrix
  if (!file.exists(ruta)) return(NULL)                   # absent before steps 4-6 have run
  bb <- readRDS(ruta); if (is.list(bb) && !is.null(bb$betas)) bb <- bb$betas   # list or matrix
  bb <- bb[, colnames(M), drop = FALSE]                  # same column order as M
  # Returns the SCZ - CNT vector.
  db <- rowMeans(bb[, scz, drop = FALSE]) - rowMeans(bb[, cnt, drop = FALSE]); rm(bb); gc(); db
}
d_dasen <- dbeta_desde(file.path(DIR_NORMALIZED, "betas_dasen.rds"))        # dasen checkpoint
d_prod  <- dbeta_desde(file.path(DIR_NORMALIZED, "betas_normalized.rds"))   # dasen + BMIQ output

alld <- c(d1, d2, d3,   # every method that produced a vector
          if (!is.null(dbeta_lumi_crudo)) list(lumi_crudo = dbeta_lumi_crudo),
          if (!is.null(dbeta_lumi))       list(lumi_qn = dbeta_lumi),
          if (!is.null(d_dasen)) list(dasen = d_dasen),
          if (!is.null(d_prod))  list(dasen_bmiq = d_prod))
alld <- alld[!vapply(alld, is.null, logical(1))]         # keep available methods only
comun <- Reduce(intersect, lapply(alld, names))          # CpGs present in every vector

# Raw (offset only) vs between-array quantile-normalized methods.
crudas <- c("off80", "off100", "naten", "minfi_illumina", "lumi_crudo")   # offset-only methods
nivel  <- ifelse(names(alld) %in% crudas, "cruda", "QN")   # label each method raw or QN

cat("\n=== Δβ: comparación FINAL de las", length(alld), "estrategias (", length(comun), "CpG) ===\n")
resumen_final <- data.frame(   # one row per method
  estrategia = names(alld), nivel = nivel,
  sd_dbeta   = sapply(alld, function(v) sd(v[comun], na.rm = TRUE)),        # spread of delta beta
  max_abs_db = sapply(alld, function(v) max(abs(v[comun]), na.rm = TRUE)), row.names = NULL)
# QN group first, larger spread first.
print(resumen_final[order(resumen_final$nivel, -resumen_final$sd_dbeta), ], row.names = FALSE)

cm_all <- cor(sapply(alld, function(v) v[comun]), use = "complete.obs")    # all methods vs all
cat("\nCorrelación de Δβ — TODAS las estrategias (alineado por CpG):\n"); print(round(cm_all, 3))

# Fair comparison: only methods that include between-array quantile normalization.
qn_fam <- intersect(c("nasen", "danes", "danet", "nanet", "betaqn", "minfi_quantile",
                      "lumi_qn", "dasen", "dasen_bmiq"), names(alld))
if (length(qn_fam) >= 2) {
  cat("\n>>> Comparación JUSTA — SOLO métodos con QN inter-array (lumi_qn vs sus pares):\n")
  print(round(cm_all[qn_fam, qn_fam], 3))
}
saveRDS(list(correlacion = cm_all, resumen = resumen_final, nivel = nivel),
        file.path(DIR_NORMALIZED, "betas_comparacion_final.rds"))
cat("\nGuardado: betas_comparacion_final.rds\n")
}

# ---- 4. Joint between-array quantile normalization (dasen) ------------------
# dasen applies background adjustment and quantile normalization across all samples at once,
# separately for type I and type II probes, on M and U. It keeps the compression of type II
# values, which BMIQ corrects in step 5. It must be joint: per-cohort normalization would leave
# each cohort on its own reference distribution.

# 4a. Probe design vector aligned to the rows of M.
anno   <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)   # 450K manifest (hg19)
onetwo <- anno[rownames(M), "Type"]                      # "I" or "II" per probe
names(onetwo) <- rownames(M)

cat("¿onetwo alineado con las filas de M?: ",
    identical(names(onetwo), rownames(M)), "\n")         # must be TRUE
cat("Sondas sin tipo asignado (NA):", sum(is.na(onetwo)), "\n")   # expected 0
cat("Reparto de sondas por diseño:\n"); print(table(onetwo))   # counts of type I and type II


# 4b. dasen with a checkpoint: reloads betas_dasen.rds when present unless FORCE_DASEN is TRUE.
FORCE_DASEN <- FALSE                                     # TRUE recomputes even if the file exists
ckpt_dasen  <- file.path(DIR_NORMALIZED, "betas_dasen.rds")   # path of the dasen checkpoint

if (!FORCE_DASEN && file.exists(ckpt_dasen)) {
  beta_qn <- readRDS(ckpt_dasen)                         # skip the slowest step on reruns
  cat("beta_qn RECARGADO del checkpoint:", ckpt_dasen, "\n")
} else {
  # The warning "Sentrix position model failed" is expected: column names are GSM IDs, so the
  # chip-position background step is skipped and the type-wise quantile normalization still runs.
  beta_qn <- dasen(M, U, onetwo)                      # fudge = 100 (default)
  saveRDS(beta_qn, ckpt_dasen)                        # checkpoint before BMIQ
  cat("dasen calculado y guardado en checkpoint:", ckpt_dasen, "\n")
}

# Sanity checks: range within [0, 1], no NA, same probes and samples as the input.
cat("Dimensiones beta_qn:", nrow(beta_qn), "×", ncol(beta_qn), "\n")
cat("Rango de betas: [", round(min(beta_qn, na.rm = TRUE), 4), ",",
    round(max(beta_qn, na.rm = TRUE), 4), "]\n")
cat("NA en beta_qn:", sum(is.na(beta_qn)), "(esperado 0)\n")
cat("¿Mismas dimnames que M?: ", identical(dimnames(beta_qn), dimnames(M)), "\n")


# 4c. Type I vs type II beta densities before and after dasen, on 60 random samples.
# The type II peaks staying inside the type I peaks is what justifies BMIQ in step 5.
set.seed(7777777)                                        # same sample subset on every run
js     <- sort(sample(ncol(M), min(60, ncol(M))))        # columns to plot
is_I   <- onetwo == "I"                                  # logical row masks by design
is_II  <- onetwo == "II"

beta_raw_sub <- M[, js] / (M[, js] + U[, js] + 100)   # raw beta, subset only

options(repr.plot.width = 12, repr.plot.height = 5)    # notebook plot size
par(mfrow = c(1, 2))                                     # before | after
yl <- c(0, 11)                                           # common y-axis
plot(density(beta_raw_sub[is_I, ], na.rm = TRUE), col = "#1b9e77", lwd = 2,
     main = "Antes de dasen (crudo)", xlab = "beta", xlim = c(0, 1), ylim = yl)
lines(density(beta_raw_sub[is_II, ], na.rm = TRUE), col = "#d95f02", lwd = 2, lty = 2)
legend("topright", c("tipo I", "tipo II"), col = c("#1b9e77", "#d95f02"),
       lwd = 2, lty = c(1, 2), bty = "n")
plot(density(beta_qn[is_I, js], na.rm = TRUE), col = "#1b9e77", lwd = 2,
     main = "Después de dasen", xlab = "beta", xlim = c(0, 1), ylim = yl)
lines(density(beta_qn[is_II, js], na.rm = TRUE), col = "#d95f02", lwd = 2, lty = 2)
legend("topright", c("tipo I", "tipo II"), col = c("#1b9e77", "#d95f02"),
       lwd = 2, lty = c(1, 2), bty = "n")
par(mfrow = c(1, 1))


# 4d. Per-sample beta densities before and after dasen, colored by cohort.
cohortes     <- sort(unique(cohorte_de_muestra))   # cohort IDs in fixed order
paleta       <- c("#1b9e77", "#d95f02", "#7570b3",
                  "#e7298a", "#66a61e", "#e6ab02")       # one color per cohort
cohorte_cols <- setNames(paleta[seq_along(cohortes)], cohortes)   # cohort -> color
coh_js       <- cohorte_de_muestra[colnames(beta_qn)[js]]   # cohort of each plotted sample

options(repr.plot.width = 12, repr.plot.height = 5)
par(mfrow = c(1, 2))
for (k in seq_along(js)) {                               # raw: one curve per sample
  d <- density(beta_raw_sub[, k], na.rm = TRUE, from = 0, to = 1)   # density on [0, 1] for one sample
  if (k == 1) plot(d, col = cohorte_cols[coh_js[k]], lwd = 0.6,
                   main = "Antes de dasen (por muestra)", xlab = "beta", ylim = c(0, 5))
  else lines(d, col = cohorte_cols[coh_js[k]], lwd = 0.6)
}
legend("topright", names(cohorte_cols), col = cohorte_cols, lwd = 2, bty = "n", cex = 0.8)
for (k in seq_along(js)) {                               # after dasen: curves should overlap
  d <- density(beta_qn[, js[k]], na.rm = TRUE, from = 0, to = 1)
  if (k == 1) plot(d, col = cohorte_cols[coh_js[k]], lwd = 0.6,
                   main = "Después de dasen (por muestra)", xlab = "beta", ylim = c(0, 5))
  else lines(d, col = cohorte_cols[coh_js[k]], lwd = 0.6)
}
par(mfrow = c(1, 1))
rm(beta_raw_sub); gc()   # free the raw subset


# 4e. PCA of the 20,000 most variable probes after dasen, colored by cohort (between-cohort
# structure before harmonization).
vars   <- matrixStats::rowVars(beta_qn)                  # variance per probe
top    <- order(vars, decreasing = TRUE)[1:20000]        # most variable probes
pca    <- prcomp(t(beta_qn[top, ]), center = TRUE, scale. = FALSE)   # samples as rows
pct    <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)   # % variance per PC for axis labels

coh_all <- cohorte_de_muestra[colnames(beta_qn)]   # cohort of every sample
options(repr.plot.width = 8, repr.plot.height = 7)
plot(pca$x[, 1], pca$x[, 2], col = cohorte_cols[coh_all], pch = 19, cex = 0.4,
     xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
     main = "PCA tras dasen — color = cohorte (pre-ComBat)")   # cohort clusters expected here
legend("topright", names(cohorte_cols), col = cohorte_cols, pch = 19, bty = "n", cex = 0.8)

rm(M, U); gc()                                           # raw signals no longer needed


# ---- 5. Within-array BMIQ ---------------------------------------------------
# BMIQ (Teschendorff 2013) fits a three-component beta mixture per sample and maps type II
# values onto the type I scale. Runs in parallel over samples; if BMIQ fails for a sample, its
# dasen beta is kept and the sample is counted below.
design.v <- ifelse(onetwo == "I", 1, 2)          # BMIQ coding: 1 = type I, 2 = type II
sealseed(c("BMIQ"))                              # same seed on every BMIQ call: reproducible

# Forked workers can each end up holding a copy of the beta matrix; 4 keeps peak memory moderate.
n_workers <- 4   # parallel BMIQ workers
bp <- MulticoreParam(workers = n_workers, tasks = ncol(beta_qn), progressbar = TRUE)   # one task per sample

res <- bplapply(seq_len(ncol(beta_qn)), function(j) {   # one BMIQ fit per sample
  out <- tryCatch(BMIQ(beta_qn[, j], design.v, plots = FALSE, pri = FALSE)$nbeta,   # normalized betas
                  error = function(e) NULL)
  if (is.null(out)) beta_qn[, j] else out         # fallback: dasen beta
}, BPPARAM = bp)

beta_bmiq <- do.call(cbind, res)                  # back to probes x samples
dimnames(beta_bmiq) <- dimnames(beta_qn)          # restore probe and sample IDs
rm(res); gc()   # free the per-sample list

# Samples left unchanged by BMIQ are fallbacks.
n_fallback <- sum(matrixStats::colSums2(beta_bmiq != beta_qn) == 0)
cat("Dimensiones beta_bmiq:", nrow(beta_bmiq), "×", ncol(beta_bmiq), "\n")
cat("NA en beta_bmiq:", sum(is.na(beta_bmiq)), "(esperado 0)\n")
cat("Muestras con BMIQ fallido (fallback a dasen):", n_fallback, "\n")


# 5b. Type I vs type II densities after dasen and after BMIQ (same 60 samples).
# After BMIQ the type II peaks should line up with the type I peaks.
options(repr.plot.width = 12, repr.plot.height = 5)
par(mfrow = c(1, 2))                                     # dasen | dasen + BMIQ
yl <- c(0, 11)   # same y-axis as 4c
plot(density(beta_qn[is_I, js], na.rm = TRUE), col = "#1b9e77", lwd = 2,
     main = "Después de dasen", xlab = "beta", xlim = c(0, 1), ylim = yl)
lines(density(beta_qn[is_II, js], na.rm = TRUE), col = "#d95f02", lwd = 2, lty = 2)
legend("topright", c("tipo I", "tipo II"), col = c("#1b9e77", "#d95f02"),
       lwd = 2, lty = c(1, 2), bty = "n")
plot(density(beta_bmiq[is_I, js], na.rm = TRUE), col = "#1b9e77", lwd = 2,
     main = "Después de BMIQ", xlab = "beta", xlim = c(0, 1), ylim = yl)
lines(density(beta_bmiq[is_II, js], na.rm = TRUE), col = "#d95f02", lwd = 2, lty = 2)
legend("topright", c("tipo I", "tipo II"), col = c("#1b9e77", "#d95f02"),
       lwd = 2, lty = c(1, 2), bty = "n")
par(mfrow = c(1, 1))


# ---- 6. Save ----------------------------------------------------------------
# betas_normalized.rds: dasen + BMIQ beta matrix (probes x samples) and the cohort of each sample.
cohorte_out <- cohorte_de_muestra[colnames(beta_bmiq)]   # cohort labels in column order
saveRDS(list(betas = beta_bmiq, cohorte = cohorte_out),
        file.path(DIR_NORMALIZED, "betas_normalized.rds"))   # input to 04 (regRCPqn)
cat("Guardado: betas_normalized.rds —", nrow(beta_bmiq), "sondas ×",
    ncol(beta_bmiq), "muestras (dasen + BMIQ, sin NA)\n")
