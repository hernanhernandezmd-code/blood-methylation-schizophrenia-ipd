# 04 — Cross-cohort harmonization (regRCPqn) and epigenetic covariates
#
# Harmonizes the dasen + BMIQ betas across cohorts with regRCPqn (Sala et al., 2020) and estimates
# per-sample covariates: blood cell fractions (EpiDISH), smoking (EpiSmokEr) and epigenetic age
# (Horvath and Hannum, methylclock). Cell fractions and smoking use the pre-harmonization betas;
# epigenetic age uses the harmonized betas. Multiple imputation of missing age, the 80/20 partition
# and ComBat are left to script 05, so everything computed here is deterministic and leakage-free.
# Paper section: Normalization and Harmonization; Covariates
# Input:  ../data_normalized_P2/betas_normalized.rds, ../data_pheno/phenodata_combinado_imputar_edad.csv
# Output: ../data_harmonized_P2/betas_harmonized.rds, ../data_pheno/covariables_epigeneticas_P2.csv
# VARIANTE_MASKING: "P1" | "P2" | "P3" (default "P2"); must match scripts 02 and 03.
# EXCLUIR_COHORTE:  "" (all six cohorts, default) or a GSE ID to leave out.
# FORCE_REGRCPQN:   FALSE (reload the regRCPqn checkpoint if present, default) | TRUE (recompute).

# ---- 1. Utility functions and packages ---------------------------------------------------------

# Installs missing CRAN/Bioconductor packages and loads all of them.
InsLoad.pks <- function(pqks, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))                  # BiocManager installs both
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  ip <- rownames(installed.packages())                                   # packages already present
  new.pkg <- pqks[!(pqks %in% ip)]                                       # packages still missing
  if (length(new.pkg) > 0)
    BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)  # no upgrades of others
  out <- sapply(pqks, require, character.only = TRUE, quietly = TRUE)    # TRUE/FALSE per package
  print(out)
  if (length(new.pkg) > 0)
    message("Instalados ", length(new.pkg), " paquetes nuevos.")
  else
    message("Todos los paquetes requeridos ya estaban instalados.")
  invisible(out)
}


# Replaces each named function with a wrapper that calls set.seed() before every call.
# Applied to DNAmAge (step 5), whose KNN imputation of missing clock CpGs is stochastic.
sealseed <- function(func_names, seed = 7777777) {
  for (func_name in func_names) {
    original_func <- get(func_name, envir = .GlobalEnv)                  # function to be wrapped
    new_func <- function(...) {
      set.seed(as.integer(get("seed", envir = environment())))           # same seed on every call
      get("original_func", envir = environment())(...)                   # then the original call
    }
    environment(new_func) <- new.env(parent = environment(original_func))  # private closure
    assign("original_func", original_func, envir = environment(new_func))  # keeps the original
    assign("seed", seed, envir = environment(new_func))                    # keeps the seed
    assign(func_name, new_func, envir = .GlobalEnv)                        # masks the original name
  }
}


pqks_BC <- c(
  "minfi", "wateRmelon",
  "EpiDISH",                                             # cell fractions from betas (no IDAT files available)
  "methylclock", "methylclockData",                      # Horvath and Hannum clocks
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",       # region and probe type used by regRCPqn
  "preprocessCore", "impute",                            # quantile normalization; KNN for clock CpGs
  "BiocParallel"
)
pqks_CRAN <- c(
  "data.table", "matrixStats",
  "readr", "stringr", "ggplot2", "remotes"               # remotes installs the GitHub packages
)

InsLoad.pks(c(pqks_BC, pqks_CRAN))

# regRCPqn and EpiSmokEr are not on CRAN or Bioconductor; they are installed from GitHub.
gh_pkgs <- c(regRCPqn = "regRCPqn/regRCPqn",
             EpiSmokEr = "sailalithabollepalli/EpiSmokEr")
for (nm in names(gh_pkgs)) {
  if (!requireNamespace(nm, quietly = TRUE))
    remotes::install_github(gh_pkgs[[nm]], upgrade = "never")             # install only if missing
  library(nm, character.only = TRUE)
}

# Stops early if a required package failed to load, instead of failing later in step 5.
criticos <- c("minfi", "EpiDISH", "methylclock", "regRCPqn", "EpiSmokEr",
              "IlluminaHumanMethylation450kanno.ilmn12.hg19", "preprocessCore")
faltan <- criticos[!criticos %in% loadedNamespaces()]                    # critical but not loaded
if (length(faltan) > 0)
  stop("No se cargaron paquetes críticos: ", paste(faltan, collapse = ", "),
       ". Reinicia el kernel y vuelve a ejecutar esta celda.")
cat("Paquetes críticos cargados OK:", paste(criticos, collapse = ", "), "\n")


# ---- 2. Paths, parameters and data -------------------------------------------------------------

set.seed(7777777)                                        # global seed for the random steps below

# Runs from the repository root or from pipeline/.
if (!file.exists("../data_normalized") && dir.exists("fase_1/pipeline")) {
  setwd("fase_1/pipeline")
}

GSE_IDS  <- c("GSE84727", "GSE80417", "GSE147221",
              "GSE116378", "GSE116379", "GSE152027")
GSE_REF  <- "GSE84727"                                  # regRCPqn reference cohort (largest, n = 847)

VARIANTE_MASKING <- "P2"                                 # probe-masking variant read from script 03
SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)  # folder suffix
EXCLUIR_COHORTE <- ""                                    # "" keeps all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""   # suffix only when a cohort is left out
GSE_IDS <- setdiff(GSE_IDS, EXCLUIR_COHORTE)             # cohorts iterated in steps 4 and 5

DIR_NORMALIZED <- paste0("../data_normalized", SUF_VAR, SUF_COH)   # output of script 03
DIR_HARMONIZED <- paste0("../data_harmonized", SUF_VAR, SUF_COH)   # output of this script
DIR_PHENO      <- "../data_pheno"
DIR_REF        <- file.path(DIR_HARMONIZED, "regRCPqn_ref")   # reference distributions written by regRCPqn

dir.create(DIR_HARMONIZED, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_REF,        showWarnings = FALSE, recursive = TRUE)

norm   <- readRDS(file.path(DIR_NORMALIZED, "betas_normalized.rds"))   # list(betas, cohorte)
betas  <- norm$betas                                    # probes x samples, columns = GSM IDs
cohorte <- norm$cohorte                                 # cohort (GSE) of each sample
rm(norm); gc()                                          # frees the duplicate copy

# Phenotype table that keeps samples with missing age; rows reordered to the beta columns.
pheno <- readr::read_csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                         show_col_types = FALSE)
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$gsm_id                         # GSM ID as row key
pheno <- pheno[colnames(betas), ]                       # same order as the beta columns

# Two ages in GSE80417 carry a spurious extra digit and are mapped to the values below.
# Any other age above 100 is set to NA and imputed in script 05.
correc_edad <- c("891" = 90, "883" = 83)                 # recorded value -> corrected age
for (v in names(correc_edad)) {
  hit <- !is.na(pheno$edad) & pheno$edad == as.numeric(v)   # samples with that recorded value
  if (any(hit)) { cat("Edad corregida:", v, "→", correc_edad[[v]],
                      "(", sum(hit), "muestra/s)\n"); pheno$edad[hit] <- correc_edad[[v]] }
}
resto <- !is.na(pheno$edad) & pheno$edad > 100           # implausible ages left after the fix
if (any(resto)) { cat("Otras edades >100 puestas a NA:", paste(sort(pheno$edad[resto]), collapse=", "),
                     "\n"); pheno$edad[resto] <- NA }

# Console checks: dimensions, alignment, cohort sizes, sex and missing age.
cat("Dimensiones betas:", nrow(betas), "sondas ×", ncol(betas), "muestras\n")
cat("¿phenodata alineado con las columnas de betas?: ",
    identical(rownames(pheno), colnames(betas)), "\n")
cat("Muestras por cohorte:\n"); print(table(cohorte))
cat("Sexo (declarado):\n");     print(table(pheno$sexo, useNA = "ifany"))
cat("Edad cronológica faltante (NA):", sum(is.na(pheno$edad)),
    sprintf("(%.1f%%)\n", 100 * mean(is.na(pheno$edad))))


# ---- 3. Per-sample covariates (pre-harmonization betas) ----------------------------------------
# Both estimators work sample by sample against their own reference, so they use the betas from
# script 03, before regRCPqn.

# 3a. Cell fractions: EpiDISH, robust partial correlation (RPC) against the centDHSbloodDMC.m
# blood reference (7 cell types). minfi::estimateCellCounts cannot be used: it needs IDAT files.
data(centDHSbloodDMC.m, package = "EpiDISH")              # blood reference matrix
epidish_out <- epidish(beta.m = betas, ref.m = centDHSbloodDMC.m, method = "RPC")
cell_frac   <- as.data.frame(epidish_out$estF)         # samples x 7 cell types
cell_frac   <- cell_frac[colnames(betas), , drop = FALSE]   # same sample order as betas

# Console checks: reference coverage, mean fractions, and fractions summing to about 1.
cat("CpGs de referencia usados:", sum(rownames(centDHSbloodDMC.m) %in% rownames(betas)),
    "de", nrow(centDHSbloodDMC.m), "\n")
cat("Fracciones celulares (resumen, medias por tipo):\n")
print(round(colMeans(cell_frac), 4))
cat("¿Filas suman ~1?: rango", round(range(rowSums(cell_frac)), 4), "\n")


# 3b. Smoking: EpiSmokEr with method = "all" returns the smoking status class (SSt) and the
# continuous smoking score (SSc). The score is the covariate used in the models.
ss <- data.frame(row.names = colnames(betas),
                 sex = ifelse(pheno$sexo == "M", 1L, 2L))              # 1 = male, 2 = female
smoke_out <- tryCatch(
  epismoker(dataset = betas, samplesheet = ss, method = "all"),
  error = function(e) epismoker(dataset = betas, method = "all"))      # for versions without samplesheet
rownames(smoke_out) <- smoke_out$SampleName                            # one row per sample
smoke_out <- smoke_out[colnames(betas), , drop = FALSE]                # same order as betas

# The name of the score column differs between EpiSmokEr versions.
sc_col <- intersect(c("smokingScore", "SmokingScore", "Smoking.Score", "methylationScore"),
                    colnames(smoke_out))[1]                            # first known name present
if (is.na(sc_col)) sc_col <- grep("score", colnames(smoke_out), ignore.case = TRUE, value = TRUE)[1]
cat("Columnas EpiSmokEr:", paste(colnames(smoke_out), collapse = ", "), "\n")
cat("Score continuo (SSc) en columna:", sc_col, "\n")
cat("Estado de tabaquismo predicho (EpiSmokEr SSt, clase):\n")
print(table(smoke_out$PredictedSmokingStatus, useNA = "ifany"))
cat("Score continuo (SSc) — resumen:\n"); print(summary(smoke_out[[sc_col]]))


# ---- 4. Cross-cohort harmonization: regRCPqn ---------------------------------------------------
# Regional quantile normalization of M-values by genomic region (Island/Shore/Shelf/OpenSea) and
# probe type (I/II). The reference cohort builds the reference distribution; the other cohorts are
# aligned to it. The method does not use diagnosis and assumes no global methylation difference
# between groups, which holds for case-control blood.

# 4a. Betas to M-values
eps    <- 1e-6                                          # keeps betas away from 0 and 1 before the logit
betas_c <- pmin(pmax(betas, eps), 1 - eps)              # clipped copy
M_all  <- log2(betas_c / (1 - betas_c))                 # M = log2(beta / (1 - beta))
rm(betas_c); gc()

cat("M-values: rango [", round(min(M_all), 2), ",", round(max(M_all), 2), "]  NA:", sum(is.na(M_all)), "\n")


# 4b. regRCPqn per cohort, with a checkpoint that is reloaded unless FORCE_REGRCPQN is TRUE.
FORCE_REGRCPQN <- FALSE                                  # TRUE recomputes even if the checkpoint exists
ckpt_harm <- file.path(DIR_HARMONIZED, "betas_harmonized_checkpoint.rds")

if (!FORCE_REGRCPQN && file.exists(ckpt_harm)) {
  ck <- readRDS(ckpt_harm)                               # list(M, betas) from an earlier run
  M_harm <- ck$M; betas_harm <- ck$betas; cpg_comunes <- rownames(M_harm)
  cat("regRCPqn RECARGADO del checkpoint:", ckpt_harm, "\n")
} else {
  M_norm_list <- vector("list", length(GSE_IDS))         # one harmonized matrix per cohort
  names(M_norm_list) <- GSE_IDS

  for (gse in c(GSE_REF, setdiff(GSE_IDS, GSE_REF))) {  # reference cohort first
    cols <- colnames(M_all)[cohorte == gse]             # samples of this cohort
    df   <- data.frame(ID_REF = rownames(M_all),       # regRCPqn expects ID_REF as the first column
                       M_all[, cols, drop = FALSE], check.names = FALSE)

    if (gse == GSE_REF) {
      # Reference cohort: normalizes and writes the reference distributions to DIR_REF.
      out <- regRCPqn(M_data = df, ref_path = DIR_REF, data_name = "metaSCZ", save_ref = TRUE)
    } else {
      # Other cohorts: aligned to the stored reference distributions.
      out <- regRCPqnREF(M_data = df, ref_path = DIR_REF, data_name = "metaSCZ")
    }

    # The output has CpGs as rownames, reordered by position; rows are aligned below by CpG name.
    M_norm_list[[gse]] <- as.matrix(out[, cols, drop = FALSE])
    message(gse, ": regRCPqn OK (", length(cols), " muestras, ", nrow(out), " CpGs)")
    rm(df, out); gc()
  }

  cpg_comunes <- Reduce(intersect, lapply(M_norm_list, rownames))   # CpGs returned for every cohort
  M_harm <- do.call(cbind, lapply(M_norm_list, function(m) m[cpg_comunes, , drop = FALSE]))
  M_harm <- M_harm[, colnames(M_all)]                  # original sample order
  rm(M_norm_list); gc()

  betas_harm <- 2^M_harm / (2^M_harm + 1)                # back to the beta scale
  saveRDS(list(M = M_harm, betas = betas_harm), ckpt_harm)
  cat("regRCPqn calculado y guardado en checkpoint:", ckpt_harm, "\n")
}

# Console check: size, NA, CpGs kept and sample order.
cat("M_harm:", nrow(M_harm), "×", ncol(M_harm), " | NA:", sum(is.na(M_harm)),
    " | CpGs conservados:", length(cpg_comunes), "de", nrow(M_all),
    " | ¿orden de muestras intacto?:", identical(colnames(M_harm), colnames(M_all)), "\n")


# 4c. Per-sample M-value densities (60 random samples) coloured by cohort, before and after regRCPqn.
# After harmonization the cohort curves should overlap more while keeping the bimodal shape.
set.seed(7777777)
js  <- sort(sample(ncol(M_all), min(60, ncol(M_all))))   # random subset of samples
pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
coh_cols <- setNames(pal[seq_along(GSE_IDS)], GSE_IDS)   # one colour per cohort
coh_js   <- cohorte[js]                                  # cohort of each sampled column

options(repr.plot.width = 12, repr.plot.height = 5)
par(mfrow = c(1, 2))                                     # left = before, right = after
started <- FALSE                                         # first curve opens the plot
for (k in seq_along(js)) {
  x <- M_all[, js[k]]; if (sum(is.finite(x)) < 2) next   # skips columns without usable values
  d <- density(x, na.rm = TRUE)
  if (!started) { plot(d, col = coh_cols[coh_js[k]], lwd = 0.6, main = "Antes de regRCPqn (M)",
                       xlab = "M-value", xlim = c(-6, 6)); started <- TRUE }
  else lines(d, col = coh_cols[coh_js[k]], lwd = 0.6)    # later curves are overlaid
}
legend("topright", names(coh_cols), col = coh_cols, lwd = 2, bty = "n", cex = 0.8)
started <- FALSE                                         # same loop on the harmonized M-values
for (k in seq_along(js)) {
  x <- M_harm[, js[k]]; if (sum(is.finite(x)) < 2) next
  d <- density(x, na.rm = TRUE)
  if (!started) { plot(d, col = coh_cols[coh_js[k]], lwd = 0.6, main = "Después de regRCPqn (M)",
                       xlab = "M-value", xlim = c(-6, 6)); started <- TRUE }
  else lines(d, col = coh_cols[coh_js[k]], lwd = 0.6)
}
par(mfrow = c(1, 1))


# 4d. PCA on the 20,000 most variable probes of the harmonized betas, coloured by cohort.
# Baseline of the cohort structure left after regRCPqn, to compare with ComBat in script 05.
vars <- matrixStats::rowVars(betas_harm)                 # variance per probe
top  <- order(vars, decreasing = TRUE)[1:20000]          # 20,000 most variable probes
pca  <- prcomp(t(betas_harm[top, ]), center = TRUE, scale. = FALSE)   # samples as rows
pct  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)     # % variance per component

options(repr.plot.width = 8, repr.plot.height = 7)
plot(pca$x[, 1], pca$x[, 2], col = coh_cols[cohorte], pch = 19, cex = 0.4,
     xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
     main = "PCA tras regRCPqn — color = cohorte (pre-ComBat)")
legend("topright", names(coh_cols), col = coh_cols, pch = 19, bty = "n", cex = 0.8)


# ---- 5. Epigenetic age (harmonized betas) ------------------------------------------------------
# Chronological-age clocks only (Horvath and Hannum); PhenoAge and GrimAge are not used because
# they are accelerated in schizophrenia and would bias the imputed age of cases. Age is computed
# after regRCPqn and before ComBat, which can distort clock CpGs. Multiple imputation and
# calibration of missing age are done in script 05.

# 5a. Horvath and Hannum with methylclock::DNAmAge.
sealseed(c("DNAmAge"))                                   # reproducible KNN imputation of clock CpGs
# min.perc is the minimum fraction of clock CpGs that must be present (package default 0.8).
# Missing clock CpGs are imputed by KNN; steps 5b-5d check whether the resulting ages are usable.
MIN_PERC_RELOJ <- 0.4
dnam <- DNAmAge(betas_harm, clocks = c("Horvath", "Hannum"),
                cell.count = FALSE, min.perc = MIN_PERC_RELOJ)  # cell fractions come from step 3a
dnam <- as.data.frame(dnam)
rownames(dnam) <- dnam$id                                # sample ID as row key
dnam <- dnam[colnames(betas_harm), , drop = FALSE]       # same order as the betas

cat("Edad epigenética (DNAm) — resumen:\n")
print(summary(dnam[, c("Horvath", "Hannum")]))


# 5b. Uncalibrated MAE and correlation with recorded age, per cohort, in samples with known age.
# The calibrated, leakage-free MAE is reported in script 05.
edad_real <- pheno$edad                                  # recorded chronological age, with NA
cat(sprintf("%-11s %6s %8s %8s %8s\n", "cohorte", "n", "MAE_Hor", "MAE_Han", "r_Hannum"))
for (gse in GSE_IDS) {
  sel <- cohorte == gse & !is.na(edad_real)              # samples of this cohort with known age
  n   <- sum(sel)
  if (n < 3) { cat(sprintf("%-11s %6d  (sin edad suficiente)\n", gse, n)); next }
  mae_h <- mean(abs(dnam$Horvath[sel] - edad_real[sel]))   # mean absolute error, Horvath
  mae_n <- mean(abs(dnam$Hannum[sel]  - edad_real[sel]))   # mean absolute error, Hannum
  r_n   <- cor(dnam$Hannum[sel], edad_real[sel])           # Pearson correlation, Hannum
  cat(sprintf("%-11s %6d %8.2f %8.2f %8.3f\n", gse, n, mae_h, mae_n, r_n))
}


# 5c. Case-control deviation of each clock: lm(DNAm age ~ group + recorded age).
# A clock usable for imputing the missing ages of cases must be blind to diagnosis.
df_bias <- data.frame(Hannum  = dnam$Hannum, Horvath = dnam$Horvath,
                      grupo   = factor(pheno$group), edad = edad_real, cohorte = cohorte)
df_bias <- df_bias[!is.na(df_bias$edad) & df_bias$grupo %in% c("SCZ", "CNT"), ]   # known age only
df_bias$grupo <- droplevels(df_bias$grupo)

for (reloj in c("Hannum", "Horvath")) {
  fit <- lm(df_bias[[reloj]] ~ grupo + edad, data = df_bias)
  co  <- summary(fit)$coefficients
  fila <- grep("^grupo", rownames(co))                   # row of the group coefficient
  cat(sprintf("%s — desviación caso-control (β grupo): %+.2f años (p = %.3f)\n",
              reloj, co[fila, 1], co[fila, 4]))
}
cat("Esperado: |β| pequeño y no significativo → reloj ciego al diagnóstico (E7).\n")


# 5d. Sensitivity: clocks recomputed on the pre-regRCPqn betas, compared with step 5a on MAE per
# cohort and on the case-control coefficient. Similar values mean regRCPqn does not distort age.
dnam_pre <- DNAmAge(betas, clocks = c("Horvath", "Hannum"), cell.count = FALSE, min.perc = MIN_PERC_RELOJ)
dnam_pre <- as.data.frame(dnam_pre); rownames(dnam_pre) <- dnam_pre$id
dnam_pre <- dnam_pre[colnames(betas), , drop = FALSE]    # same order as the betas

# MAE per cohort, pre versus post regRCPqn.
cat(sprintf("%-11s %6s | %8s %8s | %8s %8s\n", "cohorte", "n",
            "MAEh_PRE", "MAEh_POST", "MAEn_PRE", "MAEn_POST"))
for (gse in GSE_IDS) {
  sel <- cohorte == gse & !is.na(edad_real); n <- sum(sel)
  if (n < 3) { cat(sprintf("%-11s %6d  (sin edad suficiente)\n", gse, n)); next }
  cat(sprintf("%-11s %6d | %8.2f %8.2f | %8.2f %8.2f\n", gse, n,
      mean(abs(dnam_pre$Horvath[sel] - edad_real[sel])), mean(abs(dnam$Horvath[sel] - edad_real[sel])),
      mean(abs(dnam_pre$Hannum[sel]  - edad_real[sel])), mean(abs(dnam$Hannum[sel]  - edad_real[sel]))))
}

# Case-control coefficient, pre versus post regRCPqn (post reuses df_bias from step 5c).
df_pre <- data.frame(Hannum = dnam_pre$Hannum, Horvath = dnam_pre$Horvath,
                     grupo = factor(pheno$group), edad = edad_real)
df_pre <- df_pre[!is.na(df_pre$edad) & df_pre$grupo %in% c("SCZ", "CNT"), ]
df_pre$grupo <- droplevels(df_pre$grupo)
for (reloj in c("Hannum", "Horvath")) {
  b_pre  <- summary(lm(df_pre[[reloj]]  ~ grupo + edad, data = df_pre))$coefficients
  b_post <- summary(lm(df_bias[[reloj]] ~ grupo + edad, data = df_bias))$coefficients
  cat(sprintf("%s — β grupo PRE: %+.2f  | POST: %+.2f años\n", reloj,
              b_pre[grep("^grupo", rownames(b_pre)), 1],
              b_post[grep("^grupo", rownames(b_post)), 1]))
}
cat("Esperado: MAE y β grupo ~iguales pre vs post → regRCPqn no distorsiona la edad.\n")


# ---- 6. Save outputs ---------------------------------------------------------------------------

# betas_harmonized.rds: list(betas = harmonized betas, cohorte = cohort per sample); input of 05.
saveRDS(list(betas = betas_harm, cohorte = cohorte),
        file.path(DIR_HARMONIZED, "betas_harmonized.rds"))
cat("Guardado: betas_harmonized.rds —", nrow(betas_harm), "×", ncol(betas_harm), "(regRCPqn)\n")

# Covariate table, one row per sample, aligned by GSM ID.
covar <- data.frame(
  gsm_id          = colnames(betas_harm),
  gse             = cohorte,                                # cohort, fixed covariate in script 06
  group           = pheno$group,                            # SCZ or CNT
  sexo            = pheno$sexo,
  age_chronological = pheno$edad,                           # recorded age, with NA
  cell_frac,                                                # 7 EpiDISH cell fractions
  smoking_SSt     = smoke_out$PredictedSmokingStatus,       # EpiSmokEr status class
  smoking_SSc     = smoke_out[[sc_col]],                    # EpiSmokEr continuous score (model covariate)
  DNAmAge_Horvath = dnam$Horvath,                           # uncalibrated clock ages, post-regRCPqn
  DNAmAge_Hannum  = dnam$Hannum,
  DNAmAge_Horvath_pre = dnam_pre$Horvath,                   # pre-regRCPqn clock ages (step 5d)
  DNAmAge_Hannum_pre  = dnam_pre$Hannum,
  row.names = NULL, check.names = FALSE)

readr::write_csv(covar, file.path(DIR_PHENO, paste0("covariables_epigeneticas", SUF_VAR, SUF_COH, ".csv")))
cat("Guardado: covariables_epigeneticas", SUF_VAR, ".csv —", nrow(covar),
 "muestras ×", ncol(covar), "columnas\n")
print(utils::head(covar, 3))
