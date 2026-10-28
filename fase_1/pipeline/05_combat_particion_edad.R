# 05 — ComBat (quality control and sensitivity), age imputation and 80/20 partition
#
# ComBat betas are a quality-control and sensitivity object only; the primary analysis (06) uses the
# regRCPqn betas with cohort as a covariate, because two-step batch correction inflates false
# positives with unbalanced designs (Nygaard et al., 2016). Missing age is imputed with mice from the
# clocks ("norm", m = 20, primary; "pmm" as sensitivity), then a stratified 20% hold-out is assigned
# and exported. Age is imputed before the split: the clocks are pre-trained and the hold-out is only
# kept out of the differential methylation inference.
# Paper section: Covariates; Data Partition
# Input:  ../data_harmonized_P2/betas_harmonized.rds, ../data_pheno/covariables_epigeneticas_P2.csv
# Output: ../data_combat_P2/ (phenotype, age mids, ComBat betas), ../data_ml_P2/, ../data_for_thesis/
# VARIANTE_MASKING: "P1" | "P2" | "P3" (default "P2"); EXCLUIR_COHORTE: "" (default) or a GSE ID.
# FORCE_COMBAT: FALSE (reload the ComBat checkpoint if present, default) | TRUE (recompute).
# METODO_IMPUT: "norm" (primary, default) | "pmm"; the other method is fitted as the sensitivity set.

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
# Defined for consistency with script 04; not applied here (mice and the partition set seeds directly).
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
  "sva",                                                # ComBat
  "BiocParallel"
)
pqks_CRAN <- c(
  "data.table", "matrixStats",                          # fast CSV export; row variances
  "readr", "stringr", "ggplot2",
  "mice"                                                 # multiple imputation of age
)

InsLoad.pks(c(pqks_BC, pqks_CRAN))

# Stops early if a required package failed to load.
criticos <- c("sva", "mice")
faltan   <- criticos[!criticos %in% loadedNamespaces()]                  # critical but not loaded
if (length(faltan) > 0)
  stop("Paquetes críticos no cargados: ", paste(faltan, collapse = ", "),
       ". Reinicia el kernel y re-ejecuta esta celda.")
cat("Paquetes críticos cargados OK:", paste(criticos, collapse = ", "), "\n")


# ---- 2. Paths, parameters and data -------------------------------------------------------------

set.seed(7777777)                                       # global seed for the random steps below

# Runs from the repository root or from pipeline/.
if (!file.exists("../data_harmonized") && dir.exists("fase_1/pipeline")) {
  setwd("fase_1/pipeline")
}

GSE_IDS <- c("GSE84727", "GSE80417", "GSE147221",
             "GSE116378", "GSE116379", "GSE152027")

EXCLUIR_COHORTE <- ""                                   # "" keeps all six cohorts; must match 03 and 04
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""  # suffix only when a cohort is left out
GSE_IDS <- setdiff(GSE_IDS, EXCLUIR_COHORTE)            # cohorts iterated below

VARIANTE_MASKING <- "P2"                                # probe-masking variant; must match 02-04
SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)  # folder suffix

DIR_HARMONIZED <- paste0("../data_harmonized", SUF_VAR, SUF_COH)   # output of script 04
DIR_PHENO      <- "../data_pheno"
DIR_COMBAT     <- paste0("../data_combat", SUF_VAR, SUF_COH)       # output of this script

dir.create(DIR_COMBAT, showWarnings = FALSE, recursive = TRUE)

harm        <- readRDS(file.path(DIR_HARMONIZED, "betas_harmonized.rds"))   # list(betas, cohorte)
betas_harm  <- harm$betas                               # probes x samples, columns = GSM IDs
cohorte     <- harm$cohorte                             # cohort (GSE) of each sample
rm(harm); gc()                                          # frees the duplicate copy

# Covariates from script 04, rows reordered to the beta columns.
covar <- as.data.frame(readr::read_csv(file.path(DIR_PHENO,
                                              paste0("covariables_epigeneticas", SUF_VAR, SUF_COH, ".csv")),
                                       show_col_types = FALSE))
rownames(covar) <- covar$gsm_id                         # GSM ID as row key
covar <- covar[colnames(betas_harm), ]                  # same order as the beta columns

# Console checks: dimensions, alignment and group sizes.
cat("betas_harm:", nrow(betas_harm), "sondas ×", ncol(betas_harm), "muestras\n")
cat("¿covar alineado con las columnas de betas?: ", identical(rownames(covar), colnames(betas_harm)), "\n")
cat("¿cohorte alineado con las columnas?: ", length(cohorte) == ncol(betas_harm), "\n")
cat("Grupo:\n"); print(table(covar$group, useNA = "ifany"))


# ---- 3. ComBat (quality control and sensitivity only) ------------------------------------------
# Empirical Bayes correction of cohort mean and variance effects (Johnson 2007, sva) on M-values,
# which are closer to Gaussian than betas. Not used as input to the primary analysis.

# 3a. Checks that every cohort contains both diagnostic groups (batch not confounded with diagnosis).
tab_lote <- table(cohorte, covar$group)                 # samples per cohort x group
cat("Reparto cohorte × grupo:\n"); print(tab_lote)

grupos_por_cohorte <- rowSums(tab_lote > 0)             # groups present in each cohort
no_confundido <- all(grupos_por_cohorte >= 2)           # TRUE if every cohort has SCZ and CNT
cat("\n¿toda cohorte tiene ambos grupos (lote NO confundido con diagnóstico)?:", no_confundido, "\n")
if (!no_confundido)
  warning("Lote confundido con diagnóstico en alguna cohorte: el brazo ComBat-sensibilidad no podría ",
          "preservar el grupo ahí. El primario (cohorte como covariable en el 06) NO se ve afectado.")


# 3b. ComBat on M-values, with a checkpoint that is reloaded unless FORCE_COMBAT is TRUE.
# mod protects diagnosis, sex, age (Hannum clock, complete for all samples) and cell fractions;
# Eosino is omitted because the fractions sum to about 1. Cohort is the batch and is not in mod.
FORCE_COMBAT <- FALSE                                   # TRUE recomputes even if the checkpoint exists
ckpt_combat  <- file.path(DIR_COMBAT, "betas_combat_checkpoint.rds")

if (!FORCE_COMBAT && file.exists(ckpt_combat)) {
  betas_combat <- readRDS(ckpt_combat)                  # ComBat betas from an earlier run
  cat("ComBat RECARGADO del checkpoint:", ckpt_combat, "\n")
} else {
  eps <- 1e-6                                           # keeps betas away from 0 and 1 before the logit
  bc  <- pmin(pmax(betas_harm, eps), 1 - eps)           # clipped copy
  M   <- log2(bc / (1 - bc)); rm(bc); gc()              # M = log2(beta / (1 - beta))

  # Biological design that ComBat must not remove (no NA in any column).
  mod      <- model.matrix(~ group + sexo + DNAmAge_Hannum +
                             B + NK + CD4T + CD8T + Mono + Neutro,
                           data = covar)
  M_combat <- sva::ComBat(dat = M, batch = factor(cohorte),
                          mod = mod, par.prior = TRUE, prior.plots = FALSE)   # parametric priors
  rm(M); gc()

  betas_combat <- 2^M_combat / (2^M_combat + 1)         # back to the beta scale
  rm(M_combat); gc()
  saveRDS(betas_combat, ckpt_combat)
  cat("ComBat calculado y guardado en checkpoint:", ckpt_combat, "\n")
}

# Console check: size, NA, range and unchanged probe/sample names.
cat("betas_combat:", nrow(betas_combat), "×", ncol(betas_combat),
    "| NA:", sum(is.na(betas_combat)),
    "| rango:", paste(round(range(betas_combat), 4), collapse = " .. "), "\n")
cat("¿dimnames intactos?:", identical(dimnames(betas_combat), dimnames(betas_harm)), "\n")


# 3c. PCA on the 20,000 most variable probes, coloured by cohort, before (regRCPqn) and after ComBat.
# After ComBat the cohorts should overlap.
pal      <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
coh_cols <- setNames(pal[seq_along(GSE_IDS)], GSE_IDS)  # one colour per cohort

# PCA before ComBat (regRCPqn betas).
vars_a <- matrixStats::rowVars(betas_harm)              # variance per probe
top_a  <- order(vars_a, decreasing = TRUE)[1:20000]     # 20,000 most variable probes
pca_a  <- prcomp(t(betas_harm[top_a, ]), center = TRUE, scale. = FALSE)   # samples as rows
pct_a  <- round(100 * pca_a$sdev^2 / sum(pca_a$sdev^2), 1)              # % variance per component

# PCA after ComBat, same procedure.
vars_d <- matrixStats::rowVars(betas_combat)
top_d  <- order(vars_d, decreasing = TRUE)[1:20000]
pca_d  <- prcomp(t(betas_combat[top_d, ]), center = TRUE, scale. = FALSE)
pct_d  <- round(100 * pca_d$sdev^2 / sum(pca_d$sdev^2), 1)

options(repr.plot.width = 12, repr.plot.height = 5.5)
par(mfrow = c(1, 2))                                    # left = before, right = after
plot(pca_a$x[, 1], pca_a$x[, 2], col = coh_cols[cohorte], pch = 19, cex = 0.4,
     xlab = paste0("PC1 (", pct_a[1], "%)"), ylab = paste0("PC2 (", pct_a[2], "%)"),
     main = "Antes (regRCPqn) — color = cohorte")
legend("topright", names(coh_cols), col = coh_cols, pch = 19, bty = "n", cex = 0.7)
plot(pca_d$x[, 1], pca_d$x[, 2], col = coh_cols[cohorte], pch = 19, cex = 0.4,
     xlab = paste0("PC1 (", pct_d[1], "%)"), ylab = paste0("PC2 (", pct_d[2], "%)"),
     main = "Después (ComBat) — color = cohorte")
par(mfrow = c(1, 1))


# ---- 4. Age: calibration and multiple imputation (before the partition) ------------------------
# Recorded age is never overwritten; age_epigenetic_imputed holds the recorded age where it exists
# and the imputed age where it is missing.

# 4a. Calibration of recorded age on both clocks; R2 and MAE per cohort in samples with known age.
covar$age_chronological <- suppressWarnings(as.numeric(covar$age_chronological))   # numeric, NA kept
known <- !is.na(covar$age_chronological)                # samples with recorded age
cat("Edad conocida:", sum(known), "| a imputar (NA):", sum(!known), "\n")

fit_cal <- lm(age_chronological ~ DNAmAge_Hannum + DNAmAge_Horvath, data = covar[known, ])
cat("Calibración global: R²(coef. determinación) =", round(summary(fit_cal)$r.squared, 3), "\n\n")

# Per-cohort fit of the calibrated prediction: R2 and mean absolute error.
cat(sprintf("%-11s %6s %8s %8s\n", "cohorte", "n", "R2", "MAE"))
for (gse in GSE_IDS) {
  sel <- known & covar$gse == gse                       # known-age samples of this cohort
  if (sum(sel) < 5) { cat(sprintf("%-11s %6d  (n<5)\n", gse, sum(sel))); next }
  pr <- predict(fit_cal, covar[sel, ])                  # calibrated age
  cat(sprintf("%-11s %6d %8.3f %8.2f\n", gse, sum(sel),
              cor(pr, covar$age_chronological[sel])^2,
              mean(abs(pr - covar$age_chronological[sel]))))
}


# 4b. Multiple imputation with mice. Predictors: both clocks, sex, cohort, diagnosis and cell
# fractions; only age_chronological has missing values. "norm" (Bayesian linear regression) is the
# primary method and "pmm" with 10 donors is the sensitivity set (van Buuren; Morris et al., 2014).
# With a strong continuous predictor norm is more efficient, and pmm cannot extrapolate beyond the
# observed ages.
METODO_IMPUT <- "norm"                                    # primary method
DONORS_PMM   <- 10                                        # donors for the pmm sensitivity set
M_IMPUT  <- 20                                            # number of imputations
vars_imp <- c("age_chronological", "DNAmAge_Hannum", "DNAmAge_Horvath",
              "sexo", "gse", "group", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
imp_df <- covar[, vars_imp]                               # imputation data set
for (v in c("sexo", "gse", "group")) imp_df[[v]] <- factor(imp_df[[v]])   # categorical predictors

# Runs mice with the requested method and a fixed seed.
correr_mice <- function(metodo) {
  set.seed(7777777)
  if (metodo == "pmm")
    mice(imp_df, m = M_IMPUT, method = "pmm", donors = DONORS_PMM, printFlag = FALSE, seed = 7777777)
  else
    mice(imp_df, m = M_IMPUT, method = "norm", printFlag = FALSE, seed = 7777777)
}
mids      <- correr_mice(METODO_IMPUT)                    # primary set
METODO_SENS <- if (METODO_IMPUT == "norm") "pmm" else "norm"   # the other method
mids_sens <- correr_mice(METODO_SENS)                     # sensitivity set
cat("mice OK — PRIMARIO:", METODO_IMPUT, "| SENSIBILIDAD:", METODO_SENS,
    if (METODO_SENS=="pmm") paste0("(donors=",DONORS_PMM,")") else "", "| m =", mids$m, "\n")


# 4c. Point estimate per sample = mean of the m completed datasets (equals recorded age where known).
# The mids object is kept for pooling with Rubin's rules in script 06.
# n x m matrix of imputed ages; complete() keeps the row order of covar.
imp_mat <- sapply(seq_len(mids$m), function(i) complete(mids, i)$age_chronological)
covar$age_epigenetic_imputed <- rowMeans(imp_mat)        # mean over the m imputations

# Console checks: summary, recorded age preserved, no NA left.
cat("age_epigenetic_imputed — resumen:\n"); print(summary(covar$age_epigenetic_imputed))
cat("¿coincide con la edad real donde se conoce?:",
    isTRUE(all.equal(covar$age_epigenetic_imputed[known], covar$age_chronological[known])), "\n")
cat("NA restantes en la edad final:", sum(is.na(covar$age_epigenetic_imputed)), "(esperado 0)\n")


# ---- 5. Stratified 80/20 partition -------------------------------------------------------------
# 20% of samples within each cohort x diagnosis stratum are labelled "test" and excluded from the
# differential methylation inference.
PROP_TEST <- 0.20                                       # hold-out fraction
set.seed(7777777)                                       # reproducible partition
covar$set <- "train"                                    # everyone starts in train
estrato   <- interaction(covar$gse, covar$group, drop = TRUE)   # cohort x diagnosis strata
for (st in levels(estrato)) {
  idx    <- which(estrato == st)                        # samples in this stratum
  n_test <- round(PROP_TEST * length(idx))              # 20% of the stratum
  if (n_test >= 1) covar$set[sample(idx, n_test)] <- "test"   # random draw without replacement
}
# Console checks: set sizes, test share per cohort, set x group.
cat("Partición (muestras):\n"); print(table(covar$set))
cat("\n% test por cohorte:\n")
print(round(100 * prop.table(table(covar$gse, covar$set), 1)[, "test"], 1))
cat("\nReparto set × grupo:\n"); print(table(covar$set, covar$group))


# ---- 6. Save outputs ---------------------------------------------------------------------------
# phenotype_final.csv: covariates, age_epigenetic_imputed and set; used by script 06 with the
# regRCPqn betas. age_mids.rds: primary mids (norm). age_mids_sens.rds: sensitivity mids (pmm).
# betas_combat.rds: list(betas, cohorte), sensitivity object only.
# Phenotype and mids are written to both the ComBat folder and ../data_pheno.
readr::write_csv(covar, file.path(DIR_COMBAT, "phenotype_final.csv"))
cat("Guardado: phenotype_final.csv —", nrow(covar), "muestras ×", ncol(covar), "columnas\n")
readr::write_csv(covar, file.path(DIR_PHENO, "phenotype_final.csv"))
cat("Guardado: phenotype_final.csv —", nrow(covar), "muestras ×", ncol(covar), "columnas en", DIR_PHENO, "\n")

saveRDS(mids, file.path(DIR_COMBAT, "age_mids.rds"))     # primary imputations
cat("Guardado: age_mids.rds (m =", mids$m, ",", METODO_IMPUT, ")\n")
saveRDS(mids, file.path(DIR_PHENO, "age_mids.rds"))
cat("Guardado: age_mids.rds en", DIR_PHENO, "\n")
saveRDS(mids_sens, file.path(DIR_COMBAT, "age_mids_sens.rds"))   # sensitivity imputations
saveRDS(mids_sens, file.path(DIR_PHENO,  "age_mids_sens.rds"))
cat("Guardado: age_mids_sens.rds (", METODO_SENS, ") para la sensibilidad pmm↔norm del 06\n")

saveRDS(list(betas = betas_combat, cohorte = cohorte),
        file.path(DIR_COMBAT, "betas_combat.rds"))      # ComBat betas, sensitivity only
cat("Guardado: betas_combat.rds —", nrow(betas_combat), "×", ncol(betas_combat), "(sensibilidad)\n")
cat("Guardado: betas_combat.rds —", nrow(betas_combat), "×", ncol(betas_combat),
    "(sensibilidad) en", DIR_COMBAT, "\n")
cat("\nNOTA: insumo PRIMARIO del 06 = betas_harmonized.rds (regRCPqn, del 04) + phenotype_final.csv.\n")


# ---- 6.1 Export the 80% training set -----------------------------------------------------------
# regRCPqn betas (not ComBat) for the training samples, written as samples x CpGs with gsm_id as
# the first column, plus the matching phenotype. Gzip-compressed; the test set is not touched.
DIR_ML <- paste0("../data_ml", SUF_VAR, SUF_COH)        # machine-learning export folder
dir.create(DIR_ML, showWarnings = FALSE, recursive = TRUE)

# Reads the saved phenotype so this block also runs on its own.
ph_full <- as.data.frame(readr::read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"),
                                         show_col_types = FALSE))
ph_train <- ph_full[ph_full$set == "train", ]           # training samples only

if (!exists("betas_harm"))                              # reload the regRCPqn betas if not in memory
  betas_harm <- readRDS(file.path(DIR_HARMONIZED, "betas_harmonized.rds"))$betas
betas_train <- betas_harm[, ph_train$gsm_id, drop = FALSE]   # columns in phenotype order
cat("Train: ", ncol(betas_train), " muestras × ", nrow(betas_train), " CpGs | grupos: ",
    paste(names(table(ph_train$group)), table(ph_train$group), sep="=", collapse=" "), "\n", sep="")

readr::write_csv(ph_train, file.path(DIR_ML, "phenotype_train_80.csv"))   # labels and covariates
cat("Guardado: phenotype_train_80.csv (", nrow(ph_train), "×", ncol(ph_train), ")\n")

bt <- data.table::as.data.table(t(betas_train), keep.rownames = "gsm_id")   # samples as rows
data.table::fwrite(bt, file.path(DIR_ML, "betas_train_80.csv.gz"), compress = "gzip")
cat("Guardado: betas_train_80.csv.gz (", nrow(bt), "muestras ×", ncol(bt) - 1, "CpGs)\n")
cat("¿gsm_id de betas == fenotipo, mismo orden?: ",
    identical(bt$gsm_id, ph_train$gsm_id), "\n")         # betas and labels joined row by row
rm(bt); gc()


# ---- 6.2 Export the 20% test set ---------------------------------------------------------------
# Same format as 6.1 for the held-out samples, intended for final validation of a model fitted on
# the training set only.
DIR_ML <- paste0("../data_ml", SUF_VAR, SUF_COH)
dir.create(DIR_ML, showWarnings = FALSE, recursive = TRUE)

# Reads the saved phenotype so this block also runs on its own.
ph_full <- as.data.frame(readr::read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"),
                                         show_col_types = FALSE))
ph_test <- ph_full[ph_full$set == "test", ]             # held-out samples only

if (!exists("betas_harm"))                              # reload the regRCPqn betas if not in memory
  betas_harm <- readRDS(file.path(DIR_HARMONIZED, "betas_harmonized.rds"))$betas
betas_test <- betas_harm[, ph_test$gsm_id, drop = FALSE]     # columns in phenotype order
cat("Test: ", ncol(betas_test), " muestras × ", nrow(betas_test), " CpGs | grupos: ",
    paste(names(table(ph_test$group)), table(ph_test$group), sep="=", collapse=" "), "\n", sep="")

readr::write_csv(ph_test, file.path(DIR_ML, "phenotype_test_20.csv"))     # labels and covariates
cat("Guardado: phenotype_test_20.csv (", nrow(ph_test), "×", ncol(ph_test), ")\n")

bt <- data.table::as.data.table(t(betas_test), keep.rownames = "gsm_id")    # samples as rows
data.table::fwrite(bt, file.path(DIR_ML, "betas_test_20.csv.gz"), compress = "gzip")
cat("Guardado: betas_test_20.csv.gz (", nrow(bt), "muestras ×", ncol(bt) - 1, "CpGs)\n")
cat("¿gsm_id de betas == fenotipo, mismo orden?: ",
    identical(bt$gsm_id, ph_test$gsm_id), "\n")
rm(bt); gc()

# ---- 6.3 Uncompressed train/test copies --------------------------------------------------------
# Writes regRCPqn betas (samples x CpGs, GSM ID as row names) and phenotypes for the train and test
# sets as plain CSV files in ../data_for_thesis.
print("==> [7/7] Exportando datos para tesis (separados en train/test)...")

dir_tesis <- "../data_for_thesis"
dir.create(dir_tesis, showWarnings = FALSE, recursive = TRUE)

train_samples <- covar$gsm_id[covar$set == "train"]     # GSM IDs of the 80% set
test_samples  <- covar$gsm_id[covar$set == "test"]      # GSM IDs of the 20% set

# Training betas and phenotype.
print("Exportando datos de ENTRENAMIENTO...")
data.table::fwrite(as.data.frame(t(betas_harm[, train_samples])),
                   file = file.path(dir_tesis, "betas_train_for_ml.csv"), row.names = TRUE)
readr::write_csv(covar[covar$set == "train", ], file.path(dir_tesis, "pheno_train_for_ml.csv"))

# Test betas and phenotype.
print("Exportando datos de VALIDACIÓN...")
data.table::fwrite(as.data.frame(t(betas_harm[, test_samples])),
                   file = file.path(dir_tesis, "betas_test_for_ml.csv"), row.names = TRUE)
readr::write_csv(covar[covar$set == "test", ], file.path(dir_tesis, "pheno_test_for_ml.csv"))

print(paste("==> Exportación para tesis completada. Los 4 archivos se encuentran en:",
            normalizePath(dir_tesis)))
