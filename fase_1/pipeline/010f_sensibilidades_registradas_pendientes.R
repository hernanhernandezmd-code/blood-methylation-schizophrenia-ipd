# 010f — Registered covariate and batch-correction sensitivity analyses
#
# The registration pre-specified three comparisons against the primary model: with and without leukocyte
# proportions, with and without ancestry components, and before and after ComBat. This script refits the
# model of script 06 (one-stage limma, Rubin's rules over the 20 age imputations, bacon, BH) in each
# scenario, for the pre-registered arm (full panel) and the main arm (no smoking CpGs, chip row/column).
# On the reference DMPs it reports sign concordance, slope through the origin and correlation of adjusted
# effects; it also reports each scenario's DMP count and lambda, because removing a covariate can change
# the count through confounding rather than power.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds, data_combat_P2/ (phenotype_final.csv, age_mids.rds,
#         betas_combat.rds), data_pheno/covariables_epigeneticas_P2.csv,
#         data_pheno/phenodata_combinado_imputar_edad.csv (chip position),
#         data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv, refs/smoking_cpgs_exclusion_450k.txt
# Output: data_dmp_P2_anc_sens_registradas/sens_registradas_resumen.csv, por_DMP_<arm>_<scenario>.csv,
#         dmp_<arm>_<scenario>.csv (FDR < 0.05 only)
# ESCENARIOS: sin_celulas = drop the six cell proportions; sin_ancestria = drop the four ancestry PCs;
#   combat = full model on the ComBat-corrected betas of 05 (cohort kept as covariate).

# ---- 1. Packages, paths and data (loaded once) -----------------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice) })
set.seed(7777777)                                            # fixed seed for reproducibility
# Run from fase_1/pipeline/ whether started there or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

VARIANTE_MASKING <- "P2"                                     # probe-masking variant of the primary analysis
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)   # regRCPqn-harmonized betas (04)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)       # phenotype, imputations, ComBat betas (05)
DIR_PHENO  <- "../data_pheno"                                   # epigenetic covariates (04, 04b)
DIR_OUT    <- "../data_dmp_P2_anc_sens_registradas"             # outputs of this script
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

CELL_COVS  <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")    # EpiDISH leukocyte proportions
ESCENARIOS <- c("sin_celulas", "sin_ancestria", "combat")       # registered sensitivity scenarios
# Arms: reference DMP table, whether smoking CpGs leave the panel, whether chip row/column enter the model.
# primario = pre-registered model on the full panel; noSmk = main analysis (no smoking CpGs, chip position).
ARMAS <- list(primario = list(ref = "../data_dmp_P2_anc/dmp_resultados.csv",
                              quitar_tabaco = FALSE, pos = FALSE),
              noSmk    = list(ref = "../data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv",
                              quitar_tabaco = TRUE, pos = TRUE))
cpg_tabaco <- readLines("../refs/smoking_cpgs_exclusion_450k.txt")   # Joehanes 2016 + Bollepalli 2019 union
# DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp of 06); tables without it: FDR < 0.05
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05

pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                              # samples indexed by GSM
cov_ep <- as.data.frame(read_csv(file.path(DIR_PHENO,
           paste0("covariables_epigeneticas_", VARIANTE_MASKING, ".csv")), show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)   # EPISTRUCTURE ancestry PCs (04b)
i <- match(rownames(pheno), cov_ep$gsm_id)                   # align covariates to the phenotype
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[i]                   # continuous EpiSmokEr score, as in 06

# Chip position (row R01-R06, column C01-C02): technical covariate of the main analysis
pos_ph <- as.data.frame(read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv",
                                  show_col_types = FALSE))
pheno$position <- pos_ph$position[match(rownames(pheno), pos_ph$gsm_id)]
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))   # every sample has a valid position
pheno$fila <- factor(substr(pheno$position, 1, 3)); pheno$columna <- factor(substr(pheno$position, 4, 6))
mids <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))       # 20 age imputations (mice, norm method)
cat("Fenotipo OK — ancestría:", paste(ANC_COVS, collapse = ", "), "| imputaciones:", mids$m, "\n")

# ---- 2. Refit per scenario and arm ------------------------------------------------------------------
# Betas are loaded once per scenario: harmonized for sin_celulas and sin_ancestria, ComBat-corrected for
# combat. The model is that of 06; only the formula or the input matrix changes.
resumen <- data.frame()                                      # one row per arm x scenario
por_dmp <- list()                                            # per-DMP effects, reference vs scenario
for (esc in ESCENARIOS) {
  # Scenario betas, training set only
  if (esc == "combat") {
    b <- readRDS(file.path(DIR_COMBAT, "betas_combat.rds")); if (is.list(b)) b <- b$betas
  } else {
    b <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))$betas
  }
  ph <- pheno[colnames(b), ]                                 # phenotype in matrix column order
  train <- ph$set == "train"                                 # 80% training partition (05)
  b <- b[, train, drop = FALSE]; ph <- ph[train, ]; gc()
  # M-values with betas bounded away from 0 and 1: limma assumes roughly normal residuals
  eps <- 1e-6; b <- pmin(pmax(b, eps), 1 - eps); M_tr <- log2(b / (1 - b)); rm(b); gc()
  ph$group <- factor(ph$group, levels = c("CNT", "SCZ"))     # controls as reference level
  ph$sexo  <- factor(ph$sexo)

  # Scenario formula: drop cells or ancestry PCs; cohort stays as a fixed covariate in all scenarios
  covs_cel <- if (esc == "sin_celulas") "" else paste("+", paste(CELL_COVS, collapse = " + "))
  covs_anc <- if (esc == "sin_ancestria") "" else paste("+", paste(ANC_COVS, collapse = " + "))
  formula_txt <- paste("~ group + gse + sexo + age_k", covs_cel, "+ smoking_SSc", covs_anc)
  cat(sprintf("\n[%s] M_tr %d × %d | modelo: %s\n", esc, nrow(M_tr), ncol(M_tr), formula_txt))

  for (arma in names(ARMAS)) {
    ref <- as.data.frame(fread(ARMAS[[arma]]$ref))           # reference results of the arm
    ref_dmp <- ref[sel_dmp(ref), c("cpg", "logFC_M")]        # reference DMPs and adjusted effects
    filas <- rownames(M_tr)
    if (ARMAS[[arma]]$quitar_tabaco) filas <- setdiff(filas, cpg_tabaco)   # main arm: no smoking CpGs
    stopifnot(setequal(filas, ref$cpg))                     # same panel as the reference arm
    M2 <- M_tr[filas, , drop = FALSE]

    n_imp <- mids$m                                          # number of age imputations
    coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
    for (k in seq_len(n_imp)) {                             # one age imputation per iteration
      ph$age_k <- complete(mids, k)[rownames(ph), "age_chronological"]
      # Main arm adds chip row and column as technical factors
      di  <- model.matrix(as.formula(paste(formula_txt, if (ARMAS[[arma]]$pos) "+ fila + columna" else "")),
                          data = ph)
      fit <- eBayes(lmFit(M2, di))                           # moderated linear model per CpG
      coef_mat[, k] <- fit$coefficients[, "groupSCZ"]        # diagnosis effect (M scale)
      se_mat[, k]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]   # moderated SE
    }
    coef <- rowMeans(coef_mat)                               # pooled effect over imputations
    SE   <- sqrt(rowMeans(se_mat^2) + (1 + 1/n_imp) * rowVars(coef_mat))   # Rubin: within + between variance
    bcn  <- bacon(effectsizes = coef, standarderrors = SE)   # empirical-null correction of bias and inflation
    coef_bc <- as.numeric(es(bcn)); p_bc <- as.numeric(pval(bcn)); se_bc <- as.numeric(se(bcn))
    lambda  <- median((coef_bc / se_bc)^2) / qchisq(0.5, 1)  # genomic inflation after bacon
    fdr     <- p.adjust(p_bc, method = "BH")                 # Benjamini-Hochberg FDR
    d <- data.frame(cpg = rownames(M2), coef = coef_bc, p = p_bc, FDR = fdr)
    fwrite(d[d$FDR < 0.05, ], file.path(DIR_OUT, paste0("dmp_", arma, "_", esc, ".csv")))

    # Comparison on the reference DMPs
    m <- merge(data.frame(cpg = ref_dmp$cpg, coef_ref = ref_dmp$logFC_M), d, by = "cpg")
    mismo_signo <- sign(m$coef) == sign(m$coef_ref)          # same direction as the reference
    pendiente   <- unname(coef(lm(m$coef ~ 0 + m$coef_ref))[1])   # slope through the origin (1 = same size)
    clave <- paste(arma, esc, sep = "_")
    por_dmp[[clave]] <- data.frame(cpg = m$cpg, coef_ref = m$coef_ref, coef = m$coef, FDR = m$FDR)
    resumen <- rbind(resumen, data.frame(
      brazo = arma, escenario = esc, DMP_escenario = sum(fdr < 0.05), lambda = round(lambda, 3),
      DMP_ref = nrow(m), signo_igual_pct = round(100 * mean(mismo_signo), 1),
      pendiente = round(pendiente, 3), cor_efectos = round(cor(m$coef, m$coef_ref), 3),
      DMP_ref_siguen_FDR05 = sum(m$FDR < 0.05),
      pct_DMP_ref_siguen = round(100 * mean(m$FDR < 0.05), 1)))
    message(sprintf("  [%s | %s] DMP=%d λ=%.3f signo=%.1f%% pendiente=%.3f cor=%.3f siguen=%d/%d",
            arma, esc, sum(fdr < 0.05), lambda, 100 * mean(mismo_signo), pendiente,
            cor(m$coef, m$coef_ref), sum(m$FDR < 0.05), nrow(m)))
    rm(M2, coef_mat, se_mat, d, m); gc()                     # free memory before the next arm
  }
  rm(M_tr); gc()
}

# ---- 3. Summary and output --------------------------------------------------------------------------
# Sign concordance well above 50% and a slope near 1 mean the diagnosis effect does not depend on that
# covariate or on the batch correction; a large change without cells would mean leukocyte composition
# explains part of the signal, which is why the primary model includes them.
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "sens_registradas_resumen.csv"))
for (clave in names(por_dmp)) fwrite(por_dmp[[clave]], file.path(DIR_OUT, paste0("por_DMP_", clave, ".csv")))
cat("\nGuardado en", DIR_OUT,
    ": sens_registradas_resumen.csv, por_DMP_<brazo>_<escenario>.csv, dmp_<brazo>_<escenario>.csv\n")
