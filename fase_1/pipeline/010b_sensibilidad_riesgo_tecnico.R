# 010b — Technical-risk sensitivity analysis
#
# Refits the model without the cohorts judged at high technical risk of bias: GSE147221 (Dublin), and
# GSE147221 plus GSE152027 (IoPPN, high risk under the literal protocol criteria), for the pre-registered
# arm (full panel) and the main arm (smoking-associated CpGs removed, chip row/column as covariates).
# The cohorts are not excluded from the main analysis; this is only a check. Removing cohorts lowers n
# and hence the DMP count even without bias, so on the reference DMPs it reports sign concordance,
# slope through the origin and correlation of the adjusted effects.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/, data_combat_P2/ (phenotype, age_mids.rds), covariables_epigeneticas_P2.csv,
#         data_pheno/phenodata_combinado_imputar_edad.csv (chip position),
#         data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv, refs/smoking_cpgs_exclusion_450k.txt
# Output: data_dmp_P2_anc_riesgo_tecnico/riesgo_tecnico_resumen.csv, por_DMP_<arm>.csv,
#         dmp_<arm>_<scenario>.csv (FDR < 0.05 only)
# USAR_RUBIN TRUE|FALSE (default TRUE): Rubin's rules over the 20 age imputations, or a single fit.

# ---- 1. Packages, paths and data (loaded once) ------------------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice) })
set.seed(7777777)                                         # fixed seed for reproducibility
# Run from the repository root or from fase_1/pipeline
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

USAR_RUBIN <- TRUE                                        # TRUE = Rubin over 20 imputations, as in 06
VARIANTE_MASKING <- "P2"                                  # probe-masking variant
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)   # harmonized betas (script 04)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)       # phenotype, partition, age imputations
DIR_PHENO  <- "../data_pheno"                             # epigenetic covariates
DIR_OUT    <- "../data_dmp_P2_anc_riesgo_tecnico"         # results of this check
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")   # EpiDISH cell proportions
# Scenarios: cohorts removed in each refit
ESCENARIOS <- list(sin_Dublin       = c("GSE147221"),
                   sin_Dublin_IoPPN = c("GSE147221", "GSE152027"))
# Arms: reference DMP table, whether smoking-associated CpGs are removed, whether chip position enters.
# noSmk = main analysis (smoking CpGs removed, chip row/column as covariates)
ARMAS <- list(primario = list(ref = "../data_dmp_P2_anc/dmp_resultados.csv",       quitar_tabaco = FALSE, pos = FALSE),
              noSmk    = list(ref = "../data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv", quitar_tabaco = TRUE, pos = TRUE))
cpg_tabaco <- readLines("../refs/smoking_cpgs_exclusion_450k.txt")   # union of Joehanes and Bollepalli
# Main-analysis DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp written by script 06);
# a table without that column (pre-registered primary model) uses FDR < 0.05 alone
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05

# Phenotype plus ancestry PCs and smoking score
pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                           # samples indexed by GSM accession
cov_ep <- as.data.frame(read_csv(file.path(DIR_PHENO,
           paste0("covariables_epigeneticas_", VARIANTE_MASKING, ".csv")), show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)   # EPISTRUCTURE ancestry PCs
i <- match(rownames(pheno), cov_ep$gsm_id)                # align covariates to the phenotype rows
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[i]                # continuous EpiSmokEr score

# Chip position (row and column): technical covariates of the main analysis, as in script 06
pos_ph <- as.data.frame(read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv", show_col_types = FALSE))
pheno$position <- pos_ph$position[match(rownames(pheno), pos_ph$gsm_id)]
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))   # 450K slide: rows 1-6, columns 1-2
pheno$fila <- factor(substr(pheno$position, 1, 3)); pheno$columna <- factor(substr(pheno$position, 4, 6))

# Betas loaded once; every arm and scenario subsets the same matrix
harm  <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds")); betas <- harm$betas; rm(harm); gc()
pheno <- pheno[colnames(betas), ]                         # phenotype in matrix column order
mids  <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))   # mice object (norm method)

# Training set and M-values
train    <- pheno$set == "train"                          # 80% training partition only
betas_tr <- betas[, train, drop = FALSE]; pheno_tr <- pheno[train, ]; rm(betas); gc()
# M-values with betas clipped away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; bc <- pmin(pmax(betas_tr, eps), 1 - eps); M_tr <- log2(bc / (1 - bc)); rm(bc, betas_tr); gc()
pheno_tr$group <- factor(pheno_tr$group, levels = c("CNT", "SCZ"))   # controls as reference level
pheno_tr$sexo  <- factor(pheno_tr$sexo)
extra_anc <- paste("+", paste(ANC_COVS, collapse = " + "))   # ancestry terms for the formula
cat("Carga OK — M_tr:", nrow(M_tr), "×", ncol(M_tr), "| Rubin:", USAR_RUBIN, "\n")

# ---- 2. Refit per arm and scenario ------------------------------------------------------------------
# Same model as script 06: limma per imputation, Rubin pooling, bacon, BH.
resumen <- data.frame()                                   # one row per arm and scenario
for (arma in names(ARMAS)) {
  ref <- as.data.frame(fread(ARMAS[[arma]]$ref))          # full-cohort results of this arm
  ref_dmp <- ref[sel_dmp(ref), c("cpg", "logFC_M")]       # reference DMPs and their adjusted effects
  filas <- rownames(M_tr)
  if (ARMAS[[arma]]$quitar_tabaco) filas <- setdiff(filas, cpg_tabaco)   # main arm: drop smoking CpGs
  stopifnot(setequal(filas, ref$cpg))                     # same panel as the reference arm
  por_dmp <- data.frame(cpg = ref_dmp$cpg, coef_ref = ref_dmp$logFC_M)
  cat(sprintf("\n[%s] panel %d CpG | DMP de referencia %d\n", arma, length(filas), nrow(ref_dmp)))

  for (nm in names(ESCENARIOS)) {
    keep <- !(as.character(pheno_tr$gse) %in% ESCENARIOS[[nm]])   # samples of the remaining cohorts
    p2 <- pheno_tr[keep, ]; p2$gse <- droplevels(factor(p2$gse))   # no empty cohort column in the design
    M2 <- M_tr[filas, keep, drop = FALSE]
    n_imp <- if (USAR_RUBIN) mids$m else 1L               # 20 imputations or a single fit
    coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
    for (k in seq_len(n_imp)) {
      # Age of imputation k, or the point imputed age when Rubin is off
      p2$age_k <- if (USAR_RUBIN) complete(mids, k)[rownames(p2), "age_chronological"]
                  else p2$age_epigenetic_imputed
      # Design of script 06, plus chip row and column in the main arm
      di  <- model.matrix(as.formula(paste("~ group + gse + sexo + age_k +",
                          paste(CELL_COVS, collapse = " + "), "+ smoking_SSc", extra_anc,
                          if (ARMAS[[arma]]$pos) "+ fila + columna" else "")), data = p2)
      fit <- eBayes(lmFit(M2, di))                        # moderated t statistics
      coef_mat[, k] <- fit$coefficients[, "groupSCZ"]     # case-control effect
      se_mat[, k]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]   # moderated standard error
    }
    # Rubin's rules: within- plus between-imputation variance
    coef <- rowMeans(coef_mat)
    if (n_imp > 1) { SE <- sqrt(rowMeans(se_mat^2) + (1 + 1/n_imp) * rowVars(coef_mat)) } else { SE <- se_mat[, 1] }
    bcn <- bacon(effectsizes = coef, standarderrors = SE)   # bacon: corrects inflation and bias
    coef_bc <- as.numeric(es(bcn)); p_bc <- as.numeric(pval(bcn)); se_bc <- as.numeric(se(bcn))
    lambda  <- median((coef_bc / se_bc)^2) / qchisq(0.5, 1)   # genomic lambda after bacon
    fdr     <- p.adjust(p_bc, method = "BH")              # Benjamini-Hochberg FDR
    d <- data.frame(cpg = rownames(M2), coef = coef_bc, p = p_bc, FDR = fdr)
    fwrite(d[d$FDR < 0.05, ], file.path(DIR_OUT, paste0("dmp_", arma, "_", nm, ".csv")))   # FDR < 0.05

    # Comparison on the reference DMPs (binomial test of sign concordance against 50%)
    m <- merge(por_dmp[, c("cpg", "coef_ref")], d, by = "cpg")
    mismo_signo <- sign(m$coef) == sign(m$coef_ref)       # random (~50%) if removed cohorts drove it
    pendiente   <- unname(coef(lm(m$coef ~ 0 + m$coef_ref))[1])   # ~1 = effect size preserved
    por_dmp[[paste0("coef_", nm)]] <- d$coef[match(por_dmp$cpg, d$cpg)]
    por_dmp[[paste0("FDR_", nm)]]  <- d$FDR[match(por_dmp$cpg, d$cpg)]
    # Summary row; DMP_ref_siguen_FDR05 is descriptive, since it depends on power
    resumen <- rbind(resumen, data.frame(
      brazo = arma, escenario = nm, n_train = sum(keep), DMP_escenario = sum(fdr < 0.05), lambda = round(lambda, 3),
      DMP_ref = nrow(m), signo_igual_pct = round(100 * mean(mismo_signo), 1),
      p_binomial_signo = signif(binom.test(sum(mismo_signo), nrow(m), 0.5)$p.value, 3),
      pendiente = round(pendiente, 3), cor_efectos = round(cor(m$coef, m$coef_ref), 3),
      DMP_ref_siguen_FDR05 = sum(m$FDR < 0.05)))
    message(sprintf("  [%s | %s] n=%d DMP=%d λ=%.3f signo=%.1f%% pendiente=%.3f cor=%.3f",
            arma, nm, sum(keep), sum(fdr < 0.05), lambda, 100 * mean(mismo_signo), pendiente,
            cor(m$coef, m$coef_ref)))
    rm(M2, coef_mat, se_mat, d, m); gc()                  # free memory before the next scenario
  }
  fwrite(por_dmp, file.path(DIR_OUT, paste0("por_DMP_", arma, ".csv")))   # per-DMP effects and FDR
}

# ---- 3. Summary and output --------------------------------------------------------------------------
# Sign concordance well above 50% and slope near 1 mean the reference DMPs do not depend on the removed
# cohorts; a lower DMP count is expected from the smaller sample and is not, alone, evidence of bias.
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "riesgo_tecnico_resumen.csv"))
cat("\nGuardado en", DIR_OUT, ": riesgo_tecnico_resumen.csv, por_DMP_<brazo>.csv, dmp_<brazo>_<escenario>.csv\n")
