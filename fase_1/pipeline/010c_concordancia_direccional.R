# 010c — Directional concordance of each DMP across cohorts
#
# Describes between-cohort heterogeneity, which a one-stage model does not estimate per cohort: the
# script 06 model without the cohort term is fitted within each cohort (training set, single fit with
# the mean imputed age), and for each DMP it counts how many of the six cohorts have an adjusted
# coefficient with the sign of the pooled effect (Sokolov and Schiöth, 2024). The same count is given
# for the unadjusted within-cohort delta-beta. Run for the pre-registered DMP list and for the main
# analysis (smoking-associated CpGs removed, chip row/column as covariates, also within each cohort).
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds; data_combat_P2/phenotype_final.csv;
#         data_pheno/covariables_epigeneticas_P2.csv; data_pheno/phenodata_combinado_imputar_edad.csv
#         (chip position); data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv
# Output: data_dmp_P2_anc_concordancia/concordancia_por_DMP_<arm>.csv, concordancia_resumen.csv
#         (DMPs by number of concordant cohorts), efecto_ajustado_por_cohorte.csv (without chip position)
#         and efecto_ajustado_por_cohorte_pos.csv (with chip position)
# CON_POS: per arm, whether the within-cohort model adds chip row and column (TRUE for the main arm).

# ---- 1. Load data -----------------------------------------------------------------------------------
suppressMessages({ library(limma); library(data.table); library(readr) })
# Run from the repository root or from fase_1/pipeline
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

DIR_OUT <- "../data_dmp_P2_anc_concordancia"; dir.create(DIR_OUT, showWarnings = FALSE)
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")   # EpiDISH cell proportions
# Arms: pre-registered primary model (full panel) and main analysis (noSmk = no smoking CpGs + position)
ARMAS <- list(primario = "../data_dmp_P2_anc/dmp_resultados.csv",
              noSmk    = "../data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv")
CON_POS <- c(primario = FALSE, noSmk = TRUE)   # main arm keeps chip position within each cohort too

# Phenotype plus ancestry PCs and smoking score
pheno <- as.data.frame(read_csv("../data_combat_P2/phenotype_final.csv", show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                           # samples indexed by GSM accession
cov_ep <- as.data.frame(read_csv("../data_pheno/covariables_epigeneticas_P2.csv", show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)   # EPISTRUCTURE ancestry PCs
i <- match(rownames(pheno), cov_ep$gsm_id)                # align covariates to the phenotype rows
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]; pheno$smoking_SSc <- cov_ep$smoking_SSc[i]
# Chip position (row and column), as in script 06
pos_ph <- as.data.frame(read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv", show_col_types = FALSE))
pheno$position <- pos_ph$position[match(rownames(pheno), pos_ph$gsm_id)]
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))   # 450K slide: rows 1-6, columns 1-2
pheno$fila <- substr(pheno$position, 1, 3); pheno$columna <- substr(pheno$position, 4, 6)

# Main-analysis DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp written by script 06);
# a table without that column (pre-registered primary model) uses FDR < 0.05 alone
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
# DMPs of each arm; only their union is read from the beta matrix
ref <- list(); dmp_todas <- c()
for (a in names(ARMAS)) {
  r <- as.data.frame(fread(ARMAS[[a]])); ref[[a]] <- r[sel_dmp(r), c("cpg", "logFC_M", "delta_beta")]
  dmp_todas <- union(dmp_todas, ref[[a]]$cpg)
}
harm <- readRDS("../data_harmonized_P2/betas_harmonized.rds"); betas <- harm$betas[dmp_todas, ]; rm(harm); gc()
pheno <- pheno[colnames(betas), ]                         # phenotype in matrix column order
tr <- pheno$set == "train"; betas <- betas[, tr]; pheno <- pheno[tr, ]   # 80% training partition only
# M-values with betas clipped away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; bc <- pmin(pmax(betas, eps), 1 - eps); M <- log2(bc / (1 - bc))
pheno$group <- factor(pheno$group, levels = c("CNT", "SCZ")); pheno$sexo <- factor(pheno$sexo)
cohortes <- sort(unique(pheno$gse))
cat("DMP a evaluar (unión de brazos):", length(dmp_todas), "| muestras de entrenamiento:", ncol(M), "\n")

# ---- 2. Adjusted effect and unadjusted delta-beta within each cohort --------------------------------
coef_coh <- matrix(NA_real_, nrow(M), length(cohortes), dimnames = list(rownames(M), cohortes))
db_coh   <- coef_coh                                      # unadjusted within-cohort delta-beta
coef_coh_pos <- coef_coh                       # same model plus chip row and column (main arm)
for (g in cohortes) {
  s <- pheno$gse == g; p <- pheno[s, ]                    # samples of this cohort
  # Model of script 06 without the cohort term; single fit with the mean imputed age
  form <- paste("~ group + sexo + age_epigenetic_imputed +", paste(CELL_COVS, collapse = " + "),
                "+ smoking_SSc +", paste(ANC_COVS, collapse = " + "))
  if (length(unique(p$sexo)) < 2) form <- sub("sexo \\+ ", "", form)   # drop sex in single-sex cohorts
  di <- model.matrix(as.formula(form), data = p)
  fit <- lmFit(M[, s], di)                                # plain least squares; only the sign is used
  coef_coh[, g] <- fit$coefficients[, "groupSCZ"]         # adjusted case-control effect (M scale)
  form_pos <- form                             # row and column as factors, only if they vary in the cohort
  if (length(unique(p$fila)) > 1)    form_pos <- paste(form_pos, "+ factor(fila)")
  if (length(unique(p$columna)) > 1) form_pos <- paste(form_pos, "+ factor(columna)")
  fit_pos <- lmFit(M[, s], model.matrix(as.formula(form_pos), data = p))
  coef_coh_pos[, g] <- fit_pos$coefficients[, "groupSCZ"]   # adjusted effect with chip position
  # Unadjusted delta-beta (cases minus controls); may differ in sign from the adjusted effect
  db_coh[, g]   <- rowMeans(betas[, s][, p$group == "SCZ", drop = FALSE]) - rowMeans(betas[, s][, p$group == "CNT", drop = FALSE])
  cat(sprintf("  %s: n=%d (casos %d, controles %d) | columnas del diseño %d\n", g, sum(s),
              sum(p$group == "SCZ"), sum(p$group == "CNT"), ncol(di)))
}

# ---- 3. Concordance per DMP and summary -------------------------------------------------------------
# Reference sign = sign of the pooled logFC_M of the arm.
resumen <- data.frame()                                   # DMPs by number of concordant cohorts
for (a in names(ref)) {
  r <- ref[[a]]
  s_ref <- sign(r$logFC_M)                                # sign of the pooled effect
  cc <- if (CON_POS[[a]]) coef_coh_pos else coef_coh      # within-cohort model matching the arm
  conc_aj  <- rowSums(sign(cc[r$cpg, , drop = FALSE]) == s_ref)   # cohorts agreeing, adjusted effect
  conc_cru <- rowSums(sign(db_coh[r$cpg, , drop = FALSE]) == s_ref)   # cohorts agreeing, raw delta-beta
  tab <- data.frame(cpg = r$cpg, efecto_conjunto = r$logFC_M, concordancia_ajustada = conc_aj,
                    concordancia_cruda = conc_cru, cc[r$cpg, , drop = FALSE], check.names = FALSE)
  fwrite(tab, file.path(DIR_OUT, paste0("concordancia_por_DMP_", a, ".csv")))
  for (k in 6:0) {                                        # distribution from 6/6 down to 0/6
    resumen <- rbind(resumen, data.frame(brazo = a, cohortes_mismo_signo = k,
      DMP_ajustada = sum(conc_aj == k), pct_ajustada = round(100 * mean(conc_aj == k), 1),
      DMP_cruda = sum(conc_cru == k), pct_cruda = round(100 * mean(conc_cru == k), 1)))
  }
  cat(sprintf("\n[%s] %d DMP | ≥5/6 cohortes con el mismo signo: %.1f %% (ajustada), %.1f %% (cruda) | 6/6: %.1f %% | mediana %g\n",
      a, nrow(r), 100 * mean(conc_aj >= 5), 100 * mean(conc_cru >= 5), 100 * mean(conc_aj == 6), median(conc_aj)))
}
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "concordancia_resumen.csv"))
# Within-cohort adjusted effects of every evaluated DMP, without and with chip position
fwrite(data.frame(cpg = rownames(coef_coh), coef_coh, check.names = FALSE), file.path(DIR_OUT, "efecto_ajustado_por_cohorte.csv"))
fwrite(data.frame(cpg = rownames(coef_coh_pos), coef_coh_pos, check.names = FALSE), file.path(DIR_OUT, "efecto_ajustado_por_cohorte_pos.csv"))
cat("\nGuardado en", DIR_OUT, "\n")
