# 010k — Remaining registered sensitivity analyses: smoking class, epigenetic age acceleration, no regRCPqn
#
# Runs three comparisons announced in section 16 of the OSF registration (2u56r). (1) clase_tabaco: the
# model with the EpiSmokEr smoking class (smoking_SSt) in place of the continuous score (smoking_SSc);
# self-reported smoking is not available in any cohort. (2) aceleracion_edad: adds epigenetic age
# acceleration, the residual of a robust regression (MASS::rlm, as registered) of Hannum, Horvath and
# PhenoAge (Levine, computed here with methylclock) on chronological age. GrimAge is not implemented in
# methylclock and is not computed. This scenario uses only participants with recorded age, because in the
# others the imputed age comes from the same clocks and their acceleration would be about zero by
# construction; it is compared with the same subset without acceleration to isolate its effect.
# (3) sin_regRCPqn: the model on dasen + BMIQ betas without regRCPqn harmonization (output of script 03).
# Each scenario is run for the pre-registered primary model (full panel) and for the main analysis
# (smoking CpGs removed, chip row and column as covariates). Same method as script 06 (one-stage limma,
# Rubin's rules over 20 age imputations, bacon, BH); on the reference DMPs it reports sign agreement,
# slope through the origin, correlation of effects and how many keep FDR < 0.05.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds; data_normalized_P2/betas_normalized.rds;
#         data_combat_P2/ (phenotype_final.csv, age_mids.rds); data_pheno/covariables_epigeneticas_P2.csv;
#         data_pheno/phenodata_combinado_imputar_edad.csv (chip position);
#         data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv; refs/smoking_cpgs_exclusion_450k.txt
# Output: data_dmp_P2_anc_sens_restantes/ (sens_restantes_resumen.csv, dmp_<arm>_<scenario>.csv with
#         FDR < 0.05 only, por_DMP_<arm>_<scenario>_vs_<reference>.csv)
# VARIANTE_MASKING: probe mask variant of the input folders (P2 = registered mask)
# ARMAS: primario = pre-registered full panel; noSmk = main analysis (_noSmk_pos)

# ---- 1. Packages, paths and phenotype (loaded once) ----------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice); library(MASS); library(methylclock) })
set.seed(7777777)                                                   # reproducible bacon Gibbs sampler
# works whether the session starts in fase_1/pipeline or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

VARIANTE_MASKING <- "P2"                                              # registered probe mask
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)      # regRCPqn-harmonized betas (script 04)
DIR_NORM   <- paste0("../data_normalized_", VARIANTE_MASKING)     # dasen + BMIQ, before regRCPqn (script 03)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)        # phenotype, split and age imputations (05)
DIR_PHENO  <- "../data_pheno"
DIR_OUT    <- "../data_dmp_P2_anc_sens_restantes"                     # output folder of this script
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

CELL_COVS  <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")       # EpiDISH leukocyte fractions
ESCENARIOS <- c("clase_tabaco", "aceleracion_edad", "sin_regRCPqn")
# Arms: reference DMP table, whether smoking CpGs are removed, whether chip row/column enter the model
ARMAS <- list(primario = list(ref = "../data_dmp_P2_anc/dmp_resultados.csv",
                              quitar_tabaco = FALSE, pos = FALSE),          # pre-registered primary model
              noSmk    = list(ref = "../data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv",
                              quitar_tabaco = TRUE,  pos = TRUE))           # main analysis
cpg_tabaco <- readLines("../refs/smoking_cpgs_exclusion_450k.txt")  # Joehanes 2016 union Bollepalli 2019
# DMPs of script 06: es_dmp (FDR and effect threshold) where present, otherwise FDR < 0.05
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05

pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                                       # one row per sample, keyed by GSM
cov_ep <- as.data.frame(read_csv(file.path(DIR_PHENO,
                                           paste0("covariables_epigeneticas_", VARIANTE_MASKING, ".csv")),
                                 show_col_types = FALSE))                # cells, smoking, clocks, ancestry
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)    # EPISTRUCTURE ancestry PCs
i <- match(rownames(pheno), cov_ep$gsm_id)                           # align covariates by GSM
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[i]                         # EpiSmokEr score (registered covariate)
pheno$smoking_SSt <- factor(cov_ep$smoking_SSt[i])                 # EpiSmokEr class
pheno$DNAmAge_Hannum  <- cov_ep$DNAmAge_Hannum[i]                  # epigenetic ages from script 04
pheno$DNAmAge_Horvath <- cov_ep$DNAmAge_Horvath[i]
# chip position (RxxCxx) split into row and column factors for the main-analysis arm
pos_ph <- as.data.frame(read_csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                                 show_col_types = FALSE))
pheno$position <- pos_ph$position[match(rownames(pheno), pos_ph$gsm_id)]  # chip position of each sample
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))          # 450K: six rows by two columns
# row R01-R06, column C01-C02
pheno$fila <- factor(substr(pheno$position, 1, 3)); pheno$columna <- factor(substr(pheno$position, 4, 6))
mids <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))              # 20 multiple imputations of age (mice)
cat("Clase de tabaquismo:\n"); print(table(pheno$smoking_SSt[pheno$set == "train"],   # class by group, train
                                             pheno$group[pheno$set == "train"]))
cat("Entrenamiento con edad registrada:", sum(pheno$set == "train" & !is.na(pheno$age_chronological)), "\n")

# ---- 2. Refit per scenario and arm ---------------------------------------------------------------
# clase_tabaco and sin_regRCPqn use the 20 age imputations (Rubin), as script 06. aceleracion_edad uses
# participants with recorded age, so a single fit with real age, with and without the three
# accelerations computed in that subset. PhenoAge comes from methylclock on the harmonized betas of the
# subset; the clock CpGs missing from the panel (16 of 513) are imputed by KNN inside the package.
ajustar <- function(M2, ph, formula_txt, rubin) {       # returns bacon-corrected coef and p, and lambda
  if (rubin) {
    n_imp <- mids$m                                                   # number of imputations (20)
    # probes x imputations
    coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
    for (k in seq_len(n_imp)) {                          # one limma fit per imputed age
      ph$age_k <- complete(mids, k)[rownames(ph), "age_chronological"]  # imputed age for imputation k
      fit <- eBayes(lmFit(M2, model.matrix(as.formula(formula_txt), data = ph)))
      coef_mat[, k] <- fit$coefficients[, "groupSCZ"]                    # diagnosis effect in M-values
      se_mat[, k]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]   # moderated SE
    }
    # Rubin: total variance = mean within-imputation variance + inflated between-imputation variance
    coef <- rowMeans(coef_mat); SE <- sqrt(rowMeans(se_mat^2) + (1 + 1/n_imp) * rowVars(coef_mat))   # Rubin
  } else {
    ph$age_k <- ph$age_chronological                     # recorded age only, no imputation
    fit <- eBayes(lmFit(M2, model.matrix(as.formula(formula_txt), data = ph)))  # single fit with recorded age
    coef <- fit$coefficients[, "groupSCZ"]; SE <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]
  }
  bcn <- bacon(effectsizes = coef, standarderrors = SE)  # bacon: corrects bias and inflation
  # corrected effects, p and SE
  cb <- as.numeric(es(bcn)); pb <- as.numeric(pval(bcn)); sb <- as.numeric(se(bcn))
  # BH FDR and genomic inflation lambda after bacon
  data.frame(cpg = rownames(M2), coef = cb, p = pb, FDR = p.adjust(pb, method = "BH"),
             lambda = median((cb / sb)^2) / qchisq(0.5, 1))
}
comparar <- function(ref_dmp, d) {                      # on the reference DMPs
  m <- merge(data.frame(cpg = ref_dmp$cpg, coef_ref = ref_dmp$coef), d[, c("cpg", "coef", "FDR")], by = "cpg")
  # sign agreement, slope through the origin, correlation and DMPs still at FDR < 0.05
  list(m = m, signo = 100 * mean(sign(m$coef) == sign(m$coef_ref)),
       pend = unname(coef(lm(m$coef ~ 0 + m$coef_ref))[1]), cor = cor(m$coef, m$coef_ref),
       siguen = sum(m$FDR < 0.05))
}

resumen <- data.frame()  # one row per arm, scenario and comparison
for (esc in ESCENARIOS) {
  # sin_regRCPqn reads the normalized betas before harmonization; the others the harmonized ones
  b <- if (esc == "sin_regRCPqn") readRDS(file.path(DIR_NORM, "betas_normalized.rds"))$betas else
         readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))$betas
  ph <- pheno[colnames(b), ]                                          # phenotype in matrix column order
  sel <- ph$set == "train"                               # inference on the 80% training set only
  if (esc == "aceleracion_edad") sel <- sel & !is.na(ph$age_chronological)   # recorded age only
  b <- b[, sel, drop = FALSE]; ph <- ph[sel, ]; gc()                  # keep training (sub)set
  if (esc == "aceleracion_edad") {                      # PhenoAge (Levine) on the betas of the subset
    reloj <- as.data.frame(DNAmAge(b, clocks = "Levine", cell.count = FALSE))   # same call as script 04
    ph$DNAmAge_PhenoAge <- reloj$Levine[match(rownames(ph), reloj$id)]  # align clock output by GSM
    stopifnot(sum(is.na(ph$DNAmAge_PhenoAge)) == 0)                   # every subset sample gets a PhenoAge
    cat("PhenoAge: media", round(mean(ph$DNAmAge_PhenoAge), 1), "| r con la edad cronológica",
        round(cor(ph$DNAmAge_PhenoAge, ph$age_chronological), 3), "\n")
  }
  # M-values: limma assumes roughly normal residuals; betas bounded away from 0 and 1
  eps <- 1e-6; b <- pmin(pmax(b, eps), 1 - eps); M_tr <- log2(b / (1 - b)); rm(b); gc()
  ph$group <- factor(ph$group, levels = c("CNT", "SCZ")); ph$sexo <- factor(ph$sexo)   # controls as reference
  # acceleration = robust (rlm) residual of each clock on chronological age
  if (esc == "aceleracion_edad") {
    ph$acel_hannum   <- resid(rlm(DNAmAge_Hannum   ~ age_chronological, data = ph))
    ph$acel_horvath  <- resid(rlm(DNAmAge_Horvath  ~ age_chronological, data = ph))
    ph$acel_phenoage <- resid(rlm(DNAmAge_PhenoAge ~ age_chronological, data = ph))
    cat("Aceleración media SCZ − CNT (años): Hannum", round(diff(tapply(ph$acel_hannum, ph$group, mean)), 2),
        "| Horvath", round(diff(tapply(ph$acel_horvath, ph$group, mean)), 2),
        "| PhenoAge", round(diff(tapply(ph$acel_phenoage, ph$group, mean)), 2), "\n")
  }
  # covariates shared by every scenario: cohort as fixed covariate, sex, age, cells and ancestry
  base <- paste("~ group + gse + sexo + age_k +", paste(CELL_COVS, collapse = " + "), "+",
                paste(ANC_COVS, collapse = " + "))
  for (arma in names(ARMAS)) {
    ref <- as.data.frame(fread(ARMAS[[arma]]$ref))       # reference results of the arm
    # reference DMPs and their effects
    ref_dmp <- data.frame(cpg = ref$cpg[sel_dmp(ref)], coef = ref$logFC_M[sel_dmp(ref)])
    filas <- intersect(rownames(M_tr), ref$cpg)        # same panel as the reference arm
    if (ARMAS[[arma]]$quitar_tabaco) filas <- setdiff(filas, cpg_tabaco)   # main analysis: no smoking CpGs
    M2 <- M_tr[filas, , drop = FALSE]                                 # panel of the arm
    pos_txt <- if (ARMAS[[arma]]$pos) "+ fila + columna" else ""   # chip row and column in the main analysis
    if (esc == "clase_tabaco") {                         # smoking class instead of the score
      d <- ajustar(M2, ph, paste(base, "+ smoking_SSt", pos_txt), rubin = TRUE)
      comp <- list(referencia = ref_dmp)                              # compared only with its reference DMPs
    } else if (esc == "sin_regRCPqn") {                  # registered model, unharmonized betas
      d <- ajustar(M2, ph, paste(base, "+ smoking_SSc", pos_txt), rubin = TRUE)
      comp <- list(referencia = ref_dmp)                              # compared only with its reference DMPs
    } else {
      d0 <- ajustar(M2, ph, paste(base, "+ smoking_SSc", pos_txt), rubin = FALSE)   # subset, no acceleration
      d  <- ajustar(M2, ph,                                                       # subset, with acceleration
                    paste(base, "+ smoking_SSc + acel_hannum + acel_horvath + acel_phenoage", pos_txt),
                    rubin = FALSE)
      # compared with the full reference and with the same subset without acceleration
      comp <- list(referencia = ref_dmp,
                   subconjunto_sin_aceleracion = data.frame(cpg = d0$cpg[d0$cpg %in% ref_dmp$cpg],
                                                            coef = d0$coef[d0$cpg %in% ref_dmp$cpg]))
    }
    fwrite(d[d$FDR < 0.05, c("cpg", "coef", "p", "FDR")],                         # scenario DMPs only
           file.path(DIR_OUT, paste0("dmp_", arma, "_", esc, ".csv")))
    for (contra in names(comp)) {                        # one summary row per comparison
      r <- comparar(comp[[contra]], d)                                # sign, slope, correlation, survivors
      # per-DMP comparison table
      fwrite(r$m, file.path(DIR_OUT, paste0("por_DMP_", arma, "_", esc, "_vs_", contra, ".csv")))
      resumen <- rbind(resumen, data.frame(
        brazo = arma, escenario = esc, comparado_con = contra, n = ncol(M2), CpG = nrow(M2),
        DMP_escenario_FDR05 = sum(d$FDR < 0.05), lambda = round(d$lambda[1], 3), DMP_ref = nrow(r$m),
        signo_igual_pct = round(r$signo, 1), pendiente = round(r$pend, 3), cor_efectos = round(r$cor, 3),
        DMP_ref_siguen_FDR05 = r$siguen, pct_DMP_ref_siguen = round(100 * r$siguen / nrow(r$m), 1)))
      message(sprintf(
        "  [%s | %s | vs %s] n=%d DMP=%d λ=%.3f signo=%.1f%% pendiente=%.3f cor=%.3f siguen=%d/%d",
        arma, esc, contra, ncol(M2), sum(d$FDR < 0.05), d$lambda[1], r$signo, r$pend, r$cor,
        r$siguen, nrow(r$m)))
    }
    rm(M2); gc()                                                      # free memory before the next arm
  }
  rm(M_tr); gc()                                         # free the matrix before the next scenario
}

# ---- 3. Summary and output -----------------------------------------------------------------------
# Kept sign and a slope near 1 on the reference DMPs mean the diagnosis effect does not depend on how
# smoking is represented, on epigenetic age acceleration or on regRCPqn harmonization. For
# aceleracion_edad the comparison that isolates acceleration is against the same subset without it.
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "sens_restantes_resumen.csv"))
cat("\nGuardado en", DIR_OUT, "\n")
