# 010h — Sensitivity with cases and controls matched on sex and age
#
# In several cohorts cases and controls differ in sex and age. The model adjusts for both, but a linear
# adjustment does not guarantee comparability when the distributions barely overlap. This script repeats
# the analysis in a matched subset: within each cohort, each case is matched 1:1 without replacement to a
# control of the same sex and close age (the age that enters the model, maximum difference CALIPER years).
# Smoking and neutrophils are not matched: they travel with the disease and are handled by adjustment and
# by the panel without smoking CpGs. Matching is greedy and written with for loops to keep it transparent.
# Balance is reported before and after (Fisher for sex, Welch t for age). The model of 06 is fitted in the
# matched subset and in random unmatched subsamples of the same size (power control), and compared with
# the reference DMPs of the pre-registered and main arms.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds, data_combat_P2/ (phenotype_final.csv, age_mids.rds),
#         data_pheno/covariables_epigeneticas_P2.csv, data_pheno/phenodata_combinado_imputar_edad.csv,
#         data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv, refs/smoking_cpgs_exclusion_450k.txt
# Output: data_dmp_P2_anc_sens_emparejado/balance_emparejamiento.csv, sens_emparejado_resumen.csv,
#         por_DMP_<arm>_emparejado.csv, dmp_<arm>_emparejado.csv (FDR < 0.05 only)
# CALIPER: maximum age difference (years) within a matched pair (default 5).
# N_AZAR: number of random unmatched subsamples of the matched size (default 3).

# ---- 1. Packages, paths and phenotype ---------------------------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice) })
set.seed(7777777)                                            # fixed seed: matching order is random
# Run from fase_1/pipeline/ whether started there or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

VARIANTE_MASKING <- "P2"                                     # probe-masking variant of the primary analysis
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)   # regRCPqn-harmonized betas (04)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)       # phenotype and age imputations (05)
DIR_PHENO  <- "../data_pheno"                                   # epigenetic covariates (04, 04b)
DIR_OUT    <- "../data_dmp_P2_anc_sens_emparejado"              # outputs of this script
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
CALIPER    <- 5            # maximum age difference (years) between matched case and control

CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")     # EpiDISH leukocyte proportions
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
tr <- pheno[pheno$set == "train", ]                          # matching uses the training set only
cat("Entrenamiento:", nrow(tr), "| casos:", sum(tr$group == "SCZ"), "| controles:", sum(tr$group == "CNT"),
    "\n")

# ---- 2. 1:1 matching by cohort, sex and age ---------------------------------------------------------
elegidos <- character(0)                                     # GSM of matched cases and controls
balance  <- data.frame()                                     # balance before and after, per cohort
for (g in sort(unique(tr$gse))) {
  s <- tr[tr$gse == g, ]
  for (sx in sort(unique(s$sexo))) {                         # exact matching on sex
    casos <- s[s$group == "SCZ" & s$sexo == sx, ]
    ctrl  <- s[s$group == "CNT" & s$sexo == sx, ]
    libres <- rep(TRUE, nrow(ctrl))                          # controls still available (no replacement)
    orden <- sample(seq_len(nrow(casos)))                    # cases visited in random order, fixed seed
    for (j in orden) {
      if (!any(libres)) break                                # no controls left in this cohort and sex
      d <- abs(ctrl$age_epigenetic_imputed - casos$age_epigenetic_imputed[j])   # age distance
      d[!libres] <- Inf                                      # used controls cannot be chosen again
      k <- which.min(d)                                      # nearest free control
      if (d[k] <= CALIPER) {                                 # keep the pair only within the caliper
        elegidos <- c(elegidos, casos$gsm_id[j], ctrl$gsm_id[k])
        libres[k] <- FALSE
      }
    }
  }
  m <- tr[tr$gsm_id %in% elegidos & tr$gse == g, ]           # matched subset of this cohort
  # Balance tests: Fisher for sex, Welch t for age, before and after matching
  balance <- rbind(balance, data.frame(
    gse = g, casos_antes = sum(s$group == "SCZ"), controles_antes = sum(s$group == "CNT"),
    pares = sum(m$group == "SCZ"),
    p_sexo_antes = fisher.test(table(s$sexo, s$group))$p.value,
    p_edad_antes = t.test(age_epigenetic_imputed ~ group, data = s)$p.value,
    p_sexo_despues = if (length(unique(m$sexo)) > 1) fisher.test(table(m$sexo, m$group))$p.value else 1,
    p_edad_despues = t.test(age_epigenetic_imputed ~ group, data = m)$p.value,
    dif_edad_despues = mean(m$age_epigenetic_imputed[m$group == "SCZ"]) -
                       mean(m$age_epigenetic_imputed[m$group == "CNT"])))
}
print(balance, row.names = FALSE, digits = 3)
fwrite(balance, file.path(DIR_OUT, "balance_emparejamiento.csv"))
cat("Participantes emparejados:", length(elegidos), "(", length(elegidos) / 2, "pares )\n")

# ---- 2b. Power control: random subsamples of the same size ------------------------------------------
# The matched subset is smaller, so fewer DMPs reach FDR < 0.05 even if matching changes nothing. N_AZAR
# random unmatched subsamples with the same number of cases and controls per cohort separate the effect of
# size from the effect of matching: a similar share of surviving DMPs means the loss is power.
N_AZAR <- 3                                                  # number of random subsamples
MUESTRAS <- list(emparejado = elegidos)                      # samples to fit: matched plus random
for (r in seq_len(N_AZAR)) {
  set.seed(2026 + r)                                         # one fixed seed per subsample
  ids_r <- character(0)
  for (g in balance$gse) {
    k <- balance$pares[balance$gse == g]                     # pairs obtained in this cohort
    # Same number of cases and of controls per cohort as the matched subset, drawn at random
    for (gr in c("SCZ", "CNT")) ids_r <- c(ids_r, sample(tr$gsm_id[tr$gse == g & tr$group == gr], k))
  }
  MUESTRAS[[paste0("azar_", r)]] <- ids_r
}
for (nm in names(MUESTRAS)) cat(nm, ":", length(MUESTRAS[[nm]]), "participantes\n")

# ---- 3. Model of 06 in each sample (matched and random), per arm ------------------------------------
mids <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))       # 20 age imputations (mice, norm method)
todos <- unique(unlist(MUESTRAS))                            # every sample used by any fit
b <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))$betas
b <- b[, todos, drop = FALSE]; gc()
# M-values with betas bounded away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; b <- pmin(pmax(b, eps), 1 - eps); M_all <- log2(b / (1 - b)); rm(b); gc()
# Same covariates as 06; chip row/column are added for the main arm below
formula_txt <- paste("~ group + gse + sexo + age_k +", paste(CELL_COVS, collapse = " + "), "+ smoking_SSc +",
                     paste(ANC_COVS, collapse = " + "))
resumen <- data.frame()                                      # one row per arm x sample
for (arma in names(ARMAS)) {
  ref <- as.data.frame(fread(ARMAS[[arma]]$ref))             # reference results of the arm
  ref_dmp <- ref[sel_dmp(ref), c("cpg", "logFC_M")]          # reference DMPs and adjusted effects
  filas <- rownames(M_all)
  if (ARMAS[[arma]]$quitar_tabaco) filas <- setdiff(filas, cpg_tabaco)   # main arm: no smoking CpGs
  stopifnot(setequal(filas, ref$cpg))                        # same panel as the reference arm
  for (nm in names(MUESTRAS)) {
    ids <- MUESTRAS[[nm]]
    M2 <- M_all[filas, ids, drop = FALSE]
    ph <- pheno[ids, ]
    ph$group <- factor(ph$group, levels = c("CNT", "SCZ")); ph$sexo <- factor(ph$sexo)   # CNT as reference
    n_imp <- mids$m                                          # number of age imputations
    coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
    for (k in seq_len(n_imp)) {                              # one age imputation per iteration
      ph$age_k <- complete(mids, k)[rownames(ph), "age_chronological"]
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
    # Comparison on the reference DMPs
    m <- merge(data.frame(cpg = ref_dmp$cpg, coef_ref = ref_dmp$logFC_M), d, by = "cpg")
    mismo_signo <- sign(m$coef) == sign(m$coef_ref)          # same direction as the reference
    pendiente   <- unname(coef(lm(m$coef ~ 0 + m$coef_ref))[1])   # slope through the origin (1 = same size)
    if (nm == "emparejado") {                                # per-CpG tables only for the matched fit
      fwrite(d[d$FDR < 0.05, ], file.path(DIR_OUT, paste0("dmp_", arma, "_emparejado.csv")))
      fwrite(m, file.path(DIR_OUT, paste0("por_DMP_", arma, "_emparejado.csv")))
    }
    resumen <- rbind(resumen, data.frame(
      brazo = arma, muestra = nm, n = ncol(M2), DMP_escenario = sum(fdr < 0.05), lambda = round(lambda, 3),
      DMP_ref = nrow(m), signo_igual_pct = round(100 * mean(mismo_signo), 1), pendiente = round(pendiente, 3),
      cor_efectos = round(cor(m$coef, m$coef_ref), 3), DMP_ref_siguen_FDR05 = sum(m$FDR < 0.05),
      pct_DMP_ref_siguen = round(100 * mean(m$FDR < 0.05), 1)))
    message(sprintf("  [%s | %s] n=%d DMP=%d λ=%.3f signo=%.1f%% pendiente=%.3f cor=%.3f siguen=%d/%d",
            arma, nm, ncol(M2), sum(fdr < 0.05), lambda, 100 * mean(mismo_signo), pendiente,
            cor(m$coef, m$coef_ref), sum(m$FDR < 0.05), nrow(m)))
    rm(M2, coef_mat, se_mat, d, m); gc()                     # free memory before the next fit
  }
}

# ---- 4. Summary and output --------------------------------------------------------------------------
# The question is answered by whether reference DMPs keep sign and size (slope near 1) once cases and
# controls no longer differ in sex or age. The share still at FDR < 0.05 is read against the random
# subsamples of the same size: if similar, the loss is power, not matching.
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "sens_emparejado_resumen.csv"))
for (arma in names(ARMAS)) {
  z <- resumen[resumen$brazo == arma & grepl("^azar_", resumen$muestra), "pct_DMP_ref_siguen"]   # random fits
  cat(sprintf("%s: emparejado %.1f %% · al azar %.1f %% (rango %.1f a %.1f)\n", arma,
      resumen$pct_DMP_ref_siguen[resumen$brazo == arma & resumen$muestra == "emparejado"],
      mean(z), min(z), max(z)))
}
cat("\nGuardado en", DIR_OUT, "\n")
