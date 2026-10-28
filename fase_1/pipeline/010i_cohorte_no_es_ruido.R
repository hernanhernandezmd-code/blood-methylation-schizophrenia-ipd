# 010i — Do between-cohort differences turn the signal into noise?
#
# After harmonization, mean methylation of controls still differs between cohorts at almost every CpG
# (04c). This script tests with the data two reasons why that does not make the result noise, each with a
# second check. (1) Cohort differences are additive: a shifted cohort is shifted equally in its cases and
# its controls, so the cohort term of the model cancels the shift. The shift of each cohort against a
# reference is measured separately in controls and in cases over all CpGs; checks: a second reference
# cohort, and the cohort x diagnosis interaction expected by chance (diagnosis permuted within cohort).
# (2) Directional concordance is not produced by that noise: if DMPs came from between-cohort noise, their
# within-cohort direction would be as random as that of any CpG. The share of DMPs concordant in all six
# cohorts is compared with non-significant CpGs, using the within-cohort model of 010c; checks: the crude
# difference instead of the adjusted effect, and the pre-registered arm. Training set only.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds, data_combat_P2/phenotype_final.csv,
#         data_pheno/covariables_epigeneticas_P2.csv, data_pheno/phenodata_combinado_imputar_edad.csv,
#         data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv
# Output: data_dmp_P2_anc_ruido/aditividad_cohorte.csv, concordancia_frente_a_nulo.csv
# N_NULO: number of non-significant CpGs (FDR > 0.5) sampled as the concordance null (default 20000).
# N_PERM: within-cohort permutations of diagnosis for the chance interaction (default 5).

# ---- 1. Packages, phenotype and training M-values ---------------------------------------------------
suppressMessages({ library(limma); library(data.table); library(readr); library(matrixStats) })
set.seed(20261022)                                           # fixed seed: null CpGs and permutations
# Run from fase_1/pipeline/ whether started there or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")
DIR_OUT <- "../data_dmp_P2_anc_ruido"; dir.create(DIR_OUT, showWarnings = FALSE)
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")  # EpiDISH leukocyte proportions
N_NULO   <- 20000          # non-significant CpGs for the expected concordance
N_PERM   <- 5              # within-cohort permutations of diagnosis
# Arms: reference table and whether chip row/column enter the within-cohort model
ARMAS <- list(primario = list(ref = "../data_dmp_P2_anc/dmp_resultados.csv",           pos = FALSE),
              noSmk    = list(ref = "../data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv", pos = TRUE))   # main

# Phenotype with ancestry PCs, smoking score and chip position
ph <- as.data.frame(read_csv("../data_combat_P2/phenotype_final.csv", show_col_types = FALSE));
rownames(ph) <- ph$gsm_id
cv <- as.data.frame(read_csv("../data_pheno/covariables_epigeneticas_P2.csv", show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cv), value = TRUE); i <- match(rownames(ph), cv$gsm_id)
ph[, ANC_COVS] <- cv[i, ANC_COVS]; ph$smoking_SSc <- cv$smoking_SSc[i]
pos_ph <- as.data.frame(read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv",
                                  show_col_types = FALSE))
ph$position <- pos_ph$position[match(rownames(ph), pos_ph$gsm_id)]
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", ph$position)))      # every sample has a valid position
ph$fila <- substr(ph$position, 1, 3); ph$columna <- substr(ph$position, 4, 6)   # chip row and column
ph <- ph[ph$set == "train", ]                                # 80% training partition (05)
ph$group <- factor(ph$group, levels = c("CNT", "SCZ")); ph$sexo <- factor(ph$sexo)   # CNT as reference
b <- readRDS("../data_harmonized_P2/betas_harmonized.rds")$betas[, ph$gsm_id]; gc()
# M-values with betas bounded away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; b <- pmin(pmax(b, eps), 1 - eps); M <- log2(b / (1 - b)); rm(b); gc()
cohortes <- sort(unique(ph$gse))
cat("M:", nrow(M), "×", ncol(M), "| cohortes:", paste(cohortes, collapse = ", "), "\n")

# ---- 2. Reason 1: each cohort's shift is the same in cases and controls -----------------------------
# For each cohort and CpG: shift of its controls against the reference controls, and the same shift in
# cases. If additive, both agree (r near 1, slope near 1) and their difference (the cohort x diagnosis
# interaction) is small. The observed interaction is compared with the one obtained when diagnosis is
# permuted N_PERM times within each cohort (same numbers of cases and controls).
medias <- function(Mx, grp, g, k) rowMeans(Mx[, ph$gse == g & grp == k, drop = FALSE])   # per-CpG group mean
aditiv <- data.frame()                                       # one row per reference x cohort
for (ref in c("GSE84727", "GSE80417")) {                     # two reference cohorts as a check
  mc_ref <- medias(M, ph$group, ref, "CNT"); ms_ref <- medias(M, ph$group, ref, "SCZ")
  for (g in setdiff(cohortes, ref)) {
    dc <- medias(M, ph$group, g, "CNT") - mc_ref             # cohort shift in controls
    ds <- medias(M, ph$group, g, "SCZ") - ms_ref             # cohort shift in cases
    sd_perm <- numeric(N_PERM)
    for (k in seq_len(N_PERM)) {                             # diagnosis permuted within each cohort
      gp <- ph$group
      for (h in c(g, ref)) { s <- which(ph$gse == h); gp[s] <- sample(ph$group[s]) }
      dcp <- medias(M, gp, g, "CNT") - medias(M, gp, ref, "CNT")
      dsp <- medias(M, gp, g, "SCZ") - medias(M, gp, ref, "SCZ")
      sd_perm[k] <- sd(dsp - dcp)                            # interaction spread under the null
    }
    # r and slope near 1 = additive shift; interaccion_sobre_azar near 1 = interaction at chance level
    aditiv <- rbind(aditiv, data.frame(
      referencia = ref, cohorte = g, r = cor(dc, ds), pendiente = unname(coef(lm(ds ~ 0 + dc))[1]),
      sd_desplazamiento = sd(dc), sd_interaccion = sd(ds - dc), sd_interaccion_azar = mean(sd_perm),
      interaccion_sobre_azar = sd(ds - dc) / mean(sd_perm)))
  }
}
print(aditiv, row.names = FALSE, digits = 3)
fwrite(aditiv, file.path(DIR_OUT, "aditividad_cohorte.csv"))

# ---- 3. Reason 2: six-cohort concordance of DMPs against non-significant CpGs -----------------------
# Same within-cohort model as 010c (no cohort term; main arm with chip row and column). For each CpG,
# count the cohorts where the adjusted effect, and the crude difference, carry the sign of the pooled
# effect. Non-significant CpGs = N_NULO drawn at random with FDR > 0.5 in the same arm.
# DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp of 06); tables without it: FDR < 0.05
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
conc_res <- data.frame()                                     # one row per arm x concordance type
for (arma in names(ARMAS)) {
  ref <- as.data.frame(fread(ARMAS[[arma]]$ref))             # pooled results of the arm
  dmp  <- ref$cpg[sel_dmp(ref)]                              # DMPs of the arm
  nulo <- sample(ref$cpg[ref$FDR > 0.5], N_NULO)             # clearly null CpGs from the same panel
  sel  <- c(dmp, nulo)                                       # DMPs first, then null CpGs
  s_ref <- sign(ref$logFC_M[match(sel, ref$cpg)])            # sign of the pooled effect
  aj <- matrix(NA_real_, length(sel), length(cohortes), dimnames = list(sel, cohortes)); cr <- aj
  for (g in cohortes) {
    s <- ph$gse == g; p <- ph[s, ]
    form <- paste("~ group + sexo + age_epigenetic_imputed +", paste(CELL_COVS, collapse = " + "),
                  "+ smoking_SSc +", paste(ANC_COVS, collapse = " + "))
    if (length(unique(p$sexo)) < 2) form <- sub("sexo \\+ ", "", form)   # drop sex if constant in cohort
    # Chip row and column only when they vary within the cohort
    if (ARMAS[[arma]]$pos && length(unique(p$fila)) > 1)    form <- paste(form, "+ factor(fila)")
    if (ARMAS[[arma]]$pos && length(unique(p$columna)) > 1) form <- paste(form, "+ factor(columna)")
    fit <- lmFit(M[sel, s], model.matrix(as.formula(form), data = p))   # within-cohort linear model
    aj[, g] <- fit$coefficients[, "groupSCZ"]                # adjusted within-cohort effect
    # Crude case-control difference of means in M-values
    cr[, g] <- rowMeans(M[sel, s][, p$group == "SCZ", drop = FALSE]) -
               rowMeans(M[sel, s][, p$group == "CNT", drop = FALSE])
  }
  for (tipo in c("ajustada", "cruda")) {
    mm <- if (tipo == "ajustada") aj else cr
    k  <- rowSums(sign(mm) == s_ref)                         # cohorts agreeing with the pooled sign
    es_dmp <- seq_along(sel) <= length(dmp)                  # TRUE for DMPs, FALSE for null CpGs
    a6 <- sum(k[es_dmp] == 6); n6 <- sum(k[!es_dmp] == 6)    # concordant in all six cohorts
    ft <- fisher.test(matrix(c(a6, length(dmp) - a6, n6, N_NULO - n6), 2))   # DMPs vs null CpGs
    # azar_puro_pct6 = 100/64: six-of-six agreement under independent random signs (2^-6)
    conc_res <- rbind(conc_res, data.frame(
      brazo = arma, concordancia = tipo, n_dmp = length(dmp), pct6_dmp = 100 * a6 / length(dmp),
      pct56_dmp = 100 * mean(k[es_dmp] >= 5), n_nulo = N_NULO, pct6_nulo = 100 * n6 / N_NULO,
      pct56_nulo = 100 * mean(k[!es_dmp] >= 5), veces6 = (a6 / length(dmp)) / (n6 / N_NULO),
      or_fisher = unname(ft$estimate), p_fisher = ft$p.value, azar_puro_pct6 = 100 / 64))
  }
  cat(sprintf("[%s] listo: %d DMP, %d CpG nulos\n", arma, length(dmp), N_NULO))
}
print(conc_res, row.names = FALSE, digits = 3)
fwrite(conc_res, file.path(DIR_OUT, "concordancia_frente_a_nulo.csv"))

# ---- 4. Reading -------------------------------------------------------------------------------------
# Reason 1 holds if r is high, the slope is near 1 and interaccion_sobre_azar is near 1: the cohort shift
# is common to cases and controls and the cohort term absorbs it. Reason 2 holds if DMPs concordant in six
# cohorts clearly exceed non-significant CpGs, with adjusted and crude effects, in both arms.
cat("Guardado en", DIR_OUT, "\n")
