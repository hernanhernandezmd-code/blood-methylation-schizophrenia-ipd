# 010j — Direction of the DMPs that fall on smoking-associated CpGs
#
# In the pre-registered primary model, smoking-associated CpGs are over-represented among the DMPs even
# after adjustment for the EpiSmokEr score. If that excess were residual confounding by smoking, those
# DMPs should change in the same direction as smoking does (cases smoke more). The script checks the sign
# of each case-control difference against two references: the current-versus-never smoking effect of
# Joehanes et al. 2016 (Supplemental Table 2), by strata of smoking-effect strength, and the association
# of each CpG with the smoking score in the training controls (crude and fully adjusted models). Positive
# controls: the 200 CpGs most associated with the score in the controls, and cg05575921 (AHRR).
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_dmp_P2_anc/dmp_resultados.csv; refs/smoking_cpgs_exclusion_450k.txt;
#         refs/joehanes2016_supplemental_tables.xlsx (not redistributed in this repository; download it
#         as described in fase_1/refs/joehanes2016_smoking_cpgs_450k.SOURCE.txt);
#         data_combat_P2/phenotype_final.csv; data_pheno/covariables_epigeneticas_P2.csv;
#         data_harmonized_P2/betas_harmonized.rds;
#         data_dmp_P2_anc_concordancia/concordancia_por_DMP_primario.csv (from 010c)
# Output: data_dmp_P2_anc_tabaco_direccion/ (direccion_frente_a_joehanes.csv,
#         direccion_frente_a_controles.csv, tabaco_direccion_resumen.csv)
# N_TOP: number of CpGs most associated with the smoking score used as positive control (default 200)

# ---- 1. Packages, primary DMPs, smoking panel and Joehanes effects -------------------------------
suppressMessages({ library(limma); library(data.table); library(readr); library(readxl) })
# works whether the session starts in fase_1/pipeline or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")
DIR_OUT <- "../data_dmp_P2_anc_tabaco_direccion"; dir.create(DIR_OUT, showWarnings = FALSE)   # output folder
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")    # EpiDISH leukocyte fractions
N_TOP <- 200     # positive control: CpGs most associated with the score in our controls

# Primary model results: raw and adjusted delta-beta give the direction of each case-control difference
prim <- as.data.frame(fread("../data_dmp_P2_anc/dmp_resultados.csv",
                            select = c("cpg", "FDR", "delta_beta", "delta_beta_adj", "logFC_M")))
panel <- readLines("../refs/smoking_cpgs_exclusion_450k.txt")  # Joehanes 2016 union Bollepalli 2019
prim <- prim[prim$cpg %in% panel, ]                      # smoking-panel CpGs only (16,922 in our panel)
prim$dmp <- prim$FDR < 0.05                              # DMP of the pre-registered primary model
cat("CpG del panel de tabaco evaluados:", nrow(prim), "| DMP primarias entre ellos:", sum(prim$dmp), "\n")

# Joehanes 2016, sheet 02 (current vs never, FDR < 0.05): probe, effect (column 7) and Z (column 9)
j <- suppressWarnings(as.data.frame(read_excel("../refs/joehanes2016_supplemental_tables.xlsx",
                                               sheet = "02", skip = 3)))
j <- data.frame(cpg = j[[1]], efecto_joe = as.numeric(j[[7]]), z_joe = as.numeric(j[[9]]))   # by position
j <- j[grepl("^cg", j$cpg) & !is.na(j$efecto_joe), ]   # CpG probes with a usable effect
j$rango_joe <- rank(-abs(j$z_joe))                       # 1 = CpG most strongly associated with smoking
prim <- merge(prim, j, by = "cpg", all.x = TRUE)         # keep panel CpGs absent from Joehanes as NA
cat("Con efecto de Joehanes:", sum(!is.na(prim$efecto_joe)), "\n")

# ---- 2. Direction against Joehanes, by strata of smoking-effect strength -------------------------
d <- prim[prim$dmp & !is.na(prim$efecto_joe), ]          # smoking-panel DMPs with a Joehanes effect
res_joe <- data.frame()                                  # one row per stratum
estratos <- list(c(0, 100), c(100, 1000), c(1000, 5000), c(5000, Inf), c(0, Inf))   # Joehanes rank bins
for (e in estratos) {
  s <- d[d$rango_joe > e[1] & d$rango_joe <= e[2], ]     # DMPs whose Joehanes rank falls in the bin
  if (nrow(s) == 0) next                                 # empty bins are skipped
  # share of DMPs whose case-control sign matches the smoking sign (raw and adjusted delta-beta)
  res_joe <- rbind(res_joe, data.frame(
    estrato_rango_joehanes = paste0(e[1], "-", ifelse(is.finite(e[2]), e[2], "fin")), n_dmp = nrow(s),
    pct_mismo_signo_crudo = 100 * mean(sign(s$delta_beta) == sign(s$efecto_joe)),
    pct_mismo_signo_ajustado = 100 * mean(sign(s$delta_beta_adj) == sign(s$efecto_joe))))
}
print(res_joe, row.names = FALSE, digits = 3)
fwrite(res_joe, file.path(DIR_OUT, "direccion_frente_a_joehanes.csv"))

# ---- 3. Direction against the smoking-score association in our own controls ----------------------
ph <- as.data.frame(read_csv("../data_combat_P2/phenotype_final.csv",          # phenotype and 80/20 split
                             show_col_types = FALSE)); rownames(ph) <- ph$gsm_id
# ancestry PCs and the continuous smoking score come from the epigenetic covariate table
cv <- as.data.frame(read_csv("../data_pheno/covariables_epigeneticas_P2.csv", show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cv), value = TRUE); i <- match(rownames(ph), cv$gsm_id)
ph[, ANC_COVS] <- cv[i, ANC_COVS]; ph$smoking_SSc <- cv$smoking_SSc[i]          # aligned by GSM
ph <- ph[ph$set == "train" & ph$group == "CNT", ]       # training controls only: no diagnosis effect
b <- readRDS("../data_harmonized_P2/betas_harmonized.rds")$betas[prim$cpg, ph$gsm_id]; gc()  # panel x ctrl
# M-values with betas bounded away from 0 and 1
M <- log2(pmin(pmax(b, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(b, 1e-6), 1 - 1e-6))); rm(b); gc()
cat("Controles de entrenamiento:", ncol(M), "| CpG:", nrow(M), "\n")

res_ctl <- data.frame()                                  # one row per model
for (modelo in c("crudo", "ajustado")) {
  # crude: score + cohort; adjusted: also sex, age, leukocyte fractions and ancestry PCs
  f <- if (modelo == "crudo") "~ smoking_SSc + gse" else
       paste("~ smoking_SSc + gse + sexo + age_epigenetic_imputed +", paste(CELL_COVS, collapse = " + "), "+",
             paste(ANC_COVS, collapse = " + "))
  fit <- eBayes(lmFit(M, model.matrix(as.formula(f), data = ph)))    # moderated linear model per CpG
  dir_smk <- sign(fit$coefficients[, "smoking_SSc"])    # direction of smoking in our controls
  t_smk <- fit$t[, "smoking_SSc"]                        # strength of the score association
  es_dmp <- prim$dmp                                     # primary DMPs within the smoking panel
  top <- order(-abs(t_smk))[seq_len(N_TOP)]             # positive control: most score-associated CpGs
  ahrr <- which(prim$cpg == "cg05575921")                # AHRR, the best replicated smoking CpG
  # sign agreement for DMPs (raw and adjusted), the whole panel, the positive control and AHRR
  res_ctl <- rbind(res_ctl, data.frame(
    modelo = modelo,
    pct_dmp_mismo_signo = 100 * mean(sign(prim$delta_beta[es_dmp]) == dir_smk[es_dmp]),
    pct_dmp_mismo_signo_adj = 100 * mean(sign(prim$delta_beta_adj[es_dmp]) == dir_smk[es_dmp]),
    pct_todo_panel_mismo_signo = 100 * mean(sign(prim$delta_beta) == dir_smk),
    control_positivo_top200_pct = 100 * mean(sign(prim$delta_beta[top]) == dir_smk[top]),
    ahrr_mismo_signo = sign(prim$delta_beta[ahrr]) == dir_smk[ahrr]))
}
print(res_ctl, row.names = FALSE, digits = 3)
fwrite(res_ctl, file.path(DIR_OUT, "direccion_frente_a_controles.csv"))

# AHRR across the six cohorts (from 010c): does it follow smoking (lower in cases) in all of them?
conc <- as.data.frame(fread("../data_dmp_P2_anc_concordancia/concordancia_por_DMP_primario.csv"))
a <- conc[conc$cpg == "cg05575921", ]                    # number of cohorts with the pooled sign
cat("AHRR: concordancia ajustada", a$concordancia_ajustada, "de 6 | Δβ crudo",
    round(prim$delta_beta[prim$cpg == "cg05575921"], 3), "\n")

# ---- 4. Summary ----------------------------------------------------------------------------------
# If most smoking-panel DMPs do not follow the smoking direction while the positive controls do, the
# excess is not mainly smoking; residual confounding sits in the most smoking-sensitive CpGs (AHRR and
# the strongest effects), which the main analysis removes.
# one-row summary: Joehanes agreement overall and in the top 100, control agreement, positive controls
resumen <- data.frame(
  n_dmp_tabaco = sum(prim$dmp), n_con_joehanes = nrow(d),
  pct_mismo_signo_joehanes_crudo = res_joe$pct_mismo_signo_crudo[res_joe$estrato_rango_joehanes == "0-fin"],
  pct_mismo_signo_joehanes_top100 = res_joe$pct_mismo_signo_crudo[res_joe$estrato_rango_joehanes == "0-100"],
  n_top100 = res_joe$n_dmp[res_joe$estrato_rango_joehanes == "0-100"],
  pct_mismo_signo_controles_crudo = res_ctl$pct_dmp_mismo_signo[res_ctl$modelo == "crudo"],
  pct_mismo_signo_controles_ajustado = res_ctl$pct_dmp_mismo_signo[res_ctl$modelo == "ajustado"],
  control_positivo_crudo = res_ctl$control_positivo_top200_pct[res_ctl$modelo == "crudo"],
  control_positivo_ajustado = res_ctl$control_positivo_top200_pct[res_ctl$modelo == "ajustado"],
  ahrr_concordancia = a$concordancia_ajustada, ahrr_con_tabaco = all(res_ctl$ahrr_mismo_signo))
print(t(resumen))
fwrite(resumen, file.path(DIR_OUT, "tabaco_direccion_resumen.csv"))
cat("Guardado en", DIR_OUT, "\n")
