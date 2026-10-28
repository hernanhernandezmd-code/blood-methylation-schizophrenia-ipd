# 010 — Leave-one-cohort-out analysis
#
# Refits the primary differential methylation model of script 06 six times, each time leaving out one
# cohort, and compares the number of DMPs, lambda, delta-beta and the adjusted effects of the primary
# DMPs with the full model. An additional seventh scenario leaves out GSE116378 and GSE152027 together
# (the two cohorts with the largest case-control imbalance); it is computed but not reported in the paper.
# Whole cohorts are removed, never individual subjects, so the check cannot be tuned sample by sample.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds; data_combat_P2/phenotype_final.csv and age_mids.rds;
#         data_pheno/covariables_epigeneticas_P2.csv; data_dmp_P2_anc/dmp_resultados.csv
# Output: data_dmp_P2_anc_LOCO/loco_resumen.csv, loco_por_DMP_primario.csv (coefficient and FDR of
#         each primary DMP per scenario), dmp_<scenario>.csv (FDR < 0.05 only)
# USAR_RUBIN TRUE|FALSE (default TRUE): Rubin's rules over the 20 age imputations, or one fit.

# ---- 1. Packages, paths and data (loaded once) ------------------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice)
                   library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
set.seed(7777777)                                         # fixed seed for reproducibility
# Run from the repository root or from fase_1/pipeline
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

USAR_RUBIN <- TRUE                                        # TRUE = same pooling as the primary model
VARIANTE_MASKING <- "P2"                                  # probe-masking variant of the primary model
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)   # harmonized betas (script 04)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)       # phenotype, partition, age imputations
DIR_PHENO  <- "../data_pheno"                             # epigenetic covariates
DIR_PRIM   <- "../data_dmp_P2_anc"                        # full six-cohort model
DIR_OUT    <- "../data_dmp_P2_anc_LOCO"                   # leave-one-cohort-out results
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")   # EpiDISH cell proportions
DB_HICONF <- 0.02                                         # |delta-beta| flag for high-confidence DMPs

# Phenotype plus ancestry PCs and smoking score
pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                           # samples indexed by GSM accession
cov_ep <- as.data.frame(read_csv(file.path(DIR_PHENO,
           paste0("covariables_epigeneticas_", VARIANTE_MASKING, ".csv")), show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)   # EPISTRUCTURE ancestry PCs
i <- match(rownames(pheno), cov_ep$gsm_id)                # align covariates to the phenotype rows
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[i]                # continuous EpiSmokEr score

# Betas loaded once; every scenario subsets the same matrix
harm  <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds")); betas <- harm$betas; rm(harm); gc()
pheno <- pheno[colnames(betas), ]                         # phenotype in matrix column order
mids  <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))   # mice object (norm method)

# Training set and M-values
train    <- pheno$set == "train"                          # 80% training partition only
betas_tr <- betas[, train, drop = FALSE]; pheno_tr <- pheno[train, ]; rm(betas); gc()
# M-values with betas clipped away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; bc <- pmin(pmax(betas_tr, eps), 1 - eps); M_tr <- log2(bc / (1 - bc)); rm(bc); gc()
pheno_tr$group <- factor(pheno_tr$group, levels = c("CNT", "SCZ"))   # controls as reference level
pheno_tr$sexo  <- factor(pheno_tr$sexo)
cat("Carga OK — M_tr:", nrow(M_tr), "×", ncol(M_tr), "| Rubin:", USAR_RUBIN, "\n")

# Primary DMPs and high-confidence DMPs, used as reference
dmp_prim  <- as.data.frame(fread(file.path(DIR_PRIM, "dmp_resultados.csv")))
sig_prim  <- dmp_prim$cpg[dmp_prim$FDR < 0.05]           # primary DMPs
hic_prim  <- dmp_prim$cpg[dmp_prim$FDR < 0.05 & dmp_prim$hiconf]   # primary high-confidence DMPs
cat("Primario:", length(sig_prim), "DMP |", length(hic_prim), "hiconf\n")
# Reference adjusted effect (M scale) of each primary DMP; scenario columns are added in the loop
por_dmp <- data.frame(cpg = sig_prim, coef_ref = dmp_prim$logFC_M[match(sig_prim, dmp_prim$cpg)])

# ---- 2. Scenarios -----------------------------------------------------------------------------------
# Each scenario: subset, limma fit per imputation, Rubin pooling, bacon, BH FDR, delta-beta.
cohortes <- sort(unique(as.character(pheno_tr$gse)))
escenarios <- c(setNames(as.list(cohortes), paste0("sin_", cohortes)),   # one scenario per cohort
                list(sin_116378_y_152027 = c("GSE116378", "GSE152027")))   # additional scenario

extra_anc <- paste("+", paste(ANC_COVS, collapse = " + "))   # ancestry terms for the formula
resumen <- data.frame()                                   # one row per scenario

for (nm in names(escenarios)) {
  drop <- escenarios[[nm]]                                # cohort(s) left out
  keep <- !(as.character(pheno_tr$gse) %in% drop)
  p2 <- pheno_tr[keep, ]; p2$gse <- droplevels(factor(p2$gse))   # no empty cohort column in the design
  M2 <- M_tr[, keep, drop = FALSE]; b2 <- betas_tr[, keep, drop = FALSE]
  message(sprintf("  [%s] n=%d (quita %s: %d muestras)", nm, sum(keep),
                  paste(drop, collapse="+"), sum(!keep)))

  n_imp <- if (USAR_RUBIN) mids$m else 1L                 # 20 imputations or a single fit
  coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
  df_total <- NULL
  for (k in seq_len(n_imp)) {
    # Age of imputation k, or the point imputed age when Rubin is off
    p2$age_k <- if (USAR_RUBIN) complete(mids, k)[rownames(p2), "age_chronological"]
                else p2$age_epigenetic_imputed
    # Same design as script 06: cohort as fixed covariate, sex, age, cells, smoking score, ancestry
    di  <- model.matrix(as.formula(paste("~ group + gse + sexo + age_k +",
                        paste(CELL_COVS, collapse=" + "), "+ smoking_SSc", extra_anc)), data = p2)
    fit <- eBayes(lmFit(M2, di))                          # moderated t statistics
    coef_mat[, k] <- fit$coefficients[, "groupSCZ"]       # case-control effect
    se_mat[, k]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]   # moderated standard error
    df_total <- fit$df.total
  }
  # Rubin's rules: within- plus between-imputation variance
  coef <- rowMeans(coef_mat)
  if (n_imp > 1) {
    Wbar <- rowMeans(se_mat^2); B <- rowVars(coef_mat); Tvar <- Wbar + (1 + 1/n_imp) * B; SE <- sqrt(Tvar)
  } else { SE <- se_mat[, 1] }
  bcn <- bacon(effectsizes = coef, standarderrors = SE)  # bacon: corrects inflation and bias
  coef_bc <- as.numeric(es(bcn)); se_bc <- as.numeric(se(bcn)); p_bc <- as.numeric(pval(bcn))
  z_bc <- coef_bc / se_bc; lambda <- median(z_bc^2) / qchisq(0.5, 1)   # genomic lambda after bacon
  fdr <- p.adjust(p_bc, method = "BH")                    # Benjamini-Hochberg FDR
  scz <- p2$group == "SCZ"; cnt <- p2$group == "CNT"
  dbeta <- rowMeans(b2[, scz, drop=FALSE]) - rowMeans(b2[, cnt, drop=FALSE])   # unadjusted delta-beta
  hiconf <- abs(dbeta) >= DB_HICONF

  d <- data.frame(cpg = rownames(M2), coef_bc, se_bc, t = z_bc, p = p_bc, FDR = fdr,
                  delta_beta = dbeta, hiconf = hiconf, stringsAsFactors = FALSE)
  fwrite(d[d$FDR < 0.05, ], file.path(DIR_OUT, paste0("dmp_", nm, ".csv")))   # significant probes only

  # Comparison with the full model
  sig  <- d$cpg[d$FDR < 0.05]
  mm   <- merge(data.frame(cpg=d$cpg, db=d$delta_beta),   # delta-beta, scenario vs full model
                data.frame(cpg=dmp_prim$cpg, dbp=dmp_prim$delta_beta), by="cpg")
  hic_keep <- hic_prim[hic_prim %in% d$cpg]               # primary high-confidence DMPs still tested
  hic_surv <- sum(d$cpg %in% hic_keep & d$FDR < 0.05 & d$hiconf)   # and still significant and flagged
  # Adjusted effects of the primary DMPs: same sign, slope through the origin, correlation (as in 010b)
  me <- merge(por_dmp[, c("cpg", "coef_ref")], data.frame(cpg = d$cpg, coef = d$coef_bc), by = "cpg")
  mismo_signo <- sign(me$coef) == sign(me$coef_ref)
  pendiente   <- unname(coef(lm(me$coef ~ 0 + me$coef_ref))[1])   # 1 = same magnitude as full model
  por_dmp[[paste0("coef_", nm)]] <- d$coef_bc[match(por_dmp$cpg, d$cpg)]
  por_dmp[[paste0("FDR_", nm)]]  <- d$FDR[match(por_dmp$cpg, d$cpg)]
  # Summary row of this scenario
  resumen <- rbind(resumen, data.frame(
    escenario = nm, n = sum(keep), n_DMP = length(sig), n_hiconf = sum(d$FDR<0.05 & d$hiconf),
    lambda = round(lambda, 4), cor_dbeta = round(cor(mm$db, mm$dbp), 4),
    hiconf_prim_presentes = length(hic_keep), hiconf_sobreviven = hic_surv,
    DMP_prim_replican = length(intersect(sig, sig_prim)),   # primary DMPs significant here too
    signo_igual_pct = round(100 * mean(mismo_signo), 1), pendiente = round(pendiente, 3),
    cor_efectos = round(cor(me$coef, me$coef_ref), 3)))
  message(sprintf("    -> DMP=%d hiconf=%d λ=%.3f cor(Δβ)=%.4f hiconf_surv=%d/%d",
          length(sig), sum(d$FDR<0.05 & d$hiconf), lambda, cor(mm$db,mm$dbp), hic_surv, length(hic_keep)))
  rm(M2, b2, coef_mat, se_mat, d); gc()                   # free memory before the next scenario
}

# ---- 3. Summary and output --------------------------------------------------------------------------
cat("\n═══════════ LEAVE-ONE-COHORT-OUT — resumen (primario:", length(sig_prim), "DMP,",
    length(hic_prim), "hiconf) ═══════════\n")
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "loco_resumen.csv"))   # one row per scenario
fwrite(por_dmp, file.path(DIR_OUT, "loco_por_DMP_primario.csv"))   # per-DMP coefficient and FDR
cat("\nGuardado: loco_resumen.csv + dmp_<escenario>.csv (DMP significativos por corrida) en", DIR_OUT, "\n")
cat("\nLectura: si 'hiconf_sobreviven' ≈ 'hiconf_prim_presentes' y cor(Δβ)≈1 en las 7 corridas,",
    "\nningún hallazgo depende de una sola cohorte → robustez confirmada.\n")
