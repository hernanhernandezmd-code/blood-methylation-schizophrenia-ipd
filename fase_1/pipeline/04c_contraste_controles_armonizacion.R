# 04c — Between-cohort contrast among controls after harmonization
#
# Tests, CpG by CpG, whether control samples still differ between cohorts after regRCPqn, as the
# registered protocol specifies, and checks whether the DMPs of the pre-registered primary model
# fall on flagged positions. The DMP model already carries cohort as a covariate, so this is a
# quality check, not an extra correction.
# Paper section: Normalization and Harmonization
# Input:  data_harmonized_P2/betas_harmonized.rds, data_combat_P2/phenotype_final.csv,
#         data_dmp_P2_anc/dmp_resultados.csv
# Output: contraste_controles_post_regRCPqn.csv (per-CpG F tests and flag),
#         contraste_controles_resumen.csv/.rds (summary), in the harmonized folder
# Model A: one-way ANOVA F on cohort (unadjusted). Model B: F for the cohort block adjusted
# for sex, age, cell proportions, smoking score and two ancestry PCs.

# ---- 1. Paths and parameters ------------------------------------------------
DIR_HARM  <- "../data_harmonized_P2"                # regRCPqn betas from 04
DIR_PHENO <- "../data_combat_P2"                    # phenotype_final.csv with the 80/20 split (05)
DIR_DMP   <- "../data_dmp_P2_anc"                   # pre-registered primary DMP table (06)
DIR_MAN   <- "../data_ml_P2_anc"              # manifest with n_concord (directional concordance)
N_COH     <- 6                                      # cohorts; "full concordance" = sign agrees in all
DIR_OUT   <- DIR_HARM                               # results sit next to the betas they describe

FDR_FLAG  <- 0.05                              # protocol "significant" = BH FDR < 0.05
CHUNK     <- 40000L                            # rows per block in model B (memory)
set.seed(20260604)                             # project seed

# ---- 2. Phenotype: training-set controls ------------------------------------
# All inference in the study is restricted to the 80% training partition.
ph <- read.csv(file.path(DIR_PHENO, "phenotype_final.csv"), stringsAsFactors = FALSE)
ph <- ph[ph$set == "train" & ph$group == "CNT", ]   # controls only: no disease effect to confound
cat(sprintf("controles train: %d en %d cohortes\n", nrow(ph), length(unique(ph$gse))))
print(table(ph$gse))                                # controls per cohort

# ---- 3. Harmonized betas to M-values ----------------------------------------
obj <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))
if (is.list(obj) && !is.matrix(obj)) betas <- obj$betas else betas <- obj  # accept list or bare matrix
rm(obj); invisible(gc())
stopifnot(all(ph$gsm_id %in% colnames(betas)))      # every control must be in the matrix
betas <- betas[, ph$gsm_id]                         # columns in phenotype order
invisible(gc())
cat(sprintf("matriz de controles: %d CpG x %d muestras\n", nrow(betas), ncol(betas)))

# M-values for the tests: closer to normal residuals than betas.
M <- log2(pmin(pmax(betas, 1e-6), 1 - 1e-6) /  # betas clipped away from 0 and 1
          (1 - pmin(pmax(betas, 1e-6), 1 - 1e-6)))

# Cohort means on the beta scale, to report the size of the differences.
cohortes <- sort(unique(ph$gse))
med_beta <- matrix(NA_real_, nrow(betas), length(cohortes),
                   dimnames = list(rownames(betas), cohortes))  # CpGs x cohorts
for (g in cohortes) {
  med_beta[, g] <- rowMeans(betas[, ph$gse == g, drop = FALSE])  # mean beta of that cohort's controls
}
rango_beta <- apply(med_beta, 1, max) - apply(med_beta, 1, min)  # max |delta beta| between cohorts
rm(betas); invisible(gc())

# ---- 4. Model A: one-way ANOVA F on cohort, vectorized ----------------------
# Between- and within-cohort sums of squares accumulated cohort by cohort for all CpGs at once.
n_tot <- ncol(M); k <- length(cohortes)             # samples and groups
gran_media <- rowMeans(M)                           # grand mean per CpG
ssb <- numeric(nrow(M)); ssw <- numeric(nrow(M))    # between and within sums of squares
for (g in cohortes) {
  Xg  <- M[, ph$gse == g, drop = FALSE]             # controls of cohort g
  n_g <- ncol(Xg)
  m_g <- rowMeans(Xg)                               # cohort mean per CpG
  ssb <- ssb + n_g * (m_g - gran_media)^2           # weighted squared deviation of the cohort mean
  ssw <- ssw + rowSums((Xg - m_g)^2)                # spread around the cohort mean
}
names(ssb) <- names(ssw) <- rownames(M)
F_A   <- (ssb / (k - 1)) / (ssw / (n_tot - k))      # F with k - 1 and n - k df
p_A   <- pf(F_A, k - 1, n_tot - k, lower.tail = FALSE)
fdr_A <- p.adjust(p_A, method = "BH")               # Benjamini-Hochberg over the panel
flag  <- fdr_A < FDR_FLAG                           # protocol flag: residual between-cohort difference
cat(sprintf("\nMODELO A (sin ajustar): %s de %s CpG marcados (FDR<%.2f) = %.1f%%\n",
            format(sum(flag), big.mark = ","), format(length(flag), big.mark = ","),
            FDR_FLAG, 100 * mean(flag)))
cat(sprintf("  |Δβ| máx entre cohortes en los marcados: mediana %.4f · p90 %.4f · máx %.4f\n",
            median(rango_beta[flag]), quantile(rango_beta[flag], .9), max(rango_beta[flag])))

# ---- 5. Model B: adjusted F for the cohort block ----------------------------
# F from the residual sums of squares of a reduced (covariates) and a full (covariates + cohort)
# design, computed with residual-maker matrices in blocks of rows. Separates technical residue
# from real demographic differences between cohorts (age, sex, smoking, ancestry, cell mix).
covs <- c("sexo", "age_epigenetic_imputed", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro",
          "smoking_SSc", "ancestry_PC1", "ancestry_PC2")   # sex, age, six cell types, smoking, ancestry
X_red  <- model.matrix(~ ., data = ph[, covs])                  # covariates only
X_full <- model.matrix(~ . + gse, data = ph[, c(covs, "gse")])  # covariates + cohort
P_red  <- diag(n_tot) - X_red  %*% solve(crossprod(X_red))  %*% t(X_red)   # I - hat matrix, reduced
P_full <- diag(n_tot) - X_full %*% solve(crossprod(X_full)) %*% t(X_full)  # I - hat matrix, full
df1 <- ncol(X_full) - ncol(X_red)              # df of the cohort block (k - 1)
df2 <- n_tot - ncol(X_full)                    # residual df of the full model

rss_red <- numeric(nrow(M)); rss_full <- numeric(nrow(M))
for (ini in seq(1L, nrow(M), by = CHUNK)) {
  fin <- min(ini + CHUNK - 1L, nrow(M))        # last row of this block
  Y   <- M[ini:fin, , drop = FALSE]
  rss_red [ini:fin] <- rowSums((Y %*% P_red)  * Y)   # y' P y = residual sum of squares per CpG
  rss_full[ini:fin] <- rowSums((Y %*% P_full) * Y)
}
F_B   <- ((rss_red - rss_full) / df1) / (rss_full / df2)   # extra-sum-of-squares F for cohort
p_B   <- pf(F_B, df1, df2, lower.tail = FALSE)
fdr_B <- p.adjust(p_B, method = "BH")
cat(sprintf("MODELO B (ajustado): %s CpG con FDR<%.2f = %.1f%% (el residuo NO explicado por demografía)\n",
            format(sum(fdr_B < FDR_FLAG), big.mark = ","), FDR_FLAG, 100 * mean(fdr_B < FDR_FLAG)))

# ---- 6. Do the primary DMPs fall on flagged CpGs? ---------------------------
# Counts flagged DMPs (all, and those concordant in every cohort) and tests enrichment of DMPs
# among flagged CpGs relative to the panel (Fisher's exact test).
dmp <- read.csv(file.path(DIR_DMP, "dmp_resultados.csv"), stringsAsFactors = FALSE)
sig <- dmp$cpg[dmp$FDR < 0.05]                      # DMPs of the pre-registered primary model
man <- read.csv(file.path(DIR_MAN, "loci_fase2_manifest.csv"), stringsAsFactors = FALSE)
core <- intersect(sig, man$locus_id[man$n_concord == N_COH])   # same sign in all six cohorts

en_panel  <- sig[sig %in% names(F_A)]               # DMPs present in this contrast
marcadas  <- en_panel[flag[en_panel]]               # DMPs on flagged CpGs
cat(sprintf("\nSENSIBILIDAD DE LA FIRMA\n  %d DMP: %d marcadas por el contraste (%.0f%%) vs %.1f%% del panel\n",
            length(sig), length(marcadas), 100 * length(marcadas) / length(en_panel), 100 * mean(flag)))
cat(sprintf("  núcleo %d/%d (n=%d): %d marcadas\n", N_COH, N_COH, length(core), sum(flag[core])))
# 2 x 2 table: flagged / not flagged, DMP / rest of the panel.
ft <- fisher.test(matrix(c(sum(flag[en_panel]), sum(!flag[en_panel]),
                           sum(flag) - sum(flag[en_panel]),
                           sum(!flag) - sum(!flag[en_panel])), 2))
cat(sprintf("  enriquecimiento de las DMP en CpGs marcados: OR = %.2f, p = %.3g\n",
            ft$estimate, ft$p.value))
if (length(marcadas)) {                             # list flagged DMPs with gene and size
  cat("  DMP marcadas (gen | FDR contraste | rango beta entre cohortes):\n")
  for (c1 in marcadas) {
    cat(sprintf("    %s  %-12s FDR=%.3g  rango=%.4f\n", c1,
                dmp$gene[dmp$cpg == c1][1], fdr_A[c1], rango_beta[c1]))
  }
}

# ---- 7. Save ----------------------------------------------------------------
# One row per CpG: unadjusted and adjusted tests, between-cohort range and protocol flag.
out <- data.frame(cpg = names(F_A), F_crudo = F_A, p_crudo = p_A, FDR_crudo = fdr_A,
                  F_ajust = F_B, p_ajust = p_B, FDR_ajust = fdr_B,
                  rango_beta_cohortes = rango_beta, flag_protocolo = flag,
                  row.names = NULL)
write.csv(out, file.path(DIR_OUT, "contraste_controles_post_regRCPqn.csv"), row.names = FALSE)
# With several hundred controls nearly every CpG reaches significance, so the summary also reports
# the size of the between-cohort range and where the DMPs sit in its distribution.
r_dmp    <- rango_beta[en_panel]                                # between-cohort range at each DMP
pctl_dmp <- sapply(r_dmp, function(v) mean(rango_beta <= v))     # percentile of each DMP in the panel
w        <- wilcox.test(r_dmp, rango_beta)                      # DMP ranges vs panel ranges
resumen <- c(n_controles = n_tot, n_cpg = nrow(out), n_flag = sum(flag),
             pct_flag = 100 * mean(flag), pct_flag_ajust = 100 * mean(fdr_B < FDR_FLAG),
             dmp_marcadas = length(marcadas), core_marcadas = sum(flag[core]),
             p_fisher = ft$p.value,
             rango_mediano_panel = median(rango_beta), rango_p90_panel = quantile(rango_beta, .9),
             rango_mediano_dmp = median(r_dmp), pctl_mediano_dmp = 100 * median(pctl_dmp),
             p_wilcoxon_rango = w$p.value)
write.csv(data.frame(metrica = names(resumen), valor = unname(resumen)),
          file.path(DIR_OUT, "contraste_controles_resumen.csv"), row.names = FALSE)  # figures cited in text
saveRDS(resumen, file.path(DIR_OUT, "contraste_controles_resumen.rds"))
cat("\nGuardado:", file.path(DIR_OUT, "contraste_controles_post_regRCPqn.csv"), "\n")
print(round(resumen, 4))
