# 010g — Sensitivity to sample position on the chip
#
# In normalization (03), dasen did not apply its background refinement by the position of each sample on
# the chip (row R01-R06, column C01-C02), because the matrix columns were named by GSM accession. The
# position is available in the phenotype (column position, from the Sentrix code or the GEO metadata).
# This script tests whether that technical effect moves the result by adding chip row and column to the
# model of 06 as factors, as with any technical batch, for the pre-registered arm (full panel) and the arm
# without smoking CpGs. On the reference DMPs it reports sign concordance, slope through the origin,
# correlation of adjusted effects and how many remain at FDR < 0.05.
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds, data_combat_P2/ (phenotype_final.csv, age_mids.rds),
#         data_pheno/covariables_epigeneticas_P2.csv, data_pheno/phenodata_combinado_imputar_edad.csv,
#         data_dmp_P2_anc{,_noSmk}/dmp_resultados.csv (models without position),
#         refs/smoking_cpgs_exclusion_450k.txt
# Output: data_dmp_P2_anc_sens_posicion/sens_posicion_resumen.csv, por_DMP_<arm>_con_posicion.csv,
#         dmp_<arm>_con_posicion.csv (FDR < 0.05 only)

# ---- 1. Packages, paths and data (loaded once) -----------------------------------------------------
suppressMessages({ library(limma); library(bacon); library(data.table); library(readr)
                   library(matrixStats); library(mice) })
set.seed(7777777)                                            # fixed seed for reproducibility
# Run from fase_1/pipeline/ whether started there or at the repository root
if (!file.exists("../data_harmonized_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

VARIANTE_MASKING <- "P2"                                     # probe-masking variant of the primary analysis
DIR_HARM   <- paste0("../data_harmonized_", VARIANTE_MASKING)   # regRCPqn-harmonized betas (04)
DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING)       # phenotype and age imputations (05)
DIR_PHENO  <- "../data_pheno"                                   # covariates and chip position
DIR_OUT    <- "../data_dmp_P2_anc_sens_posicion"                # outputs of this script
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")     # EpiDISH leukocyte proportions
# Arms: reference table (fitted without chip position) and whether smoking CpGs leave the panel
ARMAS <- list(primario = list(ref = "../data_dmp_P2_anc/dmp_resultados.csv",       quitar_tabaco = FALSE),
              noSmk    = list(ref = "../data_dmp_P2_anc_noSmk/dmp_resultados.csv", quitar_tabaco = TRUE))
cpg_tabaco <- readLines("../refs/smoking_cpgs_exclusion_450k.txt")   # Joehanes 2016 + Bollepalli 2019 union

pheno <- as.data.frame(read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"), show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id                              # samples indexed by GSM
cov_ep <- as.data.frame(read_csv(file.path(DIR_PHENO,
           paste0("covariables_epigeneticas_", VARIANTE_MASKING, ".csv")), show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)   # EPISTRUCTURE ancestry PCs (04b)
i <- match(rownames(pheno), cov_ep$gsm_id)                   # align covariates to the phenotype
pheno[, ANC_COVS] <- cov_ep[i, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[i]                   # continuous EpiSmokEr score, as in 06

# Chip position (row and column) from the phenotype built in 01
pos <- as.data.frame(read_csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                              show_col_types = FALSE))
pheno$position <- pos$position[match(rownames(pheno), pos$gsm_id)]
stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))   # every sample has a valid position
pheno$fila    <- factor(substr(pheno$position, 1, 3))        # chip row R01-R06
pheno$columna <- factor(substr(pheno$position, 4, 6))        # chip column C01-C02
print(table(pheno$fila, pheno$group))                        # row balance between cases and controls
mids <- readRDS(file.path(DIR_COMBAT, "age_mids.rds"))       # 20 age imputations (mice, norm method)

# Harmonized betas, training set only, as M-values
b <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))$betas
ph <- pheno[colnames(b), ]                                   # phenotype in matrix column order
train <- ph$set == "train"                                   # 80% training partition (05)
b <- b[, train, drop = FALSE]; ph <- ph[train, ]; gc()
# M-values with betas bounded away from 0 and 1: limma assumes roughly normal residuals
eps <- 1e-6; b <- pmin(pmax(b, eps), 1 - eps); M_tr <- log2(b / (1 - b)); rm(b); gc()
ph$group <- factor(ph$group, levels = c("CNT", "SCZ"))       # controls as reference level
ph$sexo  <- factor(ph$sexo)
# Model of 06 plus chip row and column as technical factors
formula_txt <- paste("~ group + gse + sexo + age_k +", paste(CELL_COVS, collapse = " + "), "+ smoking_SSc +",
                     paste(ANC_COVS, collapse = " + "), "+ fila + columna")
cat("Carga OK — M_tr:", nrow(M_tr), "×", ncol(M_tr), "| modelo:", formula_txt, "\n")

# ---- 2. Refit with chip position, per arm -----------------------------------------------------------
resumen <- data.frame()                                      # one row per arm
for (arma in names(ARMAS)) {
  ref <- as.data.frame(fread(ARMAS[[arma]]$ref))             # reference results without position
  ref_dmp <- ref[ref$FDR < 0.05, c("cpg", "logFC_M")]        # reference DMPs and adjusted effects
  filas <- rownames(M_tr)
  if (ARMAS[[arma]]$quitar_tabaco) filas <- setdiff(filas, cpg_tabaco)   # arm without smoking CpGs
  stopifnot(setequal(filas, ref$cpg))                        # same panel as the reference arm
  M2 <- M_tr[filas, , drop = FALSE]
  n_imp <- mids$m                                            # number of age imputations
  coef_mat <- matrix(NA_real_, nrow(M2), n_imp); se_mat <- matrix(NA_real_, nrow(M2), n_imp)
  for (k in seq_len(n_imp)) {                                # one age imputation per iteration
    ph$age_k <- complete(mids, k)[rownames(ph), "age_chronological"]
    di  <- model.matrix(as.formula(formula_txt), data = ph)
    fit <- eBayes(lmFit(M2, di))                             # moderated linear model per CpG
    coef_mat[, k] <- fit$coefficients[, "groupSCZ"]          # diagnosis effect (M scale)
    se_mat[, k]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]   # moderated SE
  }
  coef <- rowMeans(coef_mat)                                 # pooled effect over imputations
  SE   <- sqrt(rowMeans(se_mat^2) + (1 + 1/n_imp) * rowVars(coef_mat))   # Rubin: within + between variance
  bcn  <- bacon(effectsizes = coef, standarderrors = SE)     # empirical-null correction of bias and inflation
  coef_bc <- as.numeric(es(bcn)); p_bc <- as.numeric(pval(bcn)); se_bc <- as.numeric(se(bcn))
  lambda  <- median((coef_bc / se_bc)^2) / qchisq(0.5, 1)    # genomic inflation after bacon
  fdr     <- p.adjust(p_bc, method = "BH")                   # Benjamini-Hochberg FDR
  d <- data.frame(cpg = rownames(M2), coef = coef_bc, p = p_bc, FDR = fdr)
  fwrite(d[d$FDR < 0.05, ], file.path(DIR_OUT, paste0("dmp_", arma, "_con_posicion.csv")))
  # Comparison on the reference DMPs
  m <- merge(data.frame(cpg = ref_dmp$cpg, coef_ref = ref_dmp$logFC_M), d, by = "cpg")
  mismo_signo <- sign(m$coef) == sign(m$coef_ref)            # same direction as the reference
  pendiente   <- unname(coef(lm(m$coef ~ 0 + m$coef_ref))[1])   # slope through the origin (1 = same size)
  fwrite(m, file.path(DIR_OUT, paste0("por_DMP_", arma, "_con_posicion.csv")))
  resumen <- rbind(resumen, data.frame(
    brazo = arma, escenario = "con_posicion", DMP_escenario = sum(fdr < 0.05), lambda = round(lambda, 3),
    DMP_ref = nrow(m), signo_igual_pct = round(100 * mean(mismo_signo), 1),
    pendiente = round(pendiente, 3), cor_efectos = round(cor(m$coef, m$coef_ref), 3),
    DMP_ref_siguen_FDR05 = sum(m$FDR < 0.05), pct_DMP_ref_siguen = round(100 * mean(m$FDR < 0.05), 1)))
  message(sprintf("  [%s] DMP=%d λ=%.3f signo=%.1f%% pendiente=%.3f cor=%.3f siguen=%d/%d", arma,
          sum(fdr < 0.05), lambda, 100 * mean(mismo_signo), pendiente, cor(m$coef, m$coef_ref),
          sum(m$FDR < 0.05), nrow(m)))
  rm(M2, coef_mat, se_mat, d, m); gc()                       # free memory before the next arm
}

# ---- 3. Summary and output --------------------------------------------------------------------------
# If the sign is kept and the slope is near 1, chip position does not explain the DMPs and the
# position-based refinement that dasen skipped would not have changed them.
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "sens_posicion_resumen.csv"))
cat("\nGuardado en", DIR_OUT, "\n")
