# 06c — Test statistics before bacon correction
#
# Script 06 corrects genomic inflation with bacon and saves only the corrected statistics. This script
# recomputes the Rubin-combined statistics from the same per-imputation checkpoint, with the same rules,
# and writes the p-values before bacon, so the QQ plot can show how much inflation was corrected and a
# reader can judge whether bacon was needed. It re-estimates nothing and does not re-run limma. Runs for
# the pre-registered primary model and for the main analysis (no smoking CpGs, chip row/column).
# Paper section: Differential Methylation (QQ plots, Figure 1)
# Input:  data_dmp_P2_anc/ and data_dmp_P2_anc_noSmk_pos/: dmp_perimput_checkpoint.rds and
#         dmp_resultados.csv (from 06)
# Output: dmp_prebacon.csv (cpg, t_prebacon, p_prebacon) in each of those two folders
# VARIANTE_MASKING: "P2" = pre-registered probe masking (397,194 probes); "P1" = aggressive masking
# USAR_ANCESTRIA:   TRUE = model with the four ancestry PCs (folder suffix _anc)
# EXCLUIR_COHORTE:  "" = all six cohorts; a GSE id = that cohort left out (folder suffix _5coh)

VARIANTE_MASKING <- "P2"   # probe-masking variant
USAR_ANCESTRIA   <- TRUE   # ancestry PCs in the model
EXCLUIR_COHORTE  <- ""   # empty = all six cohorts

SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)   # folder suffix: masking
SUF_ANC <- if (USAR_ANCESTRIA) "_anc" else ""   # folder suffix: ancestry
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""   # folder suffix: cohort left out
DIR_BASE <- paste0("../data_dmp", SUF_VAR, SUF_ANC, SUF_COH)   # pre-registered primary folder
# Arms: pre-registered primary, and main analysis (no smoking CpGs + chip row/column; QQ in Figure 1)
BRAZOS   <- c(primario = DIR_BASE, noSmk = paste0(DIR_BASE, "_noSmk_pos"))

for (brazo in names(BRAZOS)) {
  DIR_DMP <- BRAZOS[[brazo]]   # folder of this arm
  ck <- file.path(DIR_DMP, "dmp_perimput_checkpoint.rds")   # per-imputation fits saved by 06
  stopifnot(file.exists(ck))
  cp <- readRDS(ck)
  # Checkpoint elements: coef and se (probes x imputations), df (complete-data df per probe).
  stopifnot(all(c("coef", "se", "df") %in% names(cp)))
  coef_mat <- cp$coef; se_mat <- cp$se   # SCZ coefficient and its SE in each imputation
  m_imp    <- ncol(coef_mat)   # number of imputations (20)
  df_total <- cp$df   # residual df of the complete-data fit
  cat("\n[", brazo, "] checkpoint:", nrow(coef_mat), "sondas ×", m_imp, "imputaciones\n")

  # ---- Rubin's rules, identical to step 3b of 06 ---------------------------------------------------
  # Copied line by line: any difference would make the "before" QQ not comparable with the "after".
  coef  <- rowMeans(coef_mat)   # pooled coefficient
  Wbar  <- rowMeans(se_mat^2)   # within-imputation variance
  B     <- matrixStats::rowVars(coef_mat)   # between-imputation variance
  Tvar  <- Wbar + (1 + 1 / m_imp) * B   # total variance
  SE    <- sqrt(Tvar)   # pooled SE
  tstat <- coef / SE   # pooled t statistic

  B_safe <- pmax(B, .Machine$double.eps)   # avoids division by zero when imputations agree exactly
  r      <- (1 + 1 / m_imp) * B_safe / Wbar   # relative increase in variance due to missing age
  df_rub <- (m_imp - 1) * (1 + 1 / r)^2   # Rubin df
  df_use <- pmin(df_rub, df_total)   # capped at the complete-data df
  p_pre  <- 2 * pt(-abs(tstat), df = df_use)   # two-sided p-value before bacon

  lambda_pre <- median(tstat^2) / qchisq(0.5, 1)   # genomic inflation factor before bacon
  cat("λ ANTES de bacon:", round(lambda_pre, 4), "\n")

  # ---- Match to the 06 results by probe name -------------------------------------------------------
  # dmp_resultados.csv is sorted by FDR and the checkpoint is in genomic order; full coverage required.
  dmp <- read.csv(file.path(DIR_DMP, "dmp_resultados.csv"), stringsAsFactors = FALSE)
  idx <- match(dmp$cpg, rownames(coef_mat))   # checkpoint row of each result row
  stopifnot(!anyNA(idx), length(unique(idx)) == nrow(dmp), nrow(dmp) == nrow(coef_mat))
  tstat <- tstat[idx]; p_pre <- p_pre[idx]   # reorder to the results table

  lambda_post <- median(qnorm(dmp$p / 2)^2) / qchisq(0.5, 1)   # inflation after bacon, from the 06 p
  cat("λ DESPUÉS de bacon:", round(lambda_post, 4), "\n")
  cat("la corrección redujo la inflación en", round(lambda_pre - lambda_post, 4), "\n")

  out <- data.frame(cpg = dmp$cpg, t_prebacon = tstat, p_prebacon = p_pre)   # new file; 06 outputs untouched
  write.csv(out, file.path(DIR_DMP, "dmp_prebacon.csv"), row.names = FALSE)
  cat("guardado:", file.path(DIR_DMP, "dmp_prebacon.csv"), "\n")
  rm(cp, coef_mat, se_mat, dmp, out); gc()   # free memory before the next arm
}
