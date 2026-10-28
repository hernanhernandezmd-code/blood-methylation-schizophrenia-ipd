# 010e — Overlap with the 1,048 schizophrenia DMPs of Hannon et al. (2021)
#
# Compares DMPs with the schizophrenia DMPs of the largest previous meta-analysis (eLife 10:e58430,
# Supplementary Table 3, two-stage random-effects model): number tested here, overlap, agreement in
# direction and correlation of covariate-adjusted effects, for the pre-registered analysis and the main
# analysis (smoking-associated CpGs removed, chip row/column as covariates). Four of the six cohorts
# analyzed here (GSE80417, GSE84727, GSE147221, GSE152027) also contributed to that meta-analysis, so the
# comparison is not an independent replication. Their Twins and Sweden cohorts are not in this study; the
# direction check against them is a diagnostic only, because their 1,048 sites were selected with those
# cohorts included.
# Paper section: Robustness and Sensitivity Analyses (Results); Discussion
# Input:  ../refs/hannon2021_elife_supplementary_file1.xlsx;
#         ../data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv;
#         ../data_dmp_P2_anc_concordancia/concordancia_por_DMP_{primario,noSmk}.csv (script 010c)
# Output: ../data_dmp_P2_anc_hannon2021/hannon2021_resumen.csv, hannon2021_solape_<arm>.csv,
#         hannon2021_posiciones_nombradas.csv
# NOMBRADAS: positions named in the text (AHRR, JAK1), looked up in their list with the top-FDR genes.

suppressMessages({ library(data.table); library(readxl) })
# Run from the repository root or from fase_1/pipeline
if (!dir.exists("../refs") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

F_HANNON <- "../refs/hannon2021_elife_supplementary_file1.xlsx"   # cached supplementary file
DIR_OUT  <- "../data_dmp_P2_anc_hannon2021"
dir.create(DIR_OUT, showWarnings = FALSE)
# Arms: pre-registered primary model, and main analysis (no smoking CpGs, chip position)
BRAZOS <- c(primario = "../data_dmp_P2_anc", principal = "../data_dmp_P2_anc_noSmk_pos")
# Per-DMP directional concordance across the six cohorts, written by script 010c
CONC   <- c(primario = "../data_dmp_P2_anc_concordancia/concordancia_por_DMP_primario.csv",
            principal = "../data_dmp_P2_anc_concordancia/concordancia_por_DMP_noSmk.csv")
FDR_MAX <- 0.05
NOMBRADAS <- c(AHRR = "cg05575921", JAK1 = "cg20301308")

# ---- 1. The 1,048 DMPs of Hannon et al. (Supplementary Table 3) -------------------------------------
# Four header rows. Columns: 1 probe, 4 gene, 6 number of cohorts, 7/9 mean difference and p (random
# effects), 14-25 mean difference and SE per cohort (UCL, Aberdeen, Twins, IoPPN, Dublin, Sweden)
s3 <- suppressWarnings(as.data.frame(read_excel(F_HANNON, sheet = "Supplementary Table 3", skip = 4, col_names = FALSE,
                               .name_repair = "minimal")))   # warnings: X-chromosome row with "NA"
s3 <- s3[grepl("^cg", s3[[1]]), ]                         # keep probe rows, drop footnotes
# Differences are percentage points of methylation, covariate-adjusted; text "NA" for cohorts without data
hannon <- suppressWarnings(data.frame(cpg = s3[[1]], gen_hannon = s3[[4]], n_coh_hannon = as.numeric(s3[[6]]),
                     dif_hannon = as.numeric(s3[[7]]), p_hannon = as.numeric(s3[[9]]),
                     dif_twins = as.numeric(s3[[18]]), ee_twins = as.numeric(s3[[19]]),
                     dif_sweden = as.numeric(s3[[24]]), ee_sweden = as.numeric(s3[[25]]), stringsAsFactors = FALSE))
# Twins and Sweden combined by inverse variance (the two cohorts of Hannon et al. not included here)
w_t <- 1 / hannon$ee_twins^2; w_s <- 1 / hannon$ee_sweden^2
# One cohort alone when the other has no estimate, otherwise the weighted mean
hannon$dif_indep <- ifelse(is.na(w_t), hannon$dif_sweden, ifelse(is.na(w_s), hannon$dif_twins,
                           (w_t * hannon$dif_twins + w_s * hannon$dif_sweden) / (w_t + w_s)))
cat("Con efecto en Twins o Sweden:", sum(!is.na(hannon$dif_indep)), "\n")
cat("DMP de esquizofrenia de Hannon 2021:", nrow(hannon), "\n")
stopifnot(nrow(hannon) == 1048)                           # guard against a changed sheet layout
cat("Hipermetiladas en casos:", sum(hannon$dif_hannon > 0), sprintf("(%.1f %%)\n", 100 * mean(hannon$dif_hannon > 0)))

# ---- 2. Overlap, direction and effect size in each analysis arm -------------------------------------
# Main-analysis DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp written by script 06);
# a table without that column (pre-registered primary model) uses FDR < 0.05 alone
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
resumen <- list()                                         # one summary row per arm
for (brazo in names(BRAZOS)) {
  dmp <- fread(file.path(BRAZOS[[brazo]], "dmp_resultados.csv"))
  conc <- fread(CONC[[brazo]])
  dmp$es_dmp_sel <- sel_dmp(dmp)                          # DMPs of the arm under its own definition
  sig <- dmp$es_dmp_sel
  en_h <- dmp$cpg %in% hannon$cpg                         # tested probes on their list

  # Tested positions and overlap (2 x 2 table: DMP here yes/no by Hannon DMP yes/no)
  n_eval <- sum(en_h)
  n_sol <- sum(sig & en_h)                                # their DMPs that are DMPs here
  fis <- fisher.test(matrix(c(n_sol, sum(sig & !en_h), sum(!sig & en_h), sum(!sig & !en_h)), 2))
  pct_base <- 100 * mean(sig)                               # percentage of DMPs in the whole panel

  # Direction and size: adjusted difference here (%) against their adjusted mean difference (%)
  m <- merge(dmp[en_h, .(cpg, gene, delta_beta, delta_beta_adj, FDR, es_dmp_sel)], hannon, by = "cpg")
  m$db_pct <- 100 * m$delta_beta_adj                      # adjusted delta-beta, same units as theirs
  m$mismo_signo <- sign(m$db_pct) == sign(m$dif_hannon)
  m$dmp_nuestra <- m$es_dmp_sel
  bin_todas <- binom.test(sum(m$mismo_signo), nrow(m))    # sign agreement against 50%, all tested
  bin_dmp   <- binom.test(sum(m$mismo_signo[m$dmp_nuestra]), sum(m$dmp_nuestra))   # our DMPs only
  r_todas <- cor(m$db_pct, m$dif_hannon)                  # effect-size correlation, all tested
  r_dmp   <- cor(m$db_pct[m$dmp_nuestra], m$dif_hannon[m$dmp_nuestra])
  pend    <- unname(coef(lm(db_pct ~ 0 + dif_hannon, data = m))[1])   # slope through the origin
  # Twins + Sweden comparison. Not independent of selection: their DMPs were selected with these cohorts
  # included, by a random-effects model that penalizes heterogeneity
  mi <- m[!is.na(m$dif_indep), ]
  mi$mismo_signo_indep <- sign(mi$db_pct) == sign(mi$dif_indep)
  bin_indep <- binom.test(sum(mi$mismo_signo_indep), nrow(mi))
  ind_dmp <- mi$mismo_signo_indep[mi$dmp_nuestra]         # same check restricted to our DMPs
  fwrite(m[order(m$FDR), ], file.path(DIR_OUT, paste0("hannon2021_solape_", brazo, ".csv")))

  # DMPs concordant in all six cohorts that are among their DMPs
  c6 <- conc$cpg[conc$concordancia_ajustada == 6]
  n_c6 <- length(c6); n_c6_h <- sum(c6 %in% hannon$cpg)

  cat("\n====", brazo, "====\n")
  cat("Evaluadas:", n_eval, "de 1048 | DMP nuestras entre ellas:", n_sol, sprintf("(%.1f %%; base del panel %.2f %%; OR = %.1f; Fisher p = %.2g)\n",
      100 * n_sol / n_eval, pct_base, fis$estimate, fis$p.value))
  cat(sprintf("Mismo signo: %d de %d evaluadas (%.1f %%, p = %.2g) | %d de %d DMP nuestras (%.1f %%)\n",
      sum(m$mismo_signo), nrow(m), 100 * mean(m$mismo_signo), bin_todas$p.value,
      sum(m$mismo_signo[m$dmp_nuestra]), sum(m$dmp_nuestra), 100 * mean(m$mismo_signo[m$dmp_nuestra])))
  cat(sprintf("Correlación de efectos: r = %.2f (evaluadas), %.2f (DMP nuestras); pendiente %.2f\n", r_todas, r_dmp, pend))
  cat(sprintf("Réplica independiente (Twins + Sweden): mismo signo %d de %d evaluadas (%.1f %%, p = %.2g) | %d de %d DMP nuestras (%.1f %%)\n",
      sum(mi$mismo_signo_indep), nrow(mi), 100 * mean(mi$mismo_signo_indep), bin_indep$p.value,
      sum(ind_dmp), length(ind_dmp), 100 * mean(ind_dmp)))
  cat("Concordantes 6/6:", n_c6, "| en la lista de Hannon:", n_c6_h, "\n")
  # Summary row of this arm
  resumen[[brazo]] <- data.table(brazo = brazo, hannon_n = nrow(hannon), evaluadas = n_eval, solape = n_sol,
                                 pct_solape = 100 * n_sol / n_eval, pct_base = pct_base,
                                 fisher_or = unname(fis$estimate), fisher_p = fis$p.value,
                                 signo_igual = sum(m$mismo_signo), pct_signo_igual = 100 * mean(m$mismo_signo),
                                 binom_p = bin_todas$p.value, signo_igual_dmp = sum(m$mismo_signo[m$dmp_nuestra]),
                                 pct_signo_igual_dmp = 100 * mean(m$mismo_signo[m$dmp_nuestra]),
                                 r_efectos = r_todas, r_efectos_dmp = r_dmp, pendiente = pend,
                                 n_conc6 = n_c6, n_conc6_en_hannon = n_c6_h, indep_n = nrow(mi),
                                 indep_signo_igual = sum(mi$mismo_signo_indep),
                                 pct_indep_signo_igual = 100 * mean(mi$mismo_signo_indep), indep_binom_p = bin_indep$p.value,
                                 indep_n_dmp = length(ind_dmp), pct_indep_signo_igual_dmp = 100 * mean(ind_dmp))
}
resumen <- rbindlist(resumen)
fwrite(resumen, file.path(DIR_OUT, "hannon2021_resumen.csv"))
print(resumen)

# ---- 3. Positions named in the text -----------------------------------------------------------------
main <- fread(file.path(BRAZOS[["principal"]], "dmp_resultados.csv"))
sig_m <- main[sel_dmp(main)][order(FDR)]                  # main-analysis DMPs by increasing FDR
# Named genes only: skip clone-based identifiers, LINC and antisense transcripts
nombre_gen <- function(g) !grepl("^(AC|AL|AP|BX|CT|FP|Z)[0-9]{5,}|^LINC|-AS[0-9]", g) & g != "" & !is.na(g)
top4 <- head(sig_m[nombre_gen(sapply(strsplit(gene, ";"), `[`, 1))], 4)   # four lowest-FDR named genes
prim <- fread(file.path(BRAZOS[["primario"]], "dmp_resultados.csv"))
nombradas <- rbind(data.table(etiqueta = paste0("top_FDR_", sapply(strsplit(top4$gene, ";"), `[`, 1)), cpg = top4$cpg),
                   data.table(etiqueta = names(NOMBRADAS), cpg = unname(NOMBRADAS)))
nombradas$en_hannon <- nombradas$cpg %in% hannon$cpg      # on their list of 1,048
# Their difference and p-value, and the unadjusted delta-beta (%) of the pre-registered model
nombradas <- merge(nombradas, hannon[, c("cpg", "dif_hannon", "p_hannon")], by = "cpg", all.x = TRUE, sort = FALSE)
nombradas <- merge(nombradas, prim[, .(cpg, db_primario = 100 * delta_beta)], by = "cpg", all.x = TRUE, sort = FALSE)
print(nombradas)
fwrite(nombradas, file.path(DIR_OUT, "hannon2021_posiciones_nombradas.csv"))
