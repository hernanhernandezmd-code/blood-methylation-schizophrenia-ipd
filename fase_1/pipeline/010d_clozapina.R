# 010d — Sensitivity to clozapine
#
# Medication was not deposited for any cohort, so the model cannot adjust for clozapine. Compares DMPs
# with the CpGs that Hannon et al. (2021) associated with clozapine prescription (Supplementary Table 13,
# seven DMPs) or with treatment-resistant schizophrenia (Supplementary Table 14, schizophrenia DMPs with a
# significant treatment-resistance term), and recomputes FDR without them, for the pre-registered
# analysis and the main analysis (smoking-associated CpGs removed, chip row/column as covariates).
# Not part of the registered protocol.
# Paper section: Robustness and Sensitivity Analyses
# Input:  ../refs/hannon2021_elife_supplementary_file1.xlsx (eLife 10:e58430, Supplementary file 1);
#         ../data_dmp_P2_anc{,_noSmk_pos}/dmp_resultados.csv
# Output: ../data_dmp_P2_anc_clozapina/panel_clozapina_hannon2021.csv, clozapina_solape_<arm>.csv,
#         clozapina_resumen.csv
# FDR_MAX: FDR threshold used for the overlap and for the recomputed FDR (0.05).

suppressMessages({ library(data.table); library(readxl) })
# Run from the repository root or from fase_1/pipeline
if (!dir.exists("../refs") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

F_HANNON <- "../refs/hannon2021_elife_supplementary_file1.xlsx"   # cached supplementary file
DIR_OUT  <- "../data_dmp_P2_anc_clozapina"
dir.create(DIR_OUT, showWarnings = FALSE)
# Arms: pre-registered primary model, and main analysis (no smoking CpGs, chip position)
BRAZOS <- c(primario = "../data_dmp_P2_anc", principal = "../data_dmp_P2_anc_noSmk_pos")
FDR_MAX <- 0.05

# ---- 1. Panel from Hannon et al. (2021), Supplementary file 1 ---------------------------------------
# Supplementary Table 13: four header rows; column 1 = probe
s13 <- as.data.frame(read_excel(F_HANNON, sheet = "Supplementary Table 13", skip = 4, col_names = FALSE, .name_repair = "minimal"))
s13 <- s13[grepl("^cg", s13[[1]]), ]                      # keep probe rows, drop footnotes
cloz7 <- s13[[1]]                                         # clozapine-associated DMPs
cat("DMP de clozapina (S13):", length(cloz7), "\n")
stopifnot(length(cloz7) == 7)                             # guard against a changed sheet layout

# Supplementary Table 14: column 1 = probe, 6 = Significance, 9 = schizophrenia effect,
# 16/18 = mean difference and p-value of the treatment-resistance term (fixed effects)
s14 <- as.data.frame(read_excel(F_HANNON, sheet = "Supplementary Table 14", skip = 4, col_names = FALSE, .name_repair = "minimal"))
s14 <- s14[grepl("^cg", s14[[1]]), ]
stopifnot(nrow(s14) == 1048)                              # their 1,048 schizophrenia DMPs
print(table(s14[[6]]))                                    # sites per Significance category
s14 <- data.frame(cpg = s14[[1]], categoria = s14[[6]], efecto_scz_hannon = as.numeric(s14[[9]]),
                  efecto_trs_hannon = as.numeric(s14[[16]]), p_trs_hannon = as.numeric(s14[[18]]),
                  stringsAsFactors = FALSE)
trs <- s14$cpg[s14$categoria %in% c("TRS only", "SCZ & TRS")]   # sites with a resistance effect
cat("DMP de esquizofrenia con efecto TRS significativo (S14):", length(trs), "\n")

# Clozapine panel = Table 13 plus Table 14 sites with a significant treatment-resistance term
panel_cloz <- unique(c(cloz7, trs))
cat("Panel de clozapina (unión):", length(panel_cloz), "CpG\n")
fwrite(data.table(cpg = panel_cloz, origen = ifelse(panel_cloz %in% cloz7, "S13_clozapina", "S14_TRS")),
       file.path(DIR_OUT, "panel_clozapina_hannon2021.csv"))

# ---- 2. Overlap, enrichment and survival in each analysis arm ---------------------------------------
# Main-analysis DMP = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp written by script 06);
# a table without that column (pre-registered primary model) uses FDR < 0.05 alone
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
resumen <- list()                                         # one summary row per arm
for (brazo in names(BRAZOS)) {
  dmp <- fread(file.path(BRAZOS[[brazo]], "dmp_resultados.csv"))
  dmp$es_dmp_sel <- sel_dmp(dmp)                          # DMPs of the arm under its own definition
  sig <- dmp$es_dmp_sel
  en_panel <- dmp$cpg %in% panel_cloz                     # tested probes that belong to the panel

  # Overlap with the panel and enrichment (2 x 2 table: DMP yes/no by panel yes/no)
  n_panel_eval <- sum(en_panel)
  n_sol <- sum(sig & en_panel)                            # DMPs inside the panel
  esperado <- sum(sig) * n_panel_eval / nrow(dmp)         # overlap expected by chance
  fis <- fisher.test(matrix(c(n_sol, sum(sig & !en_panel), sum(!sig & en_panel), sum(!sig & !en_panel)), 2))

  # Direction at overlapping DMPs: crude difference against the treatment-resistance and schizophrenia effects
  sol <- merge(dmp[sig & en_panel, .(cpg, gene, delta_beta, FDR)], s14, by = "cpg", all.x = TRUE)
  sol$origen <- ifelse(sol$cpg %in% cloz7, "S13_clozapina", "S14_TRS")
  mismo_trs <- sum(sign(sol$delta_beta) == sign(sol$efecto_trs_hannon), na.rm = TRUE)
  mismo_scz <- sum(sign(sol$delta_beta) == sign(sol$efecto_scz_hannon), na.rm = TRUE)
  fwrite(sol[order(sol$FDR), ], file.path(DIR_OUT, paste0("clozapina_solape_", brazo, ".csv")))

  # Most of the panel consists of schizophrenia DMPs of that study, so overlap mostly reflects replication.
  # Within its 1,048 schizophrenia DMPs, compares replication of sites with a treatment-resistance effect
  # against sites with a schizophrenia effect only; higher replication of the former would point to clozapine.
  h <- merge(s14, dmp[, .(cpg, FDR, es_dmp_sel)], by = "cpg")   # their DMPs present in our panel
  h$nuestra_dmp <- h$es_dmp_sel                           # replicated here as a DMP
  h$grupo <- ifelse(h$categoria %in% c("TRS only", "SCZ & TRS"), "con_TRS", ifelse(h$categoria == "SCZ only", "SCZ_only", "otro"))
  tab_rep <- table(factor(h$grupo, c("con_TRS", "SCZ_only")), h$nuestra_dmp)
  rep_trs <- 100 * tab_rep["con_TRS", "TRUE"] / sum(tab_rep["con_TRS", ])     # % replicated, with TRS effect
  rep_scz <- 100 * tab_rep["SCZ_only", "TRUE"] / sum(tab_rep["SCZ_only", ])   # % replicated, SCZ only
  fis_rep <- fisher.test(tab_rep[, c("TRUE", "FALSE")])
  print(table(h$categoria, h$nuestra_dmp))

  # Clozapine DMPs (Supplementary Table 13) that are DMPs here
  s13_sol <- dmp[cpg %in% cloz7 & es_dmp_sel, .(cpg, gene, delta_beta, FDR)]

  # Removes the panel and recomputes Benjamini-Hochberg FDR over the remaining CpGs (bacon-corrected
  # p-values; the model is not refitted, since only a few CpGs are removed)
  resto <- dmp[!en_panel]
  resto$FDR_sin_cloz <- p.adjust(resto$p, method = "BH")
  n_dmp_sin <- sum(resto$FDR_sin_cloz < FDR_MAX)          # probes at FDR < 0.05 without the panel
  testables <- sum(sig & !en_panel)                       # DMPs outside the panel
  sobreviven <- sum(resto$es_dmp_sel & resto$FDR_sin_cloz < FDR_MAX)   # of those, still FDR < 0.05

  cat("\n====", brazo, "====\n")
  cat("Panel evaluado:", n_panel_eval, "de", length(panel_cloz), "| DMP:", sum(sig),
      "| solape:", n_sol, sprintf("(esperado %.1f; %.1f veces; Fisher p = %.2g)\n", esperado, n_sol / esperado, fis$p.value))
  cat("Solape con S13 (clozapina):", sum(sol$origen == "S13_clozapina"), "| con S14 (TRS):", sum(sol$origen == "S14_TRS"), "\n")
  cat("Mismo signo que el efecto TRS de Hannon:", mismo_trs, "de", sum(!is.na(sol$efecto_trs_hannon)),
      "| que su efecto SCZ:", mismo_scz, "de", sum(!is.na(sol$efecto_scz_hannon)), "\n")
  cat(sprintf("Replicación dentro de los DMP de Hannon: con efecto TRS %.1f %% vs SCZ only %.1f %% (OR = %.2f; Fisher p = %.2g)\n",
              rep_trs, rep_scz, fis_rep$estimate, fis_rep$p.value))
  cat("DMP de clozapina (S13) que son DMP nuestras:\n"); print(s13_sol)
  cat("Sin el panel:", n_dmp_sin, "DMP; de", testables, "DMP testables sobreviven", sobreviven,
      sprintf("(%.1f %%)\n", 100 * sobreviven / testables))
  # Summary row of this arm
  resumen[[brazo]] <- data.table(brazo = brazo, panel_cloz = length(panel_cloz), panel_evaluado = n_panel_eval,
                                 n_dmp = sum(sig), solape = n_sol, solape_S13 = sum(sol$origen == "S13_clozapina"),
                                 esperado = esperado, veces = n_sol / esperado, fisher_p = fis$p.value,
                                 fisher_or = unname(fis$estimate), mismo_signo_trs = mismo_trs,
                                 mismo_signo_scz = mismo_scz, n_hannon_trs_eval = sum(tab_rep["con_TRS", ]),
                                 n_hannon_scz_eval = sum(tab_rep["SCZ_only", ]), rep_trs_pct = rep_trs,
                                 rep_scz_pct = rep_scz, rep_or = unname(fis_rep$estimate), rep_p = fis_rep$p.value,
                                 s13_cpg = paste(s13_sol$cpg, collapse = ";"), n_dmp_sin_cloz = n_dmp_sin,
                                 testables = testables, sobreviven = sobreviven,
                                 pct_sobreviven = 100 * sobreviven / testables)
}
resumen <- rbindlist(resumen)
fwrite(resumen, file.path(DIR_OUT, "clozapina_resumen.csv"))
print(resumen)

# ---- 3. Positions of the interleukin-10 term in the panel -------------------------------------------
# CpGs behind the interleukin-10 enrichment term; does the clozapine panel touch them or JAK1?
il10 <- c("cg02532928", "cg00697440", "cg07316313", "cg23635663", "cg20301308", "cg24718015", "cg24452282",
          "cg21370522", "cg11484872", "cg10315334")
cat("CpG del término IL-10 en el panel de clozapina:", sum(il10 %in% panel_cloz), "de", length(il10), "\n")
cat("JAK1 (cg20301308) en el panel:", "cg20301308" %in% panel_cloz, "\n")
