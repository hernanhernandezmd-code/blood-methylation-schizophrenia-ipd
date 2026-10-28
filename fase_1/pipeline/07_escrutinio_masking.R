# 07 — Technical flags for significant DMPs ("flag and consider")
#
# Classifies each significant DMP as core (probe also in the stricter P1 panel, which added Naeem 2014
# and any-SNP masking) or recovered (probe removed by P1 and kept by P2), and adds technical flags for
# each probe instead of deleting them (Byrska-Bishop et al., 2020).
# Paper section: Robustness and Sensitivity Analyses
# Input:  data_dmp_P2_anc[_noSmk|_noAnc|_noAge][_pos]/dmp_resultados.csv (from 06);
#         data_dmp/dmp_resultados.csv (P1 masking variant); probe lists in refs/
# Output: same folder as input, escrutinio_masking_P2.csv
# EXCLUIR_CPG_TABACO: TRUE (default) reads the analysis without smoking CpGs (_noSmk, main analysis of
#   the paper); FALSE reads the pre-registered primary model on the full panel.
# EXCLUIR_CPG_ANCESTRIA / EXCLUIR_CPG_EDAD: FALSE (default); TRUE reads the analysis without ancestry /
#   clock CpGs. Only one exclusion flag may be TRUE, and all flags must match the run of 06.
# USAR_POSICION: TRUE (default) reads the model with chip row/column covariates (_pos, main analysis).

# ---- 1. Packages, probe lists and annotation ---------------------------------------------------

# Installs missing packages (CRAN or Bioconductor) and loads them.
InsLoad.pks <- function(pqks, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  new.pkg <- pqks[!(pqks %in% rownames(installed.packages()))]                # packages not yet installed
  if (length(new.pkg) > 0) BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)
  invisible(sapply(pqks, require, character.only = TRUE, quietly = TRUE))   # load all, quietly
}
InsLoad.pks(c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "maxprobes", "data.table", "readr"))

DIR_DMP_P1 <- "../data_dmp"                              # DMPs under the P1 masking variant (fixed reference)

# The folder suffix is rebuilt from the same flags used in 06. The probe panel does not depend on
# these flags, so the core/recovered split is the same in every analysis.
VARIANTE_MASKING   <- "P2"                               # pre-registered masking variant
USAR_ANCESTRIA     <- TRUE                               # ancestry PCs in the model (_anc)
EXCLUIR_CPG_TABACO <- TRUE                               # TRUE = main analysis without smoking CpGs (_noSmk)
EXCLUIR_CPG_ANCESTRIA <- FALSE                           # sensitivity without ancestry CpGs (_noAnc)
EXCLUIR_CPG_EDAD      <- FALSE                           # sensitivity without clock CpGs (_noAge)
EXCLUIR_COHORTE    <- ""                                 # empty = all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""   # suffix only when a cohort is left out
SUF_DMP <- paste0("_", VARIANTE_MASKING)                 # suffixes accumulate in the order used by 06
if (USAR_ANCESTRIA)      SUF_DMP <- paste0(SUF_DMP, "_anc")
SUF_DMP <- paste0(SUF_DMP, SUF_COH)
if (EXCLUIR_CPG_TABACO)  SUF_DMP <- paste0(SUF_DMP, "_noSmk")
if (EXCLUIR_CPG_ANCESTRIA) SUF_DMP <- paste0(SUF_DMP, "_noAnc")
if (EXCLUIR_CPG_EDAD)      SUF_DMP <- paste0(SUF_DMP, "_noAge")
USAR_POSICION <- TRUE                                    # chip row/column in the model (_pos, main analysis)
if (USAR_POSICION)        SUF_DMP <- paste0(SUF_DMP, "_pos")

DIR_DMP_P2 <- paste0("../data_dmp", SUF_DMP)             # DMP table of the active analysis
DIR_REFS   <- "../refs"                                  # cached probe lists
FDR_CUT    <- 0.05                                       # significance threshold
cat("Variante activa:", SUF_DMP, "| escrutinio sobre:", DIR_DMP_P2, "\n")

# 450K annotation and the lists that define which probes P1 removed and P2 keeps.
anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
lista_naeem   <- tryCatch(readLines(file.path(DIR_REFS, "naeem2014_problematic_450k.txt")),
                          warning = function(w) character(0))         # Naeem 2014 problematic probes
lista_zhou    <- tryCatch(readLines(file.path(DIR_REFS, "zhou2017_mask_general_450k.txt")),
                          warning = function(w) character(0))         # Zhou 2017 general mask
lista_chen    <- maxprobes::xreactive_probes(array_type = "450K")     # cross-reactive, Chen + Benton
snp_probe     <- rownames(anno)[!is.na(anno$Probe_rs)]                 # SNP anywhere in the probe
snp_cpgsbe    <- rownames(anno)[!is.na(anno$CpG_rs) | !is.na(anno$SBE_rs)]   # SNP at the CpG or SBE
# Probes removed only by P1 (Naeem or any-SNP) and kept in P2.
p1_extra_excl <- unique(c(lista_naeem, snp_probe))
cat("Listas: Naeem", length(lista_naeem), "| Zhou", length(lista_zhou), "| Chen", length(lista_chen),
    "| SNP_probe", length(snp_probe), "| SNP_CpG/SBE", length(snp_cpgsbe), "\n")

# ---- 2. Load DMPs and define recovered probes --------------------------------------------------

dmp1 <- as.data.frame(data.table::fread(file.path(DIR_DMP_P1, "dmp_resultados.csv")))   # P1 results
dmp2 <- as.data.frame(data.table::fread(file.path(DIR_DMP_P2, "dmp_resultados.csv")))   # P2 results
cat("Panel evaluado — P1:", nrow(dmp1), "sondas | P2:", nrow(dmp2), "sondas\n")

# DMP of the main analysis = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp from 06);
# tables without that column (pre-registered primary model) fall back to FDR < 0.05.
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
sig1 <- dmp1$cpg[dmp1$FDR < FDR_CUT]                    # significant probes in P1
sig2 <- dmp2$cpg[sel_dmp(dmp2)]                          # DMPs in P2 (FDR and effect threshold if present)
cat("DMP FDR<0.05 — P1:", length(sig1), "| P2:", length(sig2), "\n")

# Recovered probes: in the P2 panel but not evaluated in P1, which removed them.
recuperadas_panel <- intersect(dmp2$cpg, setdiff(p1_extra_excl, dmp1$cpg))
cat("Sondas recuperadas en el panel P2 (no evaluadas en P1):", length(recuperadas_panel), "\n")

# ---- 3. Classify DMPs and add technical flags --------------------------------------------------
# P2 already masks Chen, Zhou and CpG/SBE SNPs, so recovered DMPs usually carry only the Naeem and
# Probe_rs flags (SNP in the probe body, a lesser concern than a SNP at the measured CpG).
# en_P1_sig (significant in the P1 run) is informative only.

d <- dmp2[dmp2$FDR < FDR_CUT, ]                          # significant probes of P2 (FDR only)
d$nivel     <- ifelse(d$cpg %in% recuperadas_panel, "recuperada", "core")   # core vs recovered
d$en_P1_sig <- d$cpg %in% sig1                           # also significant in P1
d$en_panelP1 <- d$cpg %in% dmp1$cpg                      # probe present in the P1 panel

d$flag_naeem    <- d$cpg %in% lista_naeem                # resequencing: mapping or repeat issue
d$flag_probe_rs <- d$cpg %in% snp_probe                  # SNP in the probe body
d$flag_snp_cpgsbe <- d$cpg %in% snp_cpgsbe               # SNP at the CpG/SBE (high concern; rare in P2)
d$flag_chen     <- d$cpg %in% lista_chen                 # cross-reactive, Chen + Benton
d$flag_zhou     <- d$cpg %in% lista_zhou                 # mapping/design issue, Zhou
d$n_flags <- with(d, flag_naeem + flag_probe_rs + flag_snp_cpgsbe + flag_chen + flag_zhou)   # 0 = no flag
d <- d[order(-d$hiconf, d$FDR), ]                        # high-confidence first, then by FDR

# escrutinio_masking_P2.csv: significant probes with level (core/recovered) and flag columns.
readr::write_csv(d, file.path(DIR_DMP_P2, "escrutinio_masking_P2.csv"))
cat("Guardado: escrutinio_masking_P2.csv —", nrow(d), "DMPs de P2 clasificados\n")

# ---- 4. Summaries ------------------------------------------------------------------------------

cat("===== DMPs P2 (", nrow(d), ") por nivel =====\n")
print(table(nivel = d$nivel))                            # core vs recovered counts
cat("\n===== hiconf (|Δβ|≥0.02) por nivel =====\n")
print(table(nivel = d$nivel, hiconf = ifelse(d$hiconf, "sí", "no")))   # level x high-confidence

# Overlap between P1 and P2 significant probes.
cat("\n===== Consistencia P1 ↔ P2 =====\n")
persisten <- intersect(sig1, sig2)                       # significant in both variants
cat("DMPs de P1 (", length(sig1), ") que SIGUEN significativos en P2:", length(persisten),
    sprintf("(%.1f%%)\n", 100 * length(persisten) / max(length(sig1), 1)))
cat("DMPs NUEVOS en P2 (no significativos en P1):", length(setdiff(sig2, sig1)), "\n")

# Flag counts among recovered probes.
cat("\n===== 'Flag-and-consider' de las RECUPERADAS =====\n")
rec <- d[d$nivel == "recuperada", ]                      # recovered significant probes
cat("Recuperadas significativas:", nrow(rec),
    "| de ellas hiconf:", sum(rec$hiconf),
    "| con ≥1 bandera técnica:", sum(rec$n_flags > 0), "\n")
cat("Banderas entre recuperadas — Naeem:", sum(rec$flag_naeem),
    "| Probe_rs:", sum(rec$flag_probe_rs),
    "| SNP CpG/SBE (ALTA sospecha):", sum(rec$flag_snp_cpgsbe),
    "| Chen:", sum(rec$flag_chen), "| Zhou:", sum(rec$flag_zhou), "\n")

# Largest-effect recovered high-confidence probes, for one-by-one review.
cat("\n===== TOP DMPs hiconf RECUPERADOS (a escrutar uno a uno) =====\n")
top <- rec[rec$hiconf, c("cpg","gene","delta_beta","FDR","en_P1_sig",
                         "flag_naeem","flag_probe_rs","flag_snp_cpgsbe","n_flags")]
print(utils::head(top[order(-abs(top$delta_beta)), ], 20), row.names = FALSE)   # top 20 by |delta-beta|

# Core, high-confidence probes with no technical flag.
cat("\n===== Núcleo ROBUSTO: DMPs core hiconf y limpios (0 banderas) =====\n")
robusto <- d[d$nivel == "core" & d$hiconf & d$n_flags == 0, ]
cat("DMPs core + hiconf + sin banderas técnicas:", nrow(robusto), "\n")
print(utils::head(robusto[order(robusto$FDR), c("cpg","gene","delta_beta","FDR","en_P1_sig")], 15),
      row.names = FALSE)                                 # top 15 by FDR
