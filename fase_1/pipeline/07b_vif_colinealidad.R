# 07b — Collinearity of model covariates (VIF)
#
# Computes the variance inflation factor of every term in the differential methylation design of
# script 06, the correlation among cell-type proportions and the correlation of each cell type with
# diagnosis, as the protocol asks to watch for collinearity of cell proportions with diagnosis. Uses
# phenotype data only; methylation values are not read.
# Paper section: Covariates
# Input:  data_combat_P2/phenotype_final.csv (from 05); data_pheno/covariables_epigeneticas_P2.csv
#         (04/04b); data_pheno/phenodata_combinado_imputar_edad.csv (chip position)
# Output: data_dmp_P2_anc_noSmk_pos/vif_colinealidad_P2.csv and vif_celula_grupo_P2.csv
# Flags must match script 06. VARIANTE_MASKING "P2" (default); USAR_ANCESTRIA TRUE (default) adds the
#   ancestry PCs; EXCLUIR_CPG_TABACO TRUE (default) / EXCLUIR_CPG_ANCESTRIA / EXCLUIR_CPG_EDAD change
#   the probe set, not the covariates, so they only set the output folder; EXCLUIR_COHORTE "" (default)
#   keeps all six cohorts; USAR_POSICION TRUE (default) adds chip row and column (_pos, main analysis).

# ---- 1. Packages and phenotype -----------------------------------------------------------------

if (!requireNamespace("car", quietly = TRUE))
  install.packages("car", repos = "https://cloud.r-project.org")
suppressMessages(library(car)); suppressMessages(library(readr))

VARIANTE_MASKING   <- "P2"                               # pre-registered masking variant
USAR_ANCESTRIA     <- TRUE                               # EPISTRUCTURE ancestry PCs in the design
EXCLUIR_CPG_TABACO <- TRUE                               # probe-set flag: folder only (_noSmk, main analysis)
EXCLUIR_CPG_ANCESTRIA <- FALSE                           # probe-set flag: folder only (_noAnc)
EXCLUIR_CPG_EDAD      <- FALSE                           # probe-set flag: folder only (_noAge)
EXCLUIR_COHORTE    <- ""                                 # empty = all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""   # suffix only when a cohort is left out
SUF_DMP <- paste0("_", VARIANTE_MASKING)                 # suffixes accumulate in the order used by 06
if (USAR_ANCESTRIA)      SUF_DMP <- paste0(SUF_DMP, "_anc")
SUF_DMP <- paste0(SUF_DMP, SUF_COH)
if (EXCLUIR_CPG_TABACO)  SUF_DMP <- paste0(SUF_DMP, "_noSmk")
if (EXCLUIR_CPG_ANCESTRIA) SUF_DMP <- paste0(SUF_DMP, "_noAnc")
if (EXCLUIR_CPG_EDAD)      SUF_DMP <- paste0(SUF_DMP, "_noAge")
USAR_POSICION <- TRUE                                    # chip row/column in the design (_pos)
if (USAR_POSICION)        SUF_DMP <- paste0(SUF_DMP, "_pos")

DIR_COMBAT <- paste0("../data_combat_", VARIANTE_MASKING, SUF_COH)   # phenotype from script 05
DIR_PHENO  <- "../data_pheno"                                # epigenetic covariates (04/04b)
DIR_DMP    <- paste0("../data_dmp", SUF_DMP)                 # output folder of the active analysis
dir.create(DIR_DMP, showWarnings = FALSE, recursive = TRUE)
CELL_COVS  <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")   # eosinophils omitted: proportions sum to ~1
COV_SMK    <- "smoking_SSc"                                  # continuous EpiSmokEr score (pre-registered)

pheno <- as.data.frame(readr::read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"),
                                       show_col_types = FALSE))
# Smoking score and ancestry PCs come from the covariate table, joined by GSM ID as in 06.
cov_ep <- as.data.frame(readr::read_csv(
  file.path(DIR_PHENO, paste0("covariables_epigeneticas_", VARIANTE_MASKING, SUF_COH, ".csv")),
  show_col_types = FALSE))
ANC_COVS <- if (USAR_ANCESTRIA) grep("^ancestry_PC", colnames(cov_ep), value = TRUE) else character(0)
idx <- match(pheno$gsm_id, cov_ep$gsm_id)                    # row of each sample in the covariate table
if (length(ANC_COVS)) pheno[, ANC_COVS] <- cov_ep[idx, ANC_COVS, drop = FALSE]
pheno$smoking_SSc <- cov_ep$smoking_SSc[idx]                 # numeric score, not a factor

# Training set only (the samples analysed in 06).
pheno <- pheno[pheno$set == "train", ]
pheno$group <- factor(pheno$group); pheno$gse <- factor(pheno$gse); pheno$sexo <- factor(pheno$sexo)
cat("Variante activa:", SUF_DMP, "| tabaco:", COV_SMK,
    "| PCs ancestría:", length(ANC_COVS), "\n")
cat("Muestras (train) para VIF:", nrow(pheno), "\n")

# ---- 2. VIF of the script 06 design ------------------------------------------------------------
# Same formula as 06: group + cohort + sex + imputed age + 6 cell types + smoking score + ancestry PCs
# (+ chip row and column). VIF depends only on the design matrix, so a random response is used to fit
# the model. Multilevel factors use GVIF^(1/(2*df)) (Fox and Monette); its square is compared with the
# usual VIF thresholds (5 watch, 10 problematic).
extra_anc <- if (length(ANC_COVS)) paste("+", paste(ANC_COVS, collapse = " + ")) else ""
if (USAR_POSICION) {                                    # chip row and column, as in 06
  pos <- as.data.frame(readr::read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv", show_col_types = FALSE))
  pheno$position <- pos$position[match(pheno$gsm_id, pos$gsm_id)]   # position on the slide, e.g. R03C01
  stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))        # 450K slide: 6 rows x 2 columns
  pheno$fila    <- factor(substr(pheno$position, 1, 3))            # row R01-R06
  pheno$columna <- factor(substr(pheno$position, 4, 6))            # column C01-C02
  extra_anc <- paste(extra_anc, "+ fila + columna")
}
form <- as.formula(paste("~ group + gse + sexo + age_epigenetic_imputed +",
                         paste(CELL_COVS, collapse = " + "), "+", COV_SMK, extra_anc))
cat("Fórmula del VIF:", deparse1(form), "\n\n")
set.seed(7777777)                                       # fixed seed for the random response
pheno$y_dummy <- rnorm(nrow(pheno))                     # arbitrary response, only to call vif()
mod  <- lm(update(form, y_dummy ~ .), data = pheno)

v <- car::vif(mod)                                      # GVIF, Df, GVIF^(1/(2*Df))
v <- as.data.frame(v)
v$termino <- rownames(v)                                # model term
col_adj <- grep("GVIF.*1", colnames(v), value = TRUE)[1]     # column name differs across car versions
v$VIF_equivalente <- v[[col_adj]]^2                          # comparable with classical VIF thresholds
v <- v[order(-v$VIF_equivalente), c("termino", "GVIF", "Df", col_adj, "VIF_equivalente")]

cat("===== VIF de la matriz de diseño del 06 (P2) =====\n")
cat("Regla: VIF_equivalente <5 ok | 5-10 vigilar | >10 problemático\n\n")
print(v, row.names = FALSE, digits = 3)

alerta <- v$termino[v$VIF_equivalente >= 5]              # terms at or above the watch threshold
cat("\nTérminos con VIF_equivalente >= 5:",
    if (length(alerta)) paste(alerta, collapse = ", ") else "NINGUNO — sin colinealidad problemática", "\n")

readr::write_csv(v, file.path(DIR_DMP, "vif_colinealidad_P2.csv"))
cat("Guardado: vif_colinealidad_P2.csv\n")

# ---- 3. Cell-cell and cell-diagnosis correlation -----------------------------------------------
# (a) correlation among the 6 cell types in the model; (b) correlation of each cell type, eosinophils
# included, with diagnosis (SCZ = 1, CNT = 0) and its p-value.
cells <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
cat("===== Correlación entre células (6 del modelo) =====\n")
print(round(cor(pheno[, cells]), 2))                     # cell-cell Pearson correlation

cat("\n===== Correlación célula <-> diagnóstico (SCZ vs CNT) =====\n")
scz01 <- as.integer(pheno$group == "SCZ")                # diagnosis as 0/1
tab <- data.frame(celula = c(cells, "Eosino"))
tab$r_con_grupo <- sapply(tab$celula, function(cc) cor(pheno[[cc]], scz01))   # point-biserial r
tab$p <- sapply(tab$celula, function(cc) cor.test(pheno[[cc]], scz01)$p.value)   # test of r = 0
tab <- tab[order(-abs(tab$r_con_grupo)), ]               # strongest association first
print(tab, row.names = FALSE, digits = 3)
cat("\n(r pequeño = célula NO confundida con diagnóstico → el ajuste por célula es estable y no roba señal)\n")

readr::write_csv(tab, file.path(DIR_DMP, "vif_celula_grupo_P2.csv"))
cat("Guardado: vif_celula_grupo_P2.csv\n")
