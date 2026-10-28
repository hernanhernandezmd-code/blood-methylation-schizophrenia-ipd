# 06 — Differential methylation: DMPs (limma + bacon) and DMRs (dmrff)
#
# One-stage IPD model fitted with limma on M-values of the training set (80%, 2,343 samples), with cohort
# as a fixed covariate (not ComBat: two-step correction inflates false positives with unbalanced designs),
# Rubin's rules over the 20 age imputations and bacon for genomic inflation. Effects are also reported
# in beta scale, unadjusted and adjusted for the model covariates. DMRs with dmrff.
# Paper section: Differential Methylation; Robustness and Sensitivity Analyses
# Input:  data_harmonized_P2/betas_harmonized.rds; data_combat_P2/phenotype_final.csv, age_mids.rds;
#         data_pheno/covariables_epigeneticas_P2.csv (ancestry PCs, smoking score);
#         data_pheno/phenodata_combinado_imputar_edad.csv (chip position); refs/ (CpG panels, Zhou 2017)
# Output: dmp_resultados.csv, dmr_resultados.csv, dmp_perimput_checkpoint.rds in the folder of the arm:
#   data_dmp_P2_anc/            pre-registered primary: EXCLUIR_CPG_TABACO = FALSE, USAR_POSICION = FALSE
#   data_dmp_P2_anc_noSmk_pos/  main analysis (default): EXCLUIR_CPG_TABACO = TRUE, USAR_POSICION = TRUE
#   data_dmp_P2_anc_noAnc/, data_dmp_P2_anc_noAge/  sensitivities: EXCLUIR_CPG_ANCESTRIA or
#                               EXCLUIR_CPG_EDAD = TRUE, with EXCLUIR_CPG_TABACO and USAR_POSICION = FALSE
# VARIANTE_MASKING: "P2" = pre-registered probe masking (397,194 probes); must match 02-05.
# EXCLUIR_COHORTE: "" = all six cohorts; a GSE id = that cohort left out (suffix _5coh).
# USAR_ANCESTRIA: TRUE adds the EPISTRUCTURE ancestry PCs to the model (suffix _anc).
# EXCLUIR_CPG_TABACO / _ANCESTRIA / _EDAD: TRUE removes published smoking / ancestry / clock CpGs
#   before limma, keeping the matching covariate (only one may be TRUE).
# USAR_POSICION: TRUE adds chip row and column as technical covariates (suffix _pos).
# UMBRAL_DB_ADJ: |adjusted Δβ| required for es_dmp; 0.005 in the main analysis, 0 elsewhere.

# ---- 1. Helper functions and packages ------------------------------------------------------------

# Installs missing packages (CRAN/Bioconductor) and loads them.
InsLoad.pks <- function(pqks, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))                 # BiocManager is a prerequisite
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  ip <- rownames(installed.packages())                                  # packages already installed
  new.pkg <- pqks[!(pqks %in% ip)]                                      # missing ones
  if (length(new.pkg) > 0)                                              # install only what is missing
    BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)
  out <- sapply(pqks, require, character.only = TRUE, quietly = TRUE)   # load and report TRUE/FALSE
  print(out)
  if (length(new.pkg) > 0)
    message("Instalados ", length(new.pkg), " paquetes nuevos.")
  else
    message("Todos los paquetes requeridos ya estaban instalados.")
  invisible(out)
}


# Wraps functions so that set.seed() runs on every call (not called here: limma and bacon are
# deterministic, and the age imputations come fixed from 05).
sealseed <- function(func_names, seed = 7777777) {
  for (func_name in func_names) {                                       # each function to wrap
    original_func <- get(func_name, envir = .GlobalEnv)                 # keep the original
    new_func <- function(...) {
      set.seed(as.integer(get("seed", envir = environment())))         # fix the seed on each call
      get("original_func", envir = environment())(...)                  # then call the original
    }
    environment(new_func) <- new.env(parent = environment(original_func))
    assign("original_func", original_func, envir = environment(new_func))
    assign("seed", seed, envir = environment(new_func))
    assign(func_name, new_func, envir = .GlobalEnv)                     # replace it in the global env
  }
}


# dmrff region calling plus a region effect in beta scale. dmrff's estimate is standardized, so the
# region effect (meanActualDbeta) is the mean per-CpG Δβ (SCZ − CNT) of the CpGs inside each region
# (Chacón & Hernández 2025). Also reports directional concordance and a |Δβ| flag (DB_pass).
# All per-CpG arguments must share order and length; CpGs are mapped to regions by position.
dmrffDB <- function(estimate, se, p.value, methylation, chr, pos, dbeta,
                    maxgap = 500, p.cutoff = 0.05, dbeta.cutoff = 0.05) {
  reg <- dmrff::dmrff(estimate = estimate, se = se, p.value = p.value,     # region calling + p-values
                      methylation = methylation, chr = chr, pos = pos,
                      maxgap = maxgap, p.cutoff = p.cutoff)
  meanDB <- nUp <- nDown <- numeric(nrow(reg))
  for (i in seq_len(nrow(reg))) {                                          # region Δβ by overlap
    en <- chr == reg$chr[i] & pos >= reg$start[i] & pos <= reg$end[i]      # CpGs inside the region
    d  <- dbeta[en]
    meanDB[i] <- mean(d, na.rm = TRUE)                                     # mean Δβ of the region
    nUp[i]    <- sum(d > 0, na.rm = TRUE)                                  # CpGs higher in SCZ
    nDown[i]  <- sum(d < 0, na.rm = TRUE)                                  # CpGs lower in SCZ
  }
  reg$meanActualDbeta  <- meanDB
  reg$n_up <- nUp; reg$n_down <- nDown
  reg$prop_concordante <- pmax(nUp, nDown) / pmax(nUp + nDown, 1)          # share in the dominant direction
  reg$DB_pass <- abs(reg$meanActualDbeta) >= dbeta.cutoff                  # effect flag, does not filter
  reg[order(reg$p.adjust, -reg$prop_concordante), ]                        # sorted by significance
}


# Gene annotation from Zhou et al. 2017 (GENCODE v36), joined by probe ID, instead of the Illumina
# manifest. is_promoter = any transcript TSS within [prom_lo, prom_hi] bp (negative distToTSS = upstream).
anotarGenesZhou <- function(cpgs,
                            ruta    = file.path(if (exists("DIR_REFS")) DIR_REFS else "../refs",
                                                "zhou2017_gene_annotation_450k.tsv.gz"),
                            prom_lo = -1500, prom_hi = 500) {   # asymmetric promoter window
  z <- data.table::fread(ruta, sep = "\t", header = TRUE, na.strings = "NA",
                         colClasses = "character", data.table = FALSE)
  z <- z[match(cpgs, z$probeID), ]                            # align to cpgs (NA if absent)
  d_list  <- strsplit(z$distToTSS, ";", fixed = TRUE)         # TSS distances, one per transcript
  is_prom <- vapply(d_list, function(s) {                     # any TSS inside the window?
    d <- suppressWarnings(as.numeric(s)); any(d >= prom_lo & d <= prom_hi, na.rm = TRUE)
  }, logical(1))
  d_near  <- vapply(d_list, function(s) {                     # signed distance to the nearest TSS
    d <- suppressWarnings(as.numeric(s)); if (all(is.na(d))) NA_real_ else d[which.min(abs(d))]
  }, numeric(1))
  data.frame(cpg = cpgs, gene = z$genesUniq, transcriptTypes = z$transcriptTypes,
             dist_nearest_TSS = d_near, is_promoter = is_prom, island = z$CGIposition,
             stringsAsFactors = FALSE)
}


pqks_BC <- c(
  "limma", "bacon",                                     # DMP: moderated linear model + inflation correction
  "DMRcate",                                            # used by 06b
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"        # 450K manifest (hg19 coordinates)
)
pqks_CRAN <- c(
  "data.table", "matrixStats",                          # fast I/O and matrix statistics
  "readr", "stringr", "ggplot2",
  "mice", "remotes"                                     # mice for the imputations; remotes for dmrff
)

InsLoad.pks(c(pqks_BC, pqks_CRAN))

# dmrff is installed from GitHub (perishky/dmrff): inverse-variance meta-analysis of CpGs in the same
# direction, adjusted for correlation between neighbouring CpGs.
if (!requireNamespace("dmrff", quietly = TRUE))
  remotes::install_github("perishky/dmrff", upgrade = "never")
library(dmrff)

# Stop with a clear message if a critical package did not load.
criticos <- c("limma", "bacon", "dmrff", "DMRcate", "mice",
              "IlluminaHumanMethylation450kanno.ilmn12.hg19")
faltan   <- criticos[!criticos %in% loadedNamespaces()]
if (length(faltan) > 0)
  stop("Paquetes críticos no cargados: ", paste(faltan, collapse = ", "),
       ". Reinicia el kernel y re-ejecuta esta celda.")
cat("Paquetes críticos cargados OK:", paste(criticos, collapse = ", "), "\n")


# ---- 2. Paths, flags, data, training subset and design -------------------------------------------

set.seed(7777777)   # global seed

# Anchor the working directory to fase_1/pipeline/ wherever the script is launched from.
if (!file.exists("../data_harmonized") && dir.exists("fase_1/pipeline")) {
  setwd("fase_1/pipeline")
}

VARIANTE_MASKING <- "P2"                                   # probe-masking variant; must match 02-05
SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)
EXCLUIR_COHORTE <- ""                                      # empty = all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""

DIR_HARMONIZED <- paste0("../data_harmonized", SUF_VAR, SUF_COH)  # regRCPqn betas from 04
DIR_COMBAT     <- paste0("../data_combat", SUF_VAR, SUF_COH)      # final phenotype + age imputations (05)
DIR_PHENO      <- "../data_pheno"
DIR_DMP        <- paste0("../data_dmp", SUF_VAR)         # output folder; arm suffixes added below
DIR_REFS       <- "../refs"                              # cached annotations and CpG panels
dir.create(DIR_DMP, showWarnings = FALSE, recursive = TRUE)

# Cell fractions in the model; eosinophils are left out because the seven fractions sum to about 1.
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

# Ancestry PCs (EPISTRUCTURE, from 04b): they lower residual inflation that cohort alone does not absorb.
USAR_ANCESTRIA <- TRUE
ANC_COVS <- character(0)                                 # filled in 2a
if (USAR_ANCESTRIA) {
  DIR_DMP <- paste0(DIR_DMP, "_anc")
}
DIR_DMP <- paste0(DIR_DMP, SUF_COH); dir.create(DIR_DMP, showWarnings = FALSE, recursive = TRUE)

# Smoking covariate: continuous EpiSmokEr smoking score, as pre-registered.
COV_SMK <- "smoking_SSc"

# CpG-exclusion arms (protocol 13.4.10): remove a published CpG panel before limma while keeping the
# matching covariate, to test that the SCZ effect survives without the probes that carry that factor.
#   _noSmk: Joehanes 2016 + Bollepalli 2019 smoking CpGs, smoking score kept (main analysis of the paper)
#   _noAnc: EPISTRUCTURE ancestry CpGs (Rahmani), ancestry PCs kept
#   _noAge: Horvath 2013 + Hannum 2013 clock CpGs, imputed age kept
# Leave-one-cohort-out is run with EXCLUIR_COHORTE.
EXCLUIR_CPG_TABACO    <- TRUE   # TRUE = without smoking CpGs (_noSmk, main analysis)
EXCLUIR_CPG_ANCESTRIA <- FALSE  # TRUE = without ancestry CpGs (_noAnc)
EXCLUIR_CPG_EDAD      <- FALSE  # TRUE = without clock CpGs (_noAge)
stopifnot("Solo UNA sensibilidad de exclusión de CpGs (a/b/c) puede estar activa a la vez" =
          sum(EXCLUIR_CPG_TABACO, EXCLUIR_CPG_ANCESTRIA, EXCLUIR_CPG_EDAD) <= 1)
# Read the active exclusion panel into cpg_excl and add its suffix to the output folder.
cpg_excl <- character(0); sens_lbl <- "ninguna (PRIMARIO)"
if (EXCLUIR_CPG_TABACO) {
  DIR_DMP <- paste0(DIR_DMP, "_noSmk"); cpg_excl <- readLines(file.path(DIR_REFS, "smoking_cpgs_exclusion_450k.txt"))
  sens_lbl <- paste0("(a) tabaco — ", length(cpg_excl), " CpGs (Joehanes∪Bollepalli)")
} else if (EXCLUIR_CPG_ANCESTRIA) {
  DIR_DMP <- paste0(DIR_DMP, "_noAnc"); cpg_excl <- readLines(file.path(DIR_REFS, "epistructure_ancestry_cpgs_450k.txt"))
  sens_lbl <- paste0("(b) ancestría — ", length(cpg_excl), " CpGs (EPISTRUCTURE/Rahmani)")
} else if (EXCLUIR_CPG_EDAD) {
  DIR_DMP <- paste0(DIR_DMP, "_noAge"); cpg_excl <- readLines(file.path(DIR_REFS, "clock_cpgs_exclusion_450k.txt"))
  sens_lbl <- paste0("(c) edad — ", length(cpg_excl), " CpGs (Horvath∪Hannum)")
}
# Chip position: row (R01-R06) and column (C01-C02) as technical factors. Not pre-registered (declared
# deviation); they remove technical noise without changing effect direction and add power.
USAR_POSICION <- TRUE
if (USAR_POSICION) DIR_DMP <- paste0(DIR_DMP, "_pos")
# Effect-size threshold of the main analysis (declared deviation): a DMP needs FDR < 0.05 and
# |adjusted Δβ| >= 0.5%, half the mean effect of the Hannon et al. 2021 DMPs. Other arms use FDR only.
UMBRAL_DB_ADJ <- if (EXCLUIR_CPG_TABACO && USAR_POSICION) 0.005 else 0
cat("Umbral de |Δβ ajustado| para ser DMP:", UMBRAL_DB_ADJ, "\n")
dir.create(DIR_DMP, showWarnings = FALSE, recursive = TRUE)
COV_SMK <- "smoking_SSc"                                 # smoking score stays in every arm
cat("Tabaquismo — covariable:", COV_SMK, "| SENSIBILIDAD CpG activa:", sens_lbl,
    "| salida:", basename(DIR_DMP), "\n")


# ---- 2a. Load harmonized betas, phenotype and age imputations ------------------------------------

harm   <- readRDS(file.path(DIR_HARMONIZED, "betas_harmonized.rds"))   # list(betas, cohorte)
betas  <- harm$betas                                    # probes x samples (colnames = GSM IDs)
rm(harm); gc()

pheno <- as.data.frame(readr::read_csv(file.path(DIR_COMBAT, "phenotype_final.csv"),
                                       show_col_types = FALSE))
rownames(pheno) <- pheno$gsm_id

# Ancestry PCs (from 04b) are appended from the covariate file.
if (USAR_ANCESTRIA) {
  cov_anc <- as.data.frame(readr::read_csv(
    file.path(DIR_PHENO, paste0("covariables_epigeneticas", "_", VARIANTE_MASKING, SUF_COH, ".csv")),
    show_col_types = FALSE))
  ANC_COVS <- grep("^ancestry_PC", colnames(cov_anc), value = TRUE)   # ancestry_PC1..K
  if (length(ANC_COVS) == 0) stop("USAR_ANCESTRIA=TRUE pero no hay columnas ancestry_PC — corre 04b primero.")
  pcs_m <- cov_anc[match(rownames(pheno), cov_anc$gsm_id), ANC_COVS, drop = FALSE]   # matched by GSM
  pheno[, ANC_COVS] <- pcs_m
  cat("Ancestría: añadidos", length(ANC_COVS), "PCs al fenotipo (", paste(ANC_COVS, collapse = ", "), ")\n")
}

# The smoking score is not in phenotype_final.csv; it is taken from the covariate file (from 04).
if (!exists("cov_anc")) cov_anc <- as.data.frame(readr::read_csv(
  file.path(DIR_PHENO, paste0("covariables_epigeneticas", "_", VARIANTE_MASKING, SUF_COH, ".csv")),
  show_col_types = FALSE))
if (!"smoking_SSc" %in% colnames(cov_anc))
  stop("No hay columna smoking_SSc en covariables_epigeneticas_", VARIANTE_MASKING,
       ".csv — corre el 04 con EpiSmokEr method='all'.")
pheno$smoking_SSc <- cov_anc$smoking_SSc[match(rownames(pheno), cov_anc$gsm_id)]   # matched by GSM
cat("Tabaquismo: añadido puntaje continuo smoking_SSc al fenotipo | NA:",
    sum(is.na(pheno$smoking_SSc)), "\n")

# Chip position (row and column) from the 01 phenotype (Sentrix code or GEO metadata).
if (USAR_POSICION) {
  pos <- as.data.frame(readr::read_csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                                       show_col_types = FALSE))
  pheno$position <- pos$position[match(rownames(pheno), pos$gsm_id)]   # e.g. R03C01
  stopifnot(all(grepl("^R0[1-6]C0[1-2]$", pheno$position)))   # every sample has a valid position
  pheno$fila    <- factor(substr(pheno$position, 1, 3))   # row R01-R06
  pheno$columna <- factor(substr(pheno$position, 4, 6))   # column C01-C02
  cat("Posición en el chip: añadidas fila y columna al fenotipo\n"); print(table(pheno$fila, pheno$columna))
}

pheno <- pheno[colnames(betas), ]                       # align to the column order of betas

mids  <- readRDS(file.path(DIR_COMBAT, "age_mids.rds")) # mice object, 20 age imputations

cat("betas:", nrow(betas), "sondas ×", ncol(betas), "muestras\n")
cat("¿fenotipo alineado con betas?:", identical(rownames(pheno), colnames(betas)), "\n")
cat("Reparto set:\n"); print(table(pheno$set, useNA = "ifany"))
cat("Reparto set × grupo:\n"); print(table(pheno$set, pheno$group))


# ---- 2b. Training subset (80%) and M-values ------------------------------------------------------
# The 20% test set does not enter inference. limma models M-values (closer to normal and
# homoscedastic; Du et al. 2010); betas are kept to express effects as Δβ.

train <- pheno$set == "train"
cat("Muestras de entrenamiento (80%):", sum(train), "| reservadas (20%):", sum(!train), "\n")

betas_tr <- betas[, train, drop = FALSE]                # training betas (for Δβ)
pheno_tr <- pheno[train, ]
rm(betas); gc()

# CpG-exclusion arms: drop the active panel before M-values and Δβ (empty in the primary model).
if (length(cpg_excl)) {
  n0 <- nrow(betas_tr)                                  # probes before exclusion
  betas_tr <- betas_tr[!(rownames(betas_tr) %in% cpg_excl), , drop = FALSE]
  cat("Exclusión de CpGs — sensibilidad", sens_lbl, ":", n0, "→", nrow(betas_tr), "sondas (−",
      n0 - nrow(betas_tr), "del panel presentes)\n")
}

eps    <- 1e-6                                          # keeps betas away from 0/1 before the logit
bc     <- pmin(pmax(betas_tr, eps), 1 - eps)
M_tr   <- log2(bc / (1 - bc))                           # training M-values (limma input)
rm(bc); gc()

cat("M_tr:", nrow(M_tr), "×", ncol(M_tr), "| NA:", sum(is.na(M_tr)),
    "| rango:", paste(round(range(M_tr), 2), collapse = " .. "), "\n")


# ---- 2c. Design matrix and rank check ------------------------------------------------------------
# ~ group + cohort + sex + age + cell fractions + smoking score (+ ancestry PCs) (+ chip row/column).
# The coefficient of interest is groupSCZ. The point-imputed age is used here only to build and check
# the design; step 3 refits on each of the 20 imputations.

pheno_tr$group <- factor(pheno_tr$group, levels = c("CNT", "SCZ"))     # CNT = reference
pheno_tr$gse   <- factor(pheno_tr$gse)                                  # cohort as fixed covariate
pheno_tr$sexo  <- factor(pheno_tr$sexo)
# smoking_SSc enters as numeric.

extra_anc <- if (length(ANC_COVS)) paste("+", paste(ANC_COVS, collapse = " + ")) else ""   # ancestry PCs
if (USAR_POSICION) extra_anc <- paste(extra_anc, "+ fila + columna")   # chip row and column
form   <- as.formula(paste("~ group + gse + sexo + age_epigenetic_imputed +",
                           paste(CELL_COVS, collapse = " + "), "+", COV_SMK, extra_anc))
design <- model.matrix(form, data = pheno_tr)

cat("Dimensiones de la matriz de diseño:", nrow(design), "×", ncol(design), "\n")
cat("Columnas:\n"); print(colnames(design))
cat("¿rango completo?:", qr(design)$rank == ncol(design),                # no aliasing among covariates
    "(rank", qr(design)$rank, "de", ncol(design), ")\n")
cat("¿muestras del diseño = columnas de M_tr?:", nrow(design) == ncol(M_tr), "\n")


# ---- 3. DMPs: limma per age imputation, combined by Rubin's rules --------------------------------

# ---- 3a. limma fit for each of the m imputations (resumable checkpoint) --------------------------
# Only the age column changes between fits. This is the expensive step.

m_imp  <- mids$m                                         # number of imputations
ck_dmp <- file.path(DIR_DMP, "dmp_perimput_checkpoint.rds")

if (file.exists(ck_dmp)) {
  ck <- readRDS(ck_dmp); coef_mat <- ck$coef; se_mat <- ck$se; df_total <- ck$df   # reuse saved fits
  cat("Per-imputación RECARGADO del checkpoint:", ck_dmp, "\n")
} else {
  coef_mat <- matrix(NA_real_, nrow(M_tr), m_imp, dimnames = list(rownames(M_tr), NULL))
  se_mat   <- matrix(NA_real_, nrow(M_tr), m_imp)
  df_total <- NULL
  for (i in seq_len(m_imp)) {                            # one age imputation per iteration
    pheno_tr$age_imp_i <- complete(mids, i)[rownames(pheno_tr), "age_chronological"]
    di  <- model.matrix(as.formula(paste("~ group + gse + sexo + age_imp_i +",
                          paste(CELL_COVS, collapse = " + "), "+", COV_SMK, extra_anc)),
                        data = pheno_tr)
    fit <- eBayes(lmFit(M_tr, di))                                      # empirical Bayes moderation
    coef_mat[, i] <- fit$coefficients[, "groupSCZ"]                     # logFC in M-values
    se_mat[, i]   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupSCZ"]  # moderated SE
    df_total      <- fit$df.total                                       # complete-data (moderated) df
    message("  imputación ", i, "/", m_imp, " OK")
  }
  saveRDS(list(coef = coef_mat, se = se_mat, df = df_total), ck_dmp)   # also read by 06c
  cat("Per-imputación calculado y guardado:", ck_dmp, "\n")
}
cat("coef_mat:", nrow(coef_mat), "×", ncol(coef_mat), "| NA:", sum(is.na(coef_mat)), "\n")


# ---- 3b. Rubin's rules per probe -----------------------------------------------------------------
# T = Wbar + (1 + 1/m) B; degrees of freedom capped at the complete-data df.

coef  <- rowMeans(coef_mat)                              # pooled estimate (logFC in M)
Wbar  <- rowMeans(se_mat^2)                              # within-imputation variance
B     <- matrixStats::rowVars(coef_mat)                  # between-imputation variance
Tvar  <- Wbar + (1 + 1 / m_imp) * B                      # Rubin total variance
SE    <- sqrt(Tvar)
tstat <- coef / SE

B_safe <- pmax(B, .Machine$double.eps)                   # avoids infinite df when B is about 0
r      <- (1 + 1 / m_imp) * B_safe / Wbar                # relative increase in variance
df_rub <- (m_imp - 1) * (1 + 1 / r)^2                    # Rubin df
df_use <- pmin(df_rub, df_total)                         # capped at the complete-data df
p      <- 2 * pt(-abs(tstat), df = df_use)

lambda_pre <- median(tstat^2) / qchisq(0.5, 1)           # genomic inflation before bacon
cat("Rubin OK | λ (inflación) cruda:", round(lambda_pre, 3), "\n")


# ---- 3c. Check: Rubin vs a single fit with the point-imputed age ---------------------------------
# Expected correlation close to 1: multiple imputation should not change the conclusions.

di0  <- model.matrix(as.formula(paste("~ group + gse + sexo + age_epigenetic_imputed +",
                       paste(CELL_COVS, collapse = " + "), "+", COV_SMK, extra_anc)),
                     data = pheno_tr)
fit0 <- eBayes(lmFit(M_tr, di0))                         # single fit with the mean imputed age
cat("cor(coef Rubin, coef único):", round(cor(coef, fit0$coefficients[, "groupSCZ"]), 4), "\n")
cat("cor(-log10 p Rubin, -log10 p único):",
    round(cor(-log10(p), -log10(fit0$p.value[, "groupSCZ"])), 4), "\n")


# ---- 4. bacon, FDR, Δβ and DMP table -------------------------------------------------------------

# ---- 4a. Inflation and bias correction with bacon ------------------------------------------------

bc      <- bacon(effectsizes = coef, standarderrors = SE)   # one set of pooled statistics
coef_bc <- as.numeric(es(bc)); se_bc <- as.numeric(se(bc)); p_bc <- as.numeric(pval(bc))   # corrected
z_bc    <- coef_bc / se_bc                                  # corrected z
lambda_post <- median(z_bc^2) / qchisq(0.5, 1)              # inflation after bacon
cat("bacon — inflación estimada:", round(inflation(bc), 3), "| sesgo:", round(bias(bc), 3), "\n")
cat("λ genómico:  antes =", round(lambda_pre, 3), " | después de bacon =", round(lambda_post, 3), "\n")


# ---- 4b. FDR, Δβ, group means, DMP call and high-confidence flag ---------------------------------
# Significance is FDR < 0.05 (Benjamini-Hochberg). es_dmp adds |adjusted Δβ| >= UMBRAL_DB_ADJ (0.5%
# in the main analysis, no threshold elsewhere). hiconf marks |Δβ| >= 0.02, above 450K technical noise;
# it flags probes and does not filter them.

DB_HICONF <- 0.02                                       # high-confidence flag threshold
fdr <- p.adjust(p_bc, method = "BH")                    # Benjamini-Hochberg on post-bacon p
scz <- pheno_tr$group == "SCZ"; cnt <- pheno_tr$group == "CNT"
m_scz <- rowMeans(betas_tr[, scz, drop = FALSE])              # mean beta per group
m_cnt <- rowMeans(betas_tr[, cnt, drop = FALSE])
dbeta <- m_scz - m_cnt                                        # unadjusted Δβ (difference of group means)

# Adjusted Δβ: the unadjusted Δβ ignores the design, so the adjusted effect is taken from the
# bacon-corrected coefficient. It is added to control M-values and the mean change in beta is taken,
# mean_j[ilogit2(M_cnt_j + coef_bc) - ilogit2(M_cnt_j)].
ilogit2 <- function(m) 2^m / (2^m + 1)                        # M-value to beta
M_cnt   <- M_tr[, cnt, drop = FALSE]                          # control M-values
delta_beta_adj <- rowMeans(ilogit2(M_cnt + coef_bc) - ilogit2(M_cnt))   # coef_bc recycled by row

# hg19 coordinates from the Illumina manifest; gene, promoter and island from Zhou 2017.
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
a    <- anno[rownames(M_tr), c("chr", "pos")]                  # position only (hg19)
zg   <- anotarGenesZhou(rownames(M_tr))                        # gene, promoter, island (GENCODE v36)

dmp <- data.frame(cpg = rownames(M_tr), chr = a$chr, pos = a$pos,
                  gene = zg$gene, is_promoter = zg$is_promoter, island = zg$island,
                  transcriptTypes = zg$transcriptTypes, dist_TSS = zg$dist_nearest_TSS,
                  logFC_M = coef_bc, SE = se_bc, t = z_bc, p = p_bc, FDR = fdr,
                  delta_beta = dbeta, delta_beta_adj = delta_beta_adj,   # unadjusted and adjusted Δβ
                  beta_SCZ = m_scz, beta_CNT = m_cnt,
                  hiconf = abs(dbeta) >= DB_HICONF,                     # flag, does not filter
                  es_dmp = fdr < 0.05 & abs(delta_beta_adj) >= UMBRAL_DB_ADJ,   # DMP call of this arm
                  row.names = NULL)
dmp <- dmp[order(dmp$FDR, dmp$p), ]                           # sorted by FDR

n_fdr <- sum(dmp$FDR < 0.05)                                  # CpGs at FDR < 0.05
n_dmp <- sum(dmp$es_dmp)                                      # DMPs (FDR and effect threshold)
n_hi  <- sum(dmp$es_dmp & dmp$hiconf)                         # DMPs with |Δβ| >= 0.02
cat("CpG con FDR<0.05:", n_fdr, "| DMP (FDR<0.05 y |Δβ aj|>=", UMBRAL_DB_ADJ, "):", n_dmp,
    "| de ellas alta confianza (|Δβ|>=", DB_HICONF, "):", n_hi, "\n")
cat("Top 5 por FDR:\n"); print(utils::head(dmp[, c("cpg","gene","logFC_M","delta_beta","hiconf","FDR")], 5))

# dmp_resultados.csv: one row per probe, post-bacon statistics, FDR, Δβ (unadjusted and adjusted),
# group means, annotation, hiconf and es_dmp.
readr::write_csv(dmp, file.path(DIR_DMP, "dmp_resultados.csv"))
cat("Guardado: dmp_resultados.csv —", nrow(dmp), "sondas\n")


# ---- 4c. QQ and volcano plots --------------------------------------------------------------------

options(repr.plot.width = 12, repr.plot.height = 5)
par(mfrow = c(1, 2))
obs <- -log10(sort(p_bc)); exp <- -log10(ppoints(length(p_bc)))   # observed vs expected -log10 p
plot(exp, obs, pch = 19, cex = 0.3, xlab = "-log10(p) esperado", ylab = "-log10(p) observado",
     main = paste0("QQ (λ=", round(lambda_post, 2), ")")); abline(0, 1, col = "red")
plot(dbeta, -log10(p_bc), pch = 19, cex = 0.3, xlab = "Δβ (SCZ − CNT)", ylab = "-log10(p)",
     main = "Volcano"); abline(v = c(-DB_HICONF, DB_HICONF), lty = 2, col = "grey50")
par(mfrow = c(1, 1))


# ---- 5. DMRs with dmrff --------------------------------------------------------------------------
# dmrff runs on the bacon-corrected statistics, with M-values to estimate correlation between nearby
# probes. DB_pass uses the same 0.02 threshold as hiconf. Regions are annotated with Zhou 2017 genes.

library(dmrff)

dmr <- dmrffDB(estimate = coef_bc, se = se_bc, p.value = p_bc,   # bacon-corrected statistics
               methylation = M_tr, chr = a$chr, pos = a$pos,     # M-values + hg19 coordinates
               dbeta = dbeta, maxgap = 500,                      # per-probe unadjusted Δβ
               p.cutoff = 0.05, dbeta.cutoff = DB_HICONF)        # DB_pass at |mean Δβ| >= 0.02

# Unique genes of the CpGs inside each region, and whether any of them is in a promoter.
reg_gene <- reg_prom <- character(nrow(dmr))
for (i in seq_len(nrow(dmr))) {
  en  <- a$chr == dmr$chr[i] & a$pos >= dmr$start[i] & a$pos <= dmr$end[i]   # CpGs of the region
  gi  <- unique(unlist(strsplit(zg$gene[en], ";", fixed = TRUE)))            # their genes
  gi  <- gi[!is.na(gi) & nzchar(gi)]
  reg_gene[i] <- paste(gi, collapse = ";")
  reg_prom[i] <- if (any(zg$is_promoter[en], na.rm = TRUE)) "yes" else "no"  # any promoter CpG
}
dmr$genes <- reg_gene
dmr$has_promoter <- reg_prom

n_sig <- sum(dmr$p.adjust < 0.05, na.rm = TRUE)                  # significant regions
n_db  <- sum(dmr$p.adjust < 0.05 & dmr$DB_pass, na.rm = TRUE)    # significant and |mean Δβ| >= 0.02
cat("DMR con p.adjust<0.05:", n_sig, "| de ellas |Mean Dbeta|≥", DB_HICONF, ":", n_db, "\n")
cat("Top 5 DMR por p.adjust:\n")
print(utils::head(dmr[, c("chr","start","end","n","meanActualDbeta","prop_concordante","p.adjust","genes")], 5))

readr::write_csv(dmr, file.path(DIR_DMP, "dmr_resultados.csv"))
cat("Guardado: dmr_resultados.csv —", nrow(dmr), "regiones\n")
