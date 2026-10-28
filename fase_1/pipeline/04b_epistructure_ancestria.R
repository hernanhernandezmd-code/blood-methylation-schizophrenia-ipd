# 04b — Ancestry components (EPISTRUCTURE)
#
# Recovers population structure without genotypes by PCA on a panel of ancestry-informative CpGs
# (EPISTRUCTURE, Rahmani et al., 2017) and adds the leading PCs to the covariate table. The PCs are
# computed on the harmonized betas, as registered, and checked to track cohort and not diagnosis.
# Paper section: Covariates
# Input:  data_harmonized_P2/betas_harmonized.rds, refs/epistructure_ancestry_cpgs_450k.txt,
#         data_pheno/phenodata_combinado_imputar_edad.csv (diagnostic group)
# Output: ancestry_PC1..K appended to data_pheno/covariables_epigeneticas_P2.csv;
#         epistructure_pca.rds in the harmonized folder
# VARIANTE_MASKING: "P1" | "P2" | "P3" (default "P2"); selects input and output folders.
# EXCLUIR_COHORTE: "" (default, all cohorts) or a GSE ID to drop; adds a suffix to folders.
# Run after 04: 04 rewrites the covariate table and would drop the PCs added here.

# ---- 1. Packages, paths and data --------------------------------------------
suppressMessages({ library(data.table); library(readr) })

VARIANTE_MASKING <- "P2"                                # probe-masking variant (P2 = registered panel)
SUF_VAR <- if (VARIANTE_MASKING == "P1") "" else paste0("_", VARIANTE_MASKING)  # folder suffix
EXCLUIR_COHORTE <- ""                                   # "" keeps all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""  # suffix only when a cohort is dropped
DIR_HARMONIZED <- paste0("../data_harmonized", SUF_VAR, SUF_COH)  # regRCPqn betas from 04
DIR_PHENO      <- "../data_pheno"                       # phenotype and covariate tables
DIR_REFS       <- "../refs"                             # cached reference CpG lists

harm  <- readRDS(file.path(DIR_HARMONIZED, "betas_harmonized.rds"))   # list(betas, cohorte)
betas <- harm$betas; cohorte <- harm$cohorte; rm(harm); gc()          # CpGs x samples; free the list
cpgs_anc <- readLines(file.path(DIR_REFS, "epistructure_ancestry_cpgs_450k.txt"))  # EPISTRUCTURE panel
cat("Betas:", nrow(betas), "×", ncol(betas), "| CpGs ancestría (lista):", length(cpgs_anc), "\n")

# ---- 2. PCA on the ancestry panel -------------------------------------------
en_panel <- intersect(cpgs_anc, rownames(betas))        # ancestry CpGs that survived QC and masking
cat("CpGs de ancestría presentes en el panel P2:", length(en_panel),
    sprintf("(%.1f%% de %d)\n", 100 * length(en_panel) / length(cpgs_anc), length(cpgs_anc)))

X <- t(betas[en_panel, , drop = FALSE])                 # samples x CpGs
X <- X[, apply(X, 2, sd) > 0, drop = FALSE]             # zero-variance CpGs would give NaN when scaled
pca <- prcomp(X, center = TRUE, scale. = TRUE)          # standardized CpGs, as in EPISTRUCTURE
ve  <- pca$sdev^2 / sum(pca$sdev^2)                     # proportion of variance per PC
cat("Varianza explicada primeros 10 PCs:\n"); print(round(ve[1:10], 4))

# K is the scree elbow: largest second difference of explained variance, minimum 2.
dv  <- -diff(ve[1:11])                                  # drop in variance from one PC to the next
K   <- max(2, which.max(-diff(dv[1:9])) + 1)            # where that drop flattens most
cat("\nK sugerido por scree (codo):", K, "| varianza acumulada a K:", round(sum(ve[1:K]), 3), "\n")

options(repr.plot.width = 6, repr.plot.height = 4)      # plot size in the notebook
plot(ve[1:10], type = "b", pch = 19, xlab = "PC", ylab = "Varianza explicada",
     main = "Scree — PCs de ancestría (EPISTRUCTURE)"); abline(v = K, lty = 2, col = "red")  # K marked

# ---- 3. Check that the PCs track cohort and not diagnosis -------------------
# R-squared of each PC on cohort (expected to be appreciable) and on SCZ/CNT group (expected ~0).
pcs <- as.data.frame(pca$x[, 1:K, drop = FALSE]); colnames(pcs) <- paste0("ancestry_PC", 1:K)  # scores
pheno_min <- data.frame(gsm_id = rownames(pcs), cohorte = cohorte[rownames(pcs)], pcs)  # add cohort
grp_tab <- as.data.frame(readr::read_csv(file.path(DIR_PHENO, "phenodata_combinado_imputar_edad.csv"),
                                         show_col_types = FALSE))       # inclusive phenotype table
pheno_min$group <- grp_tab$group[match(pheno_min$gsm_id, grp_tab$gsm_id)]  # SCZ/CNT by GSM ID

cat("Asociación de cada PC con COHORTE (R²) y con GRUPO (R²):\n")
for (k in 1:K) {
  pc <- pheno_min[[paste0("ancestry_PC", k)]]
  r2_coh <- summary(lm(pc ~ factor(pheno_min$cohorte)))$r.squared  # variance explained by cohort
  r2_grp <- summary(lm(pc ~ factor(pheno_min$group)))$r.squared    # by diagnosis; ~0 if not confounded
  cat(sprintf("  PC%d: R²(cohorte)=%.3f  |  R²(grupo)=%.4f\n", k, r2_coh, r2_grp))
}
cat("\n(esperado: R²(cohorte) apreciable = captura ancestría; R²(grupo) ~0 = no confundido con SCZ)\n")

# ---- 4. Append PCs to the covariate table -----------------------------------
cov_file <- file.path(DIR_PHENO, paste0("covariables_epigeneticas", SUF_VAR, SUF_COH, ".csv"))  # from 04
cov <- as.data.frame(readr::read_csv(cov_file, show_col_types = FALSE))
cov <- cov[, !grepl("^ancestry_PC", colnames(cov)), drop = FALSE]   # drop PCs left by an earlier run of 04b
pcs$gsm_id <- rownames(pcs)                             # join key
cov <- merge(cov, pcs, by = "gsm_id", all.x = TRUE, sort = FALSE)   # keep every covariate row
cov <- cov[match(colnames(betas), cov$gsm_id), ]        # same sample order as the betas
readr::write_csv(cov, cov_file)                         # overwrite in place; 06 reads the PCs from here
cat("Guardado:", basename(cov_file), "— añadidas", K, "PCs de ancestría | filas:", nrow(cov), "\n")
saveRDS(list(pca = pca, ve = ve, K = K, cpgs = en_panel),
        file.path(DIR_HARMONIZED, "epistructure_pca.rds"))          # PCA object, for the scree figure
cat("Guardado: epistructure_pca.rds (para el scree y trazabilidad)\n")
