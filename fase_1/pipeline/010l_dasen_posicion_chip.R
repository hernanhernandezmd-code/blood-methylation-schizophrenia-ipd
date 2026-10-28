# 010l — dasen with and without the chip-position refinement of background adjustment
#
# In script 03, dasen (wateRmelon) runs with the signal columns named by GSM accession, so its optional
# refinement of the type I background offset by the sample's position on the chip (RxxCxx, internal
# function dfsfit) is skipped ("Sentrix position model failed, skipping"). No change is expected:
# dfsfit subtracts one offset per sample from type I probes, position only smooths that offset with a
# row + column model, and the quantile normalization that follows works on within-sample ranks, so a
# per-sample constant cancels out. The script tests this in one cohort, GSE80417 (675 samples), with the
# P2 mask and the same probes as script 03: dasen is run twice on the same intensities, without position
# (as in 03) and with the real position (roco), and the betas and the case-control t statistics of both
# versions are compared.
# Paper section: Robustness and Sensitivity Analyses (Supplementary Table S1)
# Input:  data_curated/GSE80417_curated_P2.rds (from 02); data_pheno/phenodata_combinado_imputar_edad.csv;
#         data_dmp_P2_anc/dmp_resultados.csv (probe panel)
# Output: data_dmp_P2_anc_dasen_posicion/dasen_posicion_resumen.csv
# GSE: cohort used for the test (GSE80417)

# ---- 1. Packages, data and chip position of each sample ------------------------------------------
suppressMessages({ library(wateRmelon); library(limma); library(data.table)
                   library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
# works whether the session starts in fase_1/pipeline or at the repository root
if (!file.exists("../data_curated") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")
GSE     <- "GSE80417"                                                # cohort tested
DIR_OUT <- "../data_dmp_P2_anc_dasen_posicion"                       # output folder
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

x <- readRDS(file.path("../data_curated", paste0(GSE, "_curated_P2.rds")))   # $M and $U, probes x samples
M <- x$M; U <- x$U; rm(x); gc()                                 # methylated and unmethylated signal
ph <- as.data.frame(fread("../data_pheno/phenodata_combinado_imputar_edad.csv"))   # carries chip position
ph <- ph[ph$gse == GSE, ]                                       # samples of the tested cohort
cat("Columnas de la señal (primeras):", head(colnames(M), 3), "\n")
# columns may carry the Sentrix code or the GSM; they are renamed to GSM as in script 03
if (!all(colnames(M) %in% ph$gsm_id)) {
  i <- match(colnames(M), ph$sentrix_position)
  stopifnot(!anyNA(i))                                          # every column must map to a sample
  colnames(M) <- colnames(U) <- ph$gsm_id[i]
}
ph <- ph[match(colnames(M), ph$gsm_id), ]                       # phenotype in signal column order
stopifnot(identical(ph$gsm_id, colnames(M)), all(grepl("^R0[1-6]C0[1-2]$", ph$position)))   # 6 rows x 2 cols
# same probes that script 03 normalized: the panel common to the six cohorts (397,194 CpGs)
panel <- fread("../data_dmp_P2_anc/dmp_resultados.csv", select = "cpg")$cpg
M <- M[rownames(M) %in% panel, ]; U <- U[rownames(M), ]        # M and U restricted to the panel
anno   <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)   # Illumina 450K manifest
onetwo <- anno[rownames(M), "Type"]; names(onetwo) <- rownames(M)   # probe design type, needed by dasen
stopifnot(!anyNA(onetwo))
cat(GSE, ":", nrow(M), "sondas ×", ncol(M), "muestras | tipo I:", sum(onetwo == "I"),
    "| tipo II:", sum(onetwo == "II"), "\n")
print(table(ph$group))                                          # cases and controls in the cohort

# ---- 2. dasen without and with chip position -----------------------------------------------------
b_sin <- dasen(M, U, onetwo)                                  # as in script 03: GSM columns, position unused
b_con <- dasen(M, U, onetwo, roco = ph$position)              # with each sample's real row and column
stopifnot(identical(dim(b_sin), dim(b_con)))
dif <- abs(b_con - b_sin)                                     # absolute beta difference per probe and sample
cat("Diferencia máxima absoluta entre betas:", format(max(dif, na.rm = TRUE), digits = 3), "\n")
cat("Proporción de valores idénticos:", mean(dif == 0, na.rm = TRUE), "\n")

# ---- 3. Functional check: case-control contrast with each version --------------------------------
# A diagnosis-only model on M-values; identical betas give identical t statistics.
grupo <- factor(ph$group, levels = c("CNT", "SCZ"))           # controls as reference
dis   <- model.matrix(~ grupo)                                # diagnosis only
t_de  <- list()
for (version in c("sin", "con")) {
  b <- if (version == "sin") b_sin else b_con                # betas of this version
  b <- pmin(pmax(b, 1e-6), 1 - 1e-6)                          # keep betas inside (0, 1) for the logit
  fit <- eBayes(lmFit(log2(b / (1 - b)), dis))                # M-values, moderated t
  t_de[[version]] <- fit$t[, "grupoSCZ"]                    # SCZ vs CNT statistic per probe
}
cor_t <- cor(t_de$sin, t_de$con, use = "complete.obs")       # agreement of the case-control statistics
cat("Correlación de los t (sin frente a con posición):", round(cor_t, 6), "\n")

# ---- 4. Summary and output -----------------------------------------------------------------------
# beta agreement and case-control statistic agreement in one row
resumen <- data.frame(cohorte = GSE, sondas = nrow(M), muestras = ncol(M),
                      dif_max_beta = max(dif, na.rm = TRUE), prop_identicos = mean(dif == 0, na.rm = TRUE),
                      cor_t_caso_control = cor_t)
print(resumen, row.names = FALSE)
fwrite(resumen, file.path(DIR_OUT, "dasen_posicion_resumen.csv"))
cat("Guardado en", DIR_OUT, "\n")
