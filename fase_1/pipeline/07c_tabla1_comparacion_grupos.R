# 07c — Case-control comparison for Table 1 (p-values)
#
# Compares cases and controls of the training set, per cohort and overall, for the characteristics
# shown in Table 1: sex (Fisher's exact test), age, EpiSmokEr smoking score and EpiDISH neutrophil
# proportion (Welch's t test). Age is the recorded age or, where missing, the mean of the 20
# imputations, as it enters the model. The tests are descriptive: all four variables are covariates.
# Paper section: Covariates (Table 1)
# Input:  data_combat_P2/phenotype_final.csv (from 05); data_pheno/covariables_epigeneticas_P2.csv (04)
# Output: data_pheno/tabla1_valores_p.csv

suppressMessages({ library(data.table) })
# Allows running from the repository root as well as from fase_1/pipeline.
if (!file.exists("../data_combat_P2") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

ph  <- fread("../data_combat_P2/phenotype_final.csv")              # phenotype with partition and imputed age
cov <- fread("../data_pheno/covariables_epigeneticas_P2.csv")      # holds the continuous smoking score
ph$smoking_SSc <- cov$smoking_SSc[match(ph$gsm_id, cov$gsm_id)]   # join by GSM ID
tr  <- ph[set == "train"]                                          # training set (80%) only
cat("Entrenamiento:", nrow(tr), "| casos:", sum(tr$group == "SCZ"), "| controles:", sum(tr$group == "CNT"), "\n")

# ---- Tests per cohort and for all cohorts together ---------------------------------------------
res <- data.frame()
for (g in c(sort(unique(tr$gse)), "Todas")) {
  s <- if (g == "Todas") tr else tr[gse == g]                    # cohort subset, or all samples
  cas <- s$group == "SCZ"                                         # TRUE = case
  p_sexo <- fisher.test(table(s$sexo, s$group))$p.value           # sex: Fisher's exact test
  p_edad <- t.test(s$age_epigenetic_imputed[cas], s$age_epigenetic_imputed[!cas])$p.value   # age: Welch
  p_smk  <- t.test(s$smoking_SSc[cas], s$smoking_SSc[!cas])$p.value   # smoking score: Welch
  p_neu  <- t.test(s$Neutro[cas], s$Neutro[!cas])$p.value             # neutrophils: Welch
  res <- rbind(res, data.frame(gse = g, n_cnt = sum(!cas), n_scz = sum(cas), p_sexo = p_sexo, p_edad = p_edad,
                               p_tabaquismo = p_smk, p_neutrofilos = p_neu))
}
print(res, row.names = FALSE, digits = 3)
fwrite(res, "../data_pheno/tabla1_valores_p.csv")                  # one row per cohort plus "Todas"
cat("Guardado: ../data_pheno/tabla1_valores_p.csv\n")
