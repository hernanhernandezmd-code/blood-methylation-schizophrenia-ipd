# 12 — Figures of the article
#
# Draws the main and supplementary figures of the Phase 1 article from the result tables of the
# six-cohort analysis. The main analysis shown in Figures 1 and 2 is the model without smoking-associated
# CpGs and with chip row/column covariates (_noSmk_pos); robustness (Figure 3) is measured on the DMPs of
# the pre-registered primary model. Style: Okabe-Ito palette, 9 pt theme_bw, dense point layers rasterised.
# Figures written (file name = number in the article; the internal object names differ):
#   Figure 1  Figura_1_EWAS_principal     main analysis: Manhattan, QQ before/after bacon, volcano
#   Figure 2  Figura_2_reproducibilidad   directional concordance and per-cohort forest of 8 loci
#   Figure 3  Figura_3_robustez           leave-one-cohort-out, technical risk and probe-panel removal
#   Figure S2 Figura_S2_celularidad       estimated leukocyte fractions by cohort and group
#   Figure S3 Figura_S3_ancestria         ancestry principal components by cohort and group
#   Figure S4 Figura_S4_tabaco_primario   smoking-associated CpGs among the primary-model DMPs
#   Figure S5 Figura_S5_contexto          genomic context of the main-analysis DMPs
#   Figure S6 Figura_S6_forest_menor_FDR  per-cohort forest of the 4 lowest-FDR loci
#   (Figure S1, the PRISMA-IPD flow diagram, is built in Word from the official template, not here.)
# Paper section: Figures
# Input:  ../data_dmp_P2_anc (primary), ../data_dmp_P2_anc_noSmk_pos (main), the _noSmk, _noAnc, _noAge,
#         _LOCO, _riesgo_tecnico and _concordancia result folders, ../data_harmonized_P2/betas_harmonized.rds,
#         ../data_combat_P2/phenotype_final.csv, ../data_pheno (covariates, chip position),
#         ../refs/smoking_cpgs_exclusion_450k.txt
# Output: ../figuras_SchizBull/ (vector PDF, 300 dpi TIFF and PNG per figure)
# Widths follow the journal: 88 mm (one column) or 180 mm (two columns).

# ---- 1. Packages, paths and common style --------------------------------------------------------
suppressMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(data.table); library(readr)
  library(ggrastr)          # rasterises dense point layers, axes stay vector
})

# Run from the repository root or from fase_1/pipeline: all paths are relative to the latter.
if (!file.exists("../data_dmp_P2_anc") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

DIR_PRIM  <- "../data_dmp_P2_anc"                 # pre-registered primary model, all probes
DIR_MAIN  <- "../data_dmp_P2_anc_noSmk_pos"       # main analysis: no smoking CpGs + chip row/column
DIR_LOCO  <- "../data_dmp_P2_anc_LOCO"            # leave-one-cohort-out refits (010_loco)
DIR_RT    <- "../data_dmp_P2_anc_riesgo_tecnico"  # technical-risk refits (010b)
DIR_CONC  <- "../data_dmp_P2_anc_concordancia"    # directional concordance per DMP (010c)
DIR_NOSMK <- "../data_dmp_P2_anc_noSmk"           # primary model without smoking CpGs, no chip position
DIR_NOANC <- "../data_dmp_P2_anc_noAnc"           # primary model without ancestry CpGs
DIR_NOAGE <- "../data_dmp_P2_anc_noAge"           # primary model without epigenetic-clock CpGs
DIR_HARM  <- "../data_harmonized_P2"              # harmonised betas (regRCPqn), for the forest plots
DIR_PHENO <- "../data_combat_P2"                  # phenotype with partition and imputed age
F_COV     <- "../data_pheno/covariables_epigeneticas_P2.csv"   # smoking score and ancestry PCs
F_TABACO  <- "../refs/smoking_cpgs_exclusion_450k.txt"         # smoking CpG panel (Joehanes + Bollepalli)
DIR_FIG   <- "../figuras_SchizBull"   # figure output folder
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

FDR_SIG   <- 0.05       # significance threshold
DB_HICONF <- 0.02       # effect flag above array technical noise (|delta beta| >= 2%)
DUBLIN    <- "GSE147221"                                      # cohort singled out in several panels
CELL_COVS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")   # EpiDISH leukocyte fractions

# Okabe-Ito palette (colour-blind safe). Grey first: it is the background colour of every panel.
COL <- c(gris = "#999999", azul = "#0072B2", naranja = "#E69F00", verde = "#009E73",
         rojo = "#D55E00", morado = "#CC79A7", celeste = "#56B4E9", amarillo = "#F0E442")
# Cohorts sharing the sign of the pooled effect: blue (all six) to red (two or fewer).
# concordance classes, colours in legend order
COL_CONC <- c("6/6" = unname(COL["azul"]), "5/6" = unname(COL["celeste"]), "4/6" = "grey75",
              "3/6" = unname(COL["naranja"]), "≤2/6" = unname(COL["rojo"]))

# Common theme: thin marks, recessive grid, no unnecessary boxes.
tema_pub <- theme_bw(base_size = 9) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major   = element_line(linewidth = 0.25, colour = "grey88"),   # light major grid
        panel.border       = element_rect(linewidth = 0.4, colour = "grey40"),   # thin grey panel frame
        axis.ticks         = element_line(linewidth = 0.3, colour = "grey40"),
        strip.background   = element_rect(fill = "grey92", colour = NA),   # facet header without border
        strip.text         = element_text(face = "bold", size = 8),
        plot.title         = element_text(face = "bold", size = 9),
        plot.subtitle      = element_text(size = 8, colour = "grey30"),
        plot.tag           = element_text(face = "bold", size = 11),   # panel letters A, B, C
        legend.key.size    = unit(3.2, "mm"),   # compact legend keys
        legend.title       = element_text(size = 8),
        legend.text        = element_text(size = 7.5))
theme_set(tema_pub)     # applies to every plot below

# Saves each figure as vector PDF, 300 dpi TIFF (journal format) and PNG for review.
guardar_fig <- function(plot, nombre, ancho_mm, alto_mm) {
  # cairo_pdf embeds the Unicode symbols (delta, beta, minus) used in the axis labels
  ggsave(file.path(DIR_FIG, paste0(nombre, ".pdf")), plot,
         width = ancho_mm, height = alto_mm, units = "mm", device = cairo_pdf, bg = "white")
  # LZW-compressed TIFF, lossless
  ggsave(file.path(DIR_FIG, paste0(nombre, ".tiff")), plot,
         width = ancho_mm, height = alto_mm, units = "mm", dpi = 300, compression = "lzw", bg = "white")
  ggsave(file.path(DIR_FIG, paste0(nombre, ".png")), plot,
         width = ancho_mm, height = alto_mm, units = "mm", dpi = 300, bg = "white")
  cat("  guardada:", nombre, sprintf("(%d x %d mm)\n", ancho_mm, alto_mm))
}

# Gene name used as label: the first one that is not an unnamed clone or lncRNA (AC012345.1, LINC, -AS1).
nombre_gen <- function(g) {
  out <- rep(NA_character_, length(g))   # NA when no usable name
  for (i in seq_along(g)) {
    if (is.na(g[i]) || g[i] == "") next                 # intergenic: no label
    for (x in strsplit(g[i], ";")[[1]]) {               # Zhou annotation lists genes separated by ";"
      if (!grepl("^(AC|AL|AP|BX|CT|FP|Z)[0-9]{5,}|^LINC|-AS[0-9]", x)) { out[i] <- x; break }
    }
  }
  out
}

cat("Rutas listas. Salida:", normalizePath(DIR_FIG), "\n")

# ---- 2. Load data ---------------------------------------------------------------------------------
# Everything comes from precomputed tables, except the forest plots, which read the betas of 12 CpGs
# from the harmonised matrix in section 6.
prim <- as.data.frame(fread(file.path(DIR_PRIM, "dmp_resultados.csv")))        # primary model, all probes
main <- as.data.frame(fread(file.path(DIR_MAIN, "dmp_resultados.csv")))        # main analysis, all probes
conc <- as.data.frame(fread(file.path(DIR_CONC, "concordancia_por_DMP_noSmk.csv")))   # sign count per DMP
conc_res <- as.data.frame(fread(file.path(DIR_CONC, "concordancia_resumen.csv")))     # distribution per arm
loco <- as.data.frame(fread(file.path(DIR_LOCO, "loco_resumen.csv")))          # one row per left-out cohort
rt   <- as.data.frame(fread(file.path(DIR_RT, "riesgo_tecnico_resumen.csv")))  # technical-risk scenarios
ph   <- as.data.frame(fread(file.path(DIR_PHENO, "phenotype_final.csv")))      # samples, set, covariates
cpg_tabaco <- readLines(F_TABACO)                                              # smoking CpG IDs

# Main-analysis DMP = FDR < 0.05 and |adjusted delta beta| >= 0.5% (column es_dmp written by 06).
# Tables without that column (pre-registered primary model) fall back to FDR < 0.05 alone.
# logical vector of DMPs
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
sig_prim <- prim[prim$FDR < FDR_SIG, ]                   # primary DMPs (pre-registered rule)
sig_main <- main[sel_dmp(main), ]                        # main DMPs: FDR and adjusted effect threshold
sig_main <- sig_main[order(sig_main$FDR), ]              # ranked by FDR, used for the rank bins of Figure 2B
i <- match(sig_main$cpg, conc$cpg)   # row of each DMP in the concordance table
stopifnot(!anyNA(i))                                     # every main DMP must have a concordance count
sig_main$conc_aj  <- conc$concordancia_ajustada[i]       # cohorts sharing the sign, adjusted effect
sig_main$conc_cru <- conc$concordancia_cruda[i]          # same, crude delta beta
sig_main$nombre   <- nombre_gen(sig_main$gene)           # label for plots

cat("primario: ", nrow(prim), "CpG |", nrow(sig_prim), "DMP\n")
cat("principal:", nrow(main), "CpG |", nrow(sig_main), "DMP |",
    sum(abs(sig_main$delta_beta) >= DB_HICONF), "con |Δβ| >= 2 %\n")
cat("cohortes:", paste(sort(unique(ph$gse)), collapse = ", "), "| entrenamiento:", sum(ph$set == "train"), "\n")
stopifnot(all(c("signo_igual_pct", "pendiente") %in% names(loco)))   # LOCO table carries the 010b metrics

# ---- 3. Loci shown in the forest plots and labelled in Figure 1 -----------------------------------
# Two fixed rules, so the choice is not made by hand (one locus per gene in both):
# concordant = the 8 lowest-FDR DMPs with the same adjusted sign in all six cohorts and |delta beta| >= 2%
#   (forest in Figure 2C); lowest FDR = the 4 lowest-FDR DMPs with a named gene (forest in Figure S6).
# Figure 1 labels the 4 lowest-FDR loci and the first 4 concordant ones.
con_nombre <- sig_main[!is.na(sig_main$nombre), ]                     # DMPs with a usable gene name
top_fdr <- head(con_nombre[!duplicated(con_nombre$nombre), ], 4)      # already sorted by FDR
cand <- con_nombre[con_nombre$conc_aj == 6 & abs(con_nombre$delta_beta) >= DB_HICONF, ]   # 6/6 and >= 2%
top_conc <- head(cand[!duplicated(cand$nombre), ], 8)   # first 8 genes meeting the rule
stopifnot(nrow(top_fdr) == 4, nrow(top_conc) == 8, !any(top_fdr$cpg %in% top_conc$cpg))  # disjoint sets
top_fdr$grupo  <- "Lowest FDR"
top_conc$grupo <- "Concordant in all six cohorts"
LOCI <- rbind(top_fdr, top_conc)                        # the 12 loci whose betas are read in section 6
LOCI_F1 <- rbind(top_fdr, head(top_conc, 4))           # loci labelled in Figure 1
print(LOCI[, c("grupo", "cpg", "nombre", "delta_beta", "FDR", "conc_aj", "conc_cru")], row.names = FALSE)

# ---- 4. Figure 1: main analysis -------------------------------------------------------------------
# A, Manhattan of all tested positions; B, QQ before and after bacon (pre-bacon p-values from 06c, same
# Rubin's rules); C, volcano with the crude delta beta. The 8 labelled loci: orange circle = lowest FDR,
# blue triangle = concordant in 6/6. All positions are drawn; the dense layer is rasterised at 600 dpi.

# --- A: Manhattan ----------------------------------------------------------------------------------
mh <- main[, c("cpg", "chr", "pos", "p", "FDR", "delta_beta")]   # columns needed for the Manhattan
mh$chr_n <- suppressWarnings(as.integer(gsub("chr", "", mh$chr)))   # numeric chromosome
mh <- mh[!is.na(mh$chr_n), ]                                        # drops any non-autosomal label
mh <- mh[order(mh$chr_n, mh$pos), ]   # genome order
# Cumulative genomic position: each chromosome starts after the previous one plus a 20 Mb gap.
offset <- 0; centros <- numeric(0); mh$pos_acum <- NA_real_   # running offset, label centres, x coordinate
for (k in sort(unique(mh$chr_n))) {
  idx <- mh$chr_n == k   # probes of chromosome k
  mh$pos_acum[idx] <- mh$pos[idx] + offset   # shift by the length of the chromosomes before
  centros <- c(centros, offset + (max(mh$pos[idx]) + min(mh$pos[idx])) / 2)   # x position of the label
  offset  <- offset + max(mh$pos[idx]) + 2e7
}
names(centros) <- sort(unique(mh$chr_n))   # axis labels = chromosome numbers
mh$par <- factor(mh$chr_n %% 2)                          # alternating grey shades by chromosome
umbral_fdr <- max(main$p[main$FDR < FDR_SIG])            # largest p-value with FDR < 0.05 (dashed line)

etq <- merge(mh, LOCI_F1[, c("cpg", "nombre", "grupo")], by = "cpg")   # labelled points
etq$grupo <- factor(etq$grupo, levels = c("Lowest FDR", "Concordant in all six cohorts"))

pA1 <- ggplot(mh, aes(x = pos_acum, y = -log10(p))) +
  # ~390,000 points rasterised, the rest of the panel stays vector
  # background points, two greys
  rasterise(geom_point(aes(colour = par), size = 0.35, alpha = 0.7, show.legend = FALSE), dpi = 600) +
  # FDR = 0.05 line
  geom_hline(yintercept = -log10(umbral_fdr), linetype = "dashed", linewidth = 0.35, colour = COL["rojo"]) +
  # labelled loci drawn on top, filled shape by selection rule
  geom_point(data = etq, aes(shape = grupo, fill = grupo), size = 1.8, colour = "grey20", stroke = 0.35) +
  geom_text_repel(data = etq, aes(label = nombre), size = 2.2, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.25, max.overlaps = 20, seed = 1) +
  scale_colour_manual(values = c("grey65", "grey35")) +   # alternating chromosome greys
  scale_shape_manual(values = c(21, 24), name = NULL) +                     # circle, triangle
  scale_fill_manual(values = unname(COL[c("naranja", "azul")]), name = NULL) +
  # one tick per chromosome
  scale_x_continuous(breaks = centros, labels = names(centros), expand = expansion(mult = 0.01)) +
  labs(tag = "A", x = "Chromosome", y = expression(-log[10](italic(p)))) +
  # legend inside the top-left corner on a translucent box
  theme(panel.grid.major.x = element_blank(), legend.position = c(0.01, 0.98),
        legend.justification = c(0, 1), legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# --- B: QQ before and after bacon ------------------------------------------------------------------
pre <- read.csv(file.path(DIR_MAIN, "dmp_prebacon.csv"), stringsAsFactors = FALSE)   # written by 06c
stopifnot(identical(pre$cpg, main$cpg))                                  # same probes, same order
# Expected vs observed -log10 p; below -log10 p = 2 only every 20th point is kept (visually identical).
curva_qq <- function(pv, etiqueta) {
  obs <- sort(pv)
  # uniform quantiles as expected values
  d   <- data.frame(esp = -log10(ppoints(length(obs))), obs = -log10(obs), serie = etiqueta)
  rbind(d[d$obs > 2, ], d[d$obs <= 2, ][seq(1, sum(d$obs <= 2), by = 20), ])   # thins the dense region
}
# Genomic inflation factor: median chi-square (1 df) over its expected median.
lambda_de <- function(pv) median(qchisq(1 - pv, 1)) / qchisq(0.5, 1)
lam_pre  <- lambda_de(pre$p_prebacon)   # lambda of the Rubin-combined statistics
lam_post <- lambda_de(main$p)   # lambda after bacon
qq <- rbind(curva_qq(pre$p_prebacon, "Before bacon"), curva_qq(main$p, "After bacon"))
qq$serie <- factor(qq$serie, levels = c("Before bacon", "After bacon"))

pB1 <- ggplot(qq, aes(x = esp, y = obs, colour = serie)) +
  geom_abline(slope = 1, intercept = 0, colour = COL["rojo"], linewidth = 0.35) +   # null diagonal
  rasterise(geom_point(size = 0.5), dpi = 600) +   # QQ points rasterised
  scale_colour_manual(values = c("Before bacon" = "grey65", "After bacon" = "grey20"), name = NULL) +
  # lambda before and after, written as plotmath subscripts
  annotate("text", x = 0.2, y = max(qq$obs) * 0.92, hjust = 0, size = 2.6,
           label = sprintf('lambda[before] == "%.3f"', lam_pre), parse = TRUE) +
  annotate("text", x = 0.2, y = max(qq$obs) * 0.80, hjust = 0, size = 2.6,
           label = sprintf('lambda[after] == "%.3f"', lam_post), parse = TRUE) +
  labs(tag = "B", x = expression("Expected " * -log[10](italic(p))),
       y = expression("Observed " * -log[10](italic(p)))) +
  theme(legend.position = c(0.03, 0.62), legend.justification = c(0, 1), legend.text = element_text(size = 6.5),
        legend.background = element_blank())

# --- C: volcano ------------------------------------------------------------------------------------
# Four classes: not significant; FDR < 0.05 but below the adjusted effect threshold; DMP; DMP with
# crude |delta beta| >= 2%.
vol <- main[, c("cpg", "delta_beta", "FDR")]; vol$dmp <- sel_dmp(main)   # crude delta beta, FDR and DMP flag
vol$clase <- "Not significant"   # default class
vol$clase[vol$FDR < FDR_SIG & !vol$dmp] <- "FDR < 0.05, below effect threshold"
vol$clase[vol$dmp] <- "DMP"   # passes FDR and the adjusted effect threshold
vol$clase[vol$dmp & abs(vol$delta_beta) >= DB_HICONF] <- "DMP and |Δβ| ≥ 2%"   # and a large crude effect
vol$clase <- factor(vol$clase, levels = c("Not significant", "FDR < 0.05, below effect threshold", "DMP", "DMP and |Δβ| ≥ 2%"))
etq_v <- merge(vol, LOCI_F1[, c("cpg", "nombre")], by = "cpg")   # same 8 labelled loci as panel A

pC1 <- ggplot(vol, aes(x = 100 * delta_beta, y = -log10(FDR))) +
  # dotted lines at +/- 2% delta beta, dashed line at FDR = 0.05
  geom_vline(xintercept = c(-100, 100) * DB_HICONF, linetype = "dotted", linewidth = 0.3, colour = "grey50") +
  geom_hline(yintercept = -log10(FDR_SIG), linetype = "dashed", linewidth = 0.35, colour = COL["rojo"]) +
  # colour, shape and size all encode the class, so the larger-effect DMPs stand out
  rasterise(geom_point(aes(colour = clase, shape = clase, size = clase), alpha = 0.8), dpi = 600) +
  geom_text_repel(data = etq_v, aes(label = nombre), size = 2.2, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.25, max.overlaps = 30, box.padding = 0.45, seed = 1) +
  # grey, grey, blue, orange by class
  scale_colour_manual(values = c("grey70", "grey45", unname(COL["azul"]), unname(COL["naranja"])), name = NULL) +
  scale_shape_manual(values = c(16, 16, 16, 17), name = NULL) +   # triangle for DMP with |delta beta| >= 2%
  scale_size_manual(values = c(0.35, 0.5, 0.9, 1.4), name = NULL) +   # larger marks for the larger effects
  labs(tag = "C", x = expression("Methylation difference, cases − controls (" * Delta * beta * ", %)"),
       y = expression(-log[10]("FDR"))) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),   # legend inside the top-left corner
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA))

# Layout: Manhattan across the top, QQ and volcano side by side below.
fig1 <- pA1 / (pB1 | pC1) + plot_layout(heights = c(1, 1.05))
guardar_fig(fig1, "Figura_1_EWAS_principal", 180, 150)   # two-column width
cat(sprintf("  λ antes de bacon %.3f | después %.3f | DMP hipometiladas %.1f %%\n",
            lam_pre, lam_post, 100 * mean(sig_main$delta_beta < 0)))
fig1

# ---- 5. Figure 3: robustness, measured on the primary DMPs -----------------------------------------
# One row per refit: the six leave-one-cohort-out runs (010_loco), the technical-risk scenario that
# removes GSE147221 and GSE152027 (010b; removing GSE147221 alone is the matching LOCO row) and the three
# probe-panel removals (smoking, ancestry, clocks). A, DMPs per refit (dashed line = primary); B, percent of
# primary DMPs whose adjusted effect keeps its sign (dashed line = 50%, chance); C, slope of the refit
# effects on the primary effects (1 = same size).
rob <- data.frame()   # one row per refit
for (k in seq_len(nrow(loco))) {   # leave-one-cohort-out rows
  if (loco$escenario[k] == "sin_116378_y_152027") next          # extra scenario, not described in Methods
  rob <- rbind(rob, data.frame(bloque = "Leave one cohort out", etiqueta = sub("^sin_", "", loco$escenario[k]),
                               n_DMP = loco$n_DMP[k], signo = loco$signo_igual_pct[k],
                               pendiente = loco$pendiente[k], dublin = loco$escenario[k] == paste0("sin_", DUBLIN)))
}
# Technical-risk row: primary model refitted without GSE147221 and GSE152027.
k <- which(rt$brazo == "primario" & rt$escenario == "sin_Dublin_IoPPN")   # row of the technical-risk table
rob <- rbind(rob, data.frame(bloque = "Technical risk", etiqueta = "GSE147221 + GSE152027",
                             n_DMP = rt$DMP_escenario[k], signo = rt$signo_igual_pct[k],
                             pendiente = rt$pendiente[k], dublin = TRUE))

# Panel-removal rows: the three are the same refit (primary model without one CpG panel), so the
# smoking row reads _noSmk, not the main analysis with chip position.
PANELES <- list(list(dir = DIR_NOSMK, nombre = "Smoking CpGs"), list(dir = DIR_NOANC, nombre = "Ancestry CpGs"),
                list(dir = DIR_NOAGE, nombre = "Clock CpGs"))
for (pl in PANELES) {
  s  <- as.data.frame(fread(file.path(pl$dir, "dmp_resultados.csv")))   # refit results, all retained probes
  # primary DMPs still present in the refit, adjusted effects (M scale) side by side
  me <- merge(sig_prim[, c("cpg", "logFC_M")], s[, c("cpg", "logFC_M")], by = "cpg", suffixes = c("_ref", ""))
  quitadas <- nrow(prim) - nrow(s)                       # probes removed with the panel, shown in the label
  rob <- rbind(rob, data.frame(bloque = "Probe panel removed",
                               etiqueta = sprintf("%s (%s)", pl$nombre, format(quitadas, big.mark = ",")),
                               n_DMP = sum(s$FDR < FDR_SIG),
                               signo = round(100 * mean(sign(me$logFC_M) == sign(me$logFC_M_ref)), 1),
                               # slope through the origin: refit effect = slope x primary effect
                               pendiente = round(unname(coef(lm(me$logFC_M ~ 0 + me$logFC_M_ref))[1]), 3),
                               dublin = FALSE))
}
# block order in the facets
rob$bloque <- factor(rob$bloque, levels = c("Leave one cohort out", "Technical risk", "Probe panel removed"))
rob <- rob[order(rob$bloque, rob$n_DMP), ]                   # rows sorted by DMP count within each block
rob$etiqueta <- factor(rob$etiqueta, levels = unique(rob$etiqueta))   # keeps that order on the y axis
print(rob, row.names = FALSE)   # table behind Figure 3

# Lollipop panel for one metric: segment from a baseline x0, point, value written to the right,
# one facet row per block (heights proportional to the number of rows).
lollipop <- function(var, x0, fmt, ref, tag, xlab, xlim, brks = waiver()) {
  ggplot(rob, aes(x = .data[[var]], y = etiqueta)) +
    geom_vline(xintercept = ref, linetype = "dashed", linewidth = 0.3, colour = "grey45") +   # reference
    # segment from the baseline
    geom_segment(aes(x = x0, xend = .data[[var]], yend = etiqueta), colour = "grey70", linewidth = 0.4) +
    # open circle at the value
    geom_point(shape = 21, fill = "white", size = 2.4, colour = "grey20", stroke = 0.4) +
    # DMP counts with a thousands separator, as in the text; other metrics with fmt
    geom_text(aes(label = if (var == "n_DMP") formatC(round(.data[[var]]), format = "d", big.mark = ",")
                         else sprintf(fmt, .data[[var]])), nudge_x = 0.05 * diff(xlim), hjust = 0, size = 2.3) +
    facet_grid(bloque ~ ., scales = "free_y", space = "free_y") +   # one facet row per block
    scale_x_continuous(limits = xlim, breaks = brks, expand = expansion(mult = c(0, 0))) +
    labs(tag = tag, x = xlab, y = NULL)
}
# A: DMP count, reference = number of primary DMPs; room on the right for the labels.
pA2 <- lollipop("n_DMP", 0, "%.0f", nrow(sig_prim), "A", "DMPs at\nFDR < 0.05", c(0, 1.3 * max(rob$n_DMP)),
                c(0, 1000, 2000)) +   # x breaks
  theme(strip.text.y = element_blank())
# B: sign retention, baseline and reference at 50% (chance); y labels only on panel A.
pB2 <- lollipop("signo", 50, "%.0f", 50, "B", "Primary DMPs\nkeeping the sign (%)", c(50, 118), c(50, 75, 100)) +
  theme(strip.text.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
# C: effect slope, reference at 1; block names shown only on the rightmost panel.
pC2 <- lollipop("pendiente", 0, "%.2f", 1, "C", "Slope of effects\nvs primary", c(0, 1.55), c(0, 0.5, 1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 7))

fig2 <- (pA2 | pB2 | pC2) + plot_layout(widths = c(1, 0.8, 0.8))   # object fig2 is saved as Figure 3
guardar_fig(fig2, "Figura_3_robustez", 180, 100)   # two-column width
fig2

# ---- 6. Figure 2: reproducibility across cohorts --------------------------------------------------
# A, directional concordance in the primary model and the main analysis (010c); B, the same distribution
# in the main analysis by FDR rank of the DMP; C, within-cohort effect of the 8 concordant loci with 95% CI,
# circle = cohort (same mark for all six), diamond = pooled estimate of the model.
# Per-cohort effects refit the 010c model (no cohort term; lm drops a cell fraction that is collinear
# within a cohort) and are checked against the 010c output before plotting. The forest shows the adjusted
# coefficient on the M scale, the one used to count concordance (the crude delta beta can change sign
# after adjustment; see 6b).

# --- A: concordance in both analyses ---------------------------------------------------------------
# Concordance class of a sign count; 0, 1 and 2 cohorts are pooled into one class.
cat_conc <- function(k) factor(ifelse(k <= 2, "≤2/6", paste0(k, "/6")), levels = names(COL_CONC))
cA <- data.frame()   # long table: arm x class x percent
for (b in c("primario", "noSmk")) {
  r <- conc_res[conc_res$brazo == b, ]   # summary rows of this arm
  for (k in r$cohortes_mismo_signo) {                    # one row per number of concordant cohorts
    cA <- rbind(cA, data.frame(brazo = b, cat = as.character(cat_conc(k)),
                               pct = r$pct_ajustada[r$cohortes_mismo_signo == k]))
  }
}
cA <- aggregate(pct ~ brazo + cat, data = cA, FUN = sum)      # sums 0/6, 1/6 and 2/6 into one class
cA$cat   <- factor(cA$cat, levels = rev(names(COL_CONC)))     # 6/6 stacked at the bottom
cA$brazo <- factor(cA$brazo, levels = c("primario", "noSmk"),
                   labels = c(sprintf("Pre-registered\n(%s DMPs)", format(nrow(sig_prim), big.mark = ",")),
                              sprintf("Main analysis\n(%s DMPs)", format(nrow(sig_main), big.mark = ","))))

# --- B: concordance by FDR rank (main analysis) ----------------------------------------------------
cortes <- c(0, 100, 250, 500, 1000, nrow(sig_main))           # rank bins; the last one runs to the end
sig_main$puesto <- seq_len(nrow(sig_main))                    # rank by FDR (table sorted in section 2)
sig_main$tramo  <- cut(sig_main$puesto, cortes, labels = sprintf("%d–%d", cortes[-6] + 1, cortes[-1]))
cB <- data.frame()   # long table: rank bin x class x percent
for (t in levels(sig_main$tramo)) {
  s <- sig_main[sig_main$tramo == t, ]   # DMPs in this bin
  # percent of DMPs in the bin falling in each concordance class
  for (cc in names(COL_CONC)) cB <- rbind(cB, data.frame(tramo = t, cat = cc, pct = 100 * mean(cat_conc(s$conc_aj) == cc)))
}
cB$tramo <- factor(cB$tramo, levels = levels(sig_main$tramo))
cB$cat   <- factor(cB$cat, levels = rev(names(COL_CONC)))
print(reshape(cB, idvar = "tramo", timevar = "cat", direction = "wide"), row.names = FALSE)   # wide table
# Share of the 100 lowest-FDR DMPs in which GSE147221 goes in the direction of the pooled effect.
cat("Dublín en la dirección del efecto, 100 DMP de menor FDR:",
    sprintf("%.0f %%", 100 * mean(sign(conc[match(head(sig_main$cpg, 100), conc$cpg), DUBLIN]) ==
                                   sign(head(sig_main$logFC_M, 100)))), "\n")

# Stacked 100% bars of concordance classes, shared by panels A and B.
barras_conc <- function(d, x, tag, xlab) {
  ggplot(d, aes(x = .data[[x]], y = pct, fill = cat)) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.25) +   # black outline separates the classes
    scale_fill_manual(values = COL_CONC, name = "Cohorts with\nthe same sign", drop = FALSE) +   # all classes
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(tag = tag, x = xlab, y = "DMPs (%)") +
    guides(fill = guide_legend(reverse = TRUE))          # legend reads 6/6 first
}
pA3 <- barras_conc(cA, "brazo", "A", NULL)
pB3 <- barras_conc(cB, "tramo", "B", "Rank of the DMP by FDR (main analysis)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- C: per-cohort effect of the 8 concordant loci -------------------------------------------------
harm  <- readRDS(file.path(DIR_HARM, "betas_harmonized.rds"))   # harmonised betas, all probes
b8    <- harm$betas[LOCI$cpg, , drop = FALSE]; rm(harm); gc()   # keep only the 12 loci, free the matrix
cov_ep <- as.data.frame(read_csv(F_COV, show_col_types = FALSE))
ANC_COVS <- grep("^ancestry_PC", colnames(cov_ep), value = TRUE)  # EPISTRUCTURE ancestry PCs
p <- ph[ph$set == "train", ]                             # training set (80%), as in the models
rownames(p) <- p$gsm_id   # rows named by GSM accession
i <- match(p$gsm_id, cov_ep$gsm_id)   # row of each sample in the covariate table
p[, ANC_COVS] <- cov_ep[i, ANC_COVS]; p$smoking_SSc <- cov_ep$smoking_SSc[i]   # add PCs and smoking score
# phenotype table that holds the chip position (Sentrix row and column)
pos_ph <- as.data.frame(read_csv("../data_pheno/phenodata_combinado_imputar_edad.csv", show_col_types = FALSE))
p$position <- pos_ph$position[match(p$gsm_id, pos_ph$gsm_id)]       # chip position (main analysis, as 010c)
p$fila <- substr(p$position, 1, 3); p$columna <- substr(p$position, 4, 6)   # "R01C01" -> row "R01", col "C01"
p$group <- factor(p$group, levels = c("CNT", "SCZ")); p$sexo <- factor(p$sexo)   # controls as reference
b8 <- b8[, p$gsm_id, drop = FALSE]                       # same sample order as the phenotype
# M-values with betas bounded away from 0 and 1, as in the differential methylation models
eps <- 1e-6; M8 <- log2(pmin(pmax(b8, eps), 1 - eps) / (1 - pmin(pmax(b8, eps), 1 - eps)))
# Per-cohort adjusted effects written by 010c with the chip-position model, used as a check below.
ref_coh <- as.data.frame(fread(file.path(DIR_CONC, "efecto_ajustado_por_cohorte_pos.csv")))

cohortes <- sort(unique(p$gse))   # the six cohorts
fo <- data.frame()                                       # forest data: one row per locus x cohort
for (j in seq_len(nrow(LOCI))) {
  cg <- LOCI$cpg[j]
  for (g in c(cohortes, "All cohorts")) {
    s <- if (g == "All cohorts") rep(TRUE, nrow(p)) else p$gse == g   # samples of the cohort, or all
    pp <- p[s, ]; pp$y <- M8[cg, s]   # M-values of this locus
    if (g == "All cohorts") {                              # pooled estimate = the main analysis itself
      k <- match(cg, main$cpg)
      m_est <- main$logFC_M[k]; m_se <- main$SE[k]   # adjusted effect and SE from the main model
    } else {
      # within-cohort model: same covariates as the main analysis, without the cohort term
      form <- paste("y ~ group + sexo + age_epigenetic_imputed +", paste(CELL_COVS, collapse = " + "),
                    "+ smoking_SSc +", paste(ANC_COVS, collapse = " + "))
      if (length(unique(pp$sexo)) < 2) form <- sub("sexo \\+ ", "", form)    # single-sex cohort
      if (length(unique(pp$fila)) > 1)    form <- paste(form, "+ factor(fila)")      # chip row, as 010c
      if (length(unique(pp$columna)) > 1) form <- paste(form, "+ factor(columna)")   # chip column
      ct <- summary(lm(as.formula(form), data = pp))$coefficients   # per-cohort OLS fit
      # group coefficient and its SE
      m_est <- ct["groupSCZ", "Estimate"]; m_se <- ct["groupSCZ", "Std. Error"]
      stopifnot(abs(m_est - ref_coh[ref_coh$cpg == cg, g]) < 1e-6)   # must match the 010c output
    }
    # crude delta beta and its Welch standard error, kept for section 6b
    x <- b8[cg, s][pp$group == "SCZ"]; z <- b8[cg, s][pp$group == "CNT"]   # betas of cases and controls
    db <- mean(x) - mean(z); db_se <- sqrt(var(x) / length(x) + var(z) / length(z))
    # 95% CI on the M scale and on the delta beta scale (%)
    fo <- rbind(fo, data.frame(cpg = cg, cohorte = g, n = sum(s), m = m_est, m_lo = m_est - 1.96 * m_se,
                               m_hi = m_est + 1.96 * m_se, db = 100 * db, db_lo = 100 * (db - 1.96 * db_se),
                               db_hi = 100 * (db + 1.96 * db_se)))
  }
}
fo$tipo <- factor(ifelse(fo$cohorte == "All cohorts", "Pooled estimate (six cohorts)", "Cohort"),
                  levels = c("Cohort", "Pooled estimate (six cohorts)"))
# y-axis labels with the training n of each cohort; pooled row at the bottom
n_coh <- table(p$gse)
etq_coh <- c(sprintf("%s (n = %d)", names(n_coh), as.integer(n_coh)), "All cohorts")
names(etq_coh) <- c(names(n_coh), "All cohorts")
fo$cohorte_lab <- factor(etq_coh[fo$cohorte], levels = rev(etq_coh))
cat("efectos por cohorte verificados contra el 010c para", nrow(LOCI), "loci\n")

# Forest panel for a set of loci: one facet per locus (header = gene and concordance), free x scales.
forest <- function(cpgs, tag) {
  d <- fo[fo$cpg %in% cpgs, ]   # rows of the requested loci
  L <- LOCI[match(cpgs, LOCI$cpg), ]   # locus metadata in the requested order
  d$locus <- factor(d$cpg, levels = L$cpg, labels = sprintf("%s\n%d/6 concordant", L$nombre, L$conc_aj))
  ggplot(d, aes(x = m, y = cohorte_lab)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +            # no difference
    geom_linerange(aes(xmin = m_lo, xmax = m_hi), linewidth = 0.4, colour = "grey45") +   # 95% CI
    # cohort point or pooled diamond
    geom_point(aes(fill = tipo, shape = tipo), size = 2, colour = "grey20", stroke = 0.4) +
    facet_wrap(~ locus, nrow = 1, scales = "free_x") +   # one panel per locus
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +   # about three ticks per panel
    # white circle = cohort, grey diamond = pooled estimate
    scale_fill_manual(values = setNames(c("white", "grey35"), levels(fo$tipo)), name = NULL) +
    scale_shape_manual(values = setNames(c(21, 23), levels(fo$tipo)), name = NULL) +
    labs(tag = tag, x = "Adjusted effect, cases − controls (M-value)", y = NULL) +
    theme(panel.spacing.x = unit(4, "mm"))
}

# Two rows of four concordant loci; legend and x title only on the second row.
fC1 <- forest(top_conc$cpg[1:4], "C") + theme(legend.position = "none", axis.title.x = element_blank())
fC2 <- forest(top_conc$cpg[5:8], "") + theme(legend.position = "bottom")
# A and B share one concordance legend; the two forest rows go below.
fig3 <- ((pA3 | pB3) + plot_layout(guides = "collect")) / fC1 / fC2 + plot_layout(heights = c(0.8, 1, 1))
guardar_fig(fig3, "Figura_2_reproducibilidad", 180, 200)   # object fig3 is saved as Figure 2

# Figure S6: the 4 lowest-FDR loci, same layout (object figS5)
figS5 <- forest(top_fdr$cpg, NULL) + theme(legend.position = "bottom")
guardar_fig(figS5, "Figura_S6_forest_menor_FDR", 180, 95)
fig3

# ---- 6b. Crude versus adjusted effects ------------------------------------------------------------
# For each locus and cohort, three numbers: crude delta beta (difference of means), adjusted delta beta
# (the 010c model fitted on betas instead of M-values) and the adjusted M-value effect drawn in the forest.
# Crude vs adjusted beta disagreeing points to the covariate adjustment; adjusted beta vs adjusted M
# disagreeing points to the scale. Console output only, no figure.
comp <- data.frame()   # one row per locus x cohort
for (j in seq_len(nrow(LOCI))) {
  cg <- LOCI$cpg[j]
  for (g in cohortes) {
    s <- p$gse == g; pp <- p[s, ]; pp$y <- b8[cg, s]     # outcome on the beta scale here
    form <- paste("y ~ group + sexo + age_epigenetic_imputed +", paste(CELL_COVS, collapse = " + "),
                  "+ smoking_SSc +", paste(ANC_COVS, collapse = " + "))
    if (length(unique(pp$sexo)) < 2) form <- sub("sexo \\+ ", "", form)
    if (length(unique(pp$fila)) > 1)    form <- paste(form, "+ factor(fila)")
    if (length(unique(pp$columna)) > 1) form <- paste(form, "+ factor(columna)")
    db_aj <- unname(coef(lm(as.formula(form), data = pp))["groupSCZ"])   # adjusted delta beta
    f <- fo[fo$cpg == cg & fo$cohorte == g, ]            # crude delta beta and M effect from the forest data
    # signo = signs of crude beta, adjusted beta and adjusted M, in that order (e.g. "+-+")
    comp <- rbind(comp, data.frame(gen = LOCI$nombre[j], cohorte = g, n = f$n,
                                   dbeta_crudo_pct = round(f$db, 2), dbeta_ajustado_pct = round(100 * db_aj, 2),
                                   M_ajustado = round(f$m, 3),
                                   signo = paste0(ifelse(f$db > 0, "+", "-"), ifelse(db_aj > 0, "+", "-"),
                                                  ifelse(f$m > 0, "+", "-"))))
  }
}
print(comp[comp$gen %in% c("JARID2", "NAV1"), ], row.names = FALSE)   # two loci shown as examples
# Counts of sign changes due to adjustment (1st vs 2nd character) and to scale (2nd vs 3rd).
cat("\nLoci × cohortes en que el signo crudo difiere del ajustado en β:", sum(substr(comp$signo, 1, 1) != substr(comp$signo, 2, 2)),
    "| en que el ajustado en β difiere del ajustado en M:", sum(substr(comp$signo, 2, 2) != substr(comp$signo, 3, 3)),
    "de", nrow(comp), "\n")
cat("\nGSE147221 en los 4 loci de menor FDR (Δβ crudo y ajustado, %):\n")
print(comp[comp$cohorte == DUBLIN & comp$gen %in% LOCI$nombre[LOCI$grupo == "Lowest FDR"], ], row.names = FALSE)
# Covariate balance between groups in the smallest cohort (training samples).
cat("\nCovariables por grupo en GSE116378 (entrenamiento):\n")
q <- p[p$gse == "GSE116378", ]   # GSE116378 training samples
print(data.frame(grupo = c("CNT", "SCZ"), n = as.integer(table(q$group)),
                 mujeres = c(sum(q$sexo[q$group == "CNT"] == "F"), sum(q$sexo[q$group == "SCZ"] == "F")),
                 edad = round(tapply(q$age_epigenetic_imputed, q$group, mean), 1),
                 neutrofilos_pct = round(100 * tapply(q$Neutro, q$group, mean), 1),
                 tabaco_SSc = round(tapply(q$smoking_SSc, q$group, mean), 2)), row.names = FALSE)

# ---- 7. Supplementary figures ---------------------------------------------------------------------
# S4, smoking CpGs in the pre-registered primary model: Manhattan with smoking-panel DMPs highlighted
# (AHRR labelled) and share of smoking CpGs among tested positions and DMPs. S5, genomic context of the
# main-analysis DMPs. S3, ancestry components by cohort and group. S2, estimated leukocyte composition by
# cohort and group. (S6 is drawn with Figure 2.) Internal object names figS1..figS4 differ from the
# numbers in the file names.

# --- S4: smoking CpGs in the primary model ---------------------------------------------------------
mp <- prim[, c("cpg", "chr", "pos", "p", "FDR", "gene")]   # columns for the primary Manhattan
mp$chr_n <- suppressWarnings(as.integer(gsub("chr", "", mp$chr)))   # numeric chromosome
mp <- mp[!is.na(mp$chr_n), ]; mp <- mp[order(mp$chr_n, mp$pos), ]   # drop unplaced, sort by position
# Cumulative genomic position, same construction as Figure 1A.
offset <- 0; centros_p <- numeric(0); mp$pos_acum <- NA_real_
for (k in sort(unique(mp$chr_n))) {
  idx <- mp$chr_n == k
  mp$pos_acum[idx] <- mp$pos[idx] + offset
  centros_p <- c(centros_p, offset + (max(mp$pos[idx]) + min(mp$pos[idx])) / 2)
  offset <- offset + max(mp$pos[idx]) + 2e7
}
names(centros_p) <- sort(unique(mp$chr_n))
mp$par <- factor(mp$chr_n %% 2)   # alternating greys
mp$tabaco_dmp <- mp$cpg %in% cpg_tabaco & mp$FDR < FDR_SIG          # primary DMPs in the smoking panel
ahrr <- mp[mp$cpg == "cg05575921", ]; ahrr$lab <- "AHRR (cg05575921)"   # canonical smoking CpG, labelled

pS1a <- ggplot(mp, aes(x = pos_acum, y = -log10(p))) +
  # background points
  rasterise(geom_point(aes(colour = par), size = 0.35, alpha = 0.7, show.legend = FALSE), dpi = 600) +
  # smoking-panel DMPs in red on top of the grey background
  rasterise(geom_point(data = mp[mp$tabaco_dmp, ], size = 0.8, colour = unname(COL["rojo"])), dpi = 600) +
  # dashed line: largest p-value with FDR < 0.05 in the primary model
  geom_hline(yintercept = -log10(max(prim$p[prim$FDR < FDR_SIG])), linetype = "dashed", linewidth = 0.35, colour = "grey30") +
  geom_point(data = ahrr, size = 2, shape = 23, fill = COL["rojo"], colour = "grey20") +   # AHRR diamond
  geom_text_repel(data = ahrr, aes(label = lab), size = 2.3, fontface = "italic", min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = c("grey70", "grey45")) +   # alternating chromosome greys
  scale_x_continuous(breaks = centros_p, labels = names(centros_p), expand = expansion(mult = 0.01)) +
  labs(tag = "A", x = "Chromosome", y = expression(-log[10](italic(p)))) +
  theme(panel.grid.major.x = element_blank())

en_panel <- prim$cpg %in% cpg_tabaco; en_dmp <- sig_prim$cpg %in% cpg_tabaco   # smoking CpG flags
# Wilson score interval for a proportion, returned in percent.
ic_wilson <- function(exitos, n, z = 1.96) {
  ph_ <- exitos / n; den <- 1 + z^2 / n   # point estimate and correction term
  # centre and half-width of the interval
  ce <- (ph_ + z^2 / (2 * n)) / den; h <- z * sqrt(ph_ * (1 - ph_) / n + z^2 / (4 * n^2)) / den
  c(max(0, ce - h), min(1, ce + h)) * 100
}
ic <- ic_wilson(sum(en_dmp), length(en_dmp))             # CI only for the DMP share (panel share is exact)
# one bar for the panel, one for the DMPs
s1b <- data.frame(grupo = factor(c("Tested\npositions", "DMPs"), levels = c("Tested\npositions", "DMPs")),
                  pct = c(100 * mean(en_panel), 100 * mean(en_dmp)), lo = c(NA, ic[1]), hi = c(NA, ic[2]))
pS1b <- ggplot(s1b, aes(x = grupo, y = pct, fill = grupo)) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +   # bars with a dark outline
  # Wilson 95% CI on the DMP bar
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.35, na.rm = TRUE) +
  # value above the bar, or above the CI when there is one
  geom_text(aes(label = sprintf("%.1f", pct), y = ifelse(is.na(hi), pct, hi)), vjust = -0.6, size = 2.4) +
  # fold change of the DMP share over the panel share
  annotate("text", x = 1.5, y = max(s1b$hi, na.rm = TRUE) + 5, size = 2.3,
           label = sprintf("FC = %.1f", s1b$pct[2] / s1b$pct[1])) +
  # grey panel bar, red DMP bar
  scale_fill_manual(values = c("grey80", unname(COL["rojo"])), guide = "none") +
  scale_y_continuous(limits = c(0, max(s1b$hi, na.rm = TRUE) + 9), expand = expansion(mult = c(0, 0))) +
  labs(tag = "B", x = NULL, y = "Smoking-associated CpGs (%)")
figS1 <- pS1a + pS1b + plot_layout(widths = c(4, 1))     # Manhattan wide, bar panel narrow
guardar_fig(figS1, "Figura_S4_tabaco_primario", 180, 75)
cat(sprintf("  CpG de tabaco: %.1f %% del panel, %d de %d DMP (%.1f %%), FC %.1f\n",
            100 * mean(en_panel), sum(en_dmp), length(en_dmp), 100 * mean(en_dmp), mean(en_dmp) / mean(en_panel)))

# --- S5: genomic context (main analysis) -----------------------------------------------------------
# Relation to CpG islands (Zhou CGIposition) and to genes (promoter = TSS -1500/+500, gene body, intergenic).
clas_isla <- function(x) ifelse(is.na(x) | x == "", "Open sea", ifelse(x == "Island", "CpG island", ifelse(grepl("Shore", x), "Island shore",
                         ifelse(grepl("Shelf", x), "Shelf", "Open sea"))))
clas_gen  <- function(gene, prom) ifelse(prom, "Promoter", ifelse(is.na(gene) | gene == "", "Intergenic", "Gene body"))
main$ctx_isla <- clas_isla(main$island); main$ctx_gen <- clas_gen(main$gene, main$is_promoter)   # background
sig_main$ctx_isla <- clas_isla(sig_main$island); sig_main$ctx_gen <- clas_gen(sig_main$gene, sig_main$is_promoter)

# Per category: percent among tested positions and among DMPs, fold change, chi-square contribution
# (observed vs expected DMP count) and Wilson CI of the DMP percent.
dist_ctx <- data.frame()   # one row per category
for (dim in c("ctx_isla", "ctx_gen")) {   # island context, then gene context
  cats <- if (dim == "ctx_isla") c("CpG island", "Island shore", "Shelf", "Open sea") else c("Promoter", "Gene body", "Intergenic")
  for (k in cats) {
    # background percent and DMP count
    p_fondo <- 100 * mean(main[[dim]] == k); n_dmp <- sum(sig_main[[dim]] == k)
    p_dmp <- 100 * n_dmp / nrow(sig_main); esp <- nrow(sig_main) * p_fondo / 100   # expected DMP count
    ic <- ic_wilson(n_dmp, nrow(sig_main))   # Wilson CI of the DMP percent
    dist_ctx <- rbind(dist_ctx, data.frame(dim = dim, cat = k, Tested = p_fondo, DMP = p_dmp, FC = p_dmp / p_fondo,
                                           chi = (n_dmp - esp)^2 / esp, lo = ic[1], hi = ic[2]))
  }
}
print(dist_ctx, row.names = FALSE)

# Paired bars (tested vs DMP) per category, one facet per category, FC and chi-square written above.
barras_ctx <- function(d, tag, titulo) {
  # long format: one row per category x group; CI only on the DMP bar
  largo <- rbind(data.frame(cat = d$cat, grupo = "Tested positions", pct = d$Tested, lo = NA_real_, hi = NA_real_),
                 data.frame(cat = d$cat, grupo = "DMP", pct = d$DMP, lo = d$lo, hi = d$hi))
  largo$cat <- factor(largo$cat, levels = d$cat)
  largo$grupo <- factor(largo$grupo, levels = c("Tested positions", "DMP"))
  anot <- data.frame(cat = factor(d$cat, levels = d$cat), lab = sprintf("FC = %.2f\nchi2 = %.1f", d$FC, d$chi),
                     y = pmax(d$Tested, d$hi) + 11)      # annotation placed above the taller bar
  ggplot(largo, aes(x = grupo, y = pct, fill = grupo, colour = grupo)) +
    geom_col(width = 0.62, linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.35, colour = "grey20", na.rm = TRUE) +
    geom_text(aes(label = sprintf("%.1f", pct), y = ifelse(is.na(hi), pct, hi)), vjust = -0.6, size = 2.3, colour = "grey20") +
    geom_text(data = anot, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 2.1, colour = "grey25") +
    facet_wrap(~ cat, nrow = 1) +
    scale_fill_manual(values = c("Tested positions" = "grey80", "DMP" = unname(COL["azul"])), name = NULL) +
    scale_colour_manual(values = c("Tested positions" = "grey40", "DMP" = "grey20"), name = NULL) +
    scale_y_continuous(limits = c(0, max(80, max(anot$y) + 8)), expand = expansion(mult = c(0, 0.02))) +
    labs(tag = tag, title = titulo, x = NULL, y = "Percentage of probes") +
    # groups identified by the legend, so x labels are hidden
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
}
# A, CpG island context; B, gene context; one shared legend at the bottom.
figS2 <- barras_ctx(dist_ctx[dist_ctx$dim == "ctx_isla", ], "A", "Relation to CpG islands") /
         barras_ctx(dist_ctx[dist_ctx$dim == "ctx_gen", ], "B", "Location relative to genes") +
         plot_layout(guides = "collect") & theme(legend.position = "bottom")
guardar_fig(figS2, "Figura_S5_contexto", 180, 130)

# --- S3: ancestry and S2: cell composition (training set) ------------------------------------------
tr <- ph[ph$set == "train", ]
pal_coh <- unname(COL[c("azul", "naranja", "verde", "rojo", "morado", "celeste")])   # one colour per cohort
# S3: first two EPISTRUCTURE components; colour = cohort, shape = group.
figS3 <- ggplot(tr, aes(x = ancestry_PC1, y = ancestry_PC2)) +
  geom_point(aes(colour = gse, shape = group), size = 1, alpha = 0.85, stroke = 0.35) +
  scale_colour_manual(values = pal_coh, name = "Cohort") +
  scale_shape_manual(values = c(CNT = 1, SCZ = 3), name = "Group") +   # open circle, cross
  labs(x = "Ancestry principal component 1", y = "Ancestry principal component 2")
guardar_fig(figS3, "Figura_S3_ancestria", 140, 105)

# S2: EpiDISH fractions in long format, one facet per cell type, boxes by cohort and group.
celulas <- c("Neutro", "CD4T", "CD8T", "NK", "Mono", "B")
larga <- data.frame()
for (cc in celulas) larga <- rbind(larga, data.frame(gse = tr$gse, group = tr$group, celula = cc, valor = 100 * tr[[cc]]))
larga$celula <- factor(larga$celula, levels = celulas,
                       labels = c("Neutrophils", "CD4+ T", "CD8+ T", "NK", "Monocytes", "B cells"))
figS4 <- ggplot(larga, aes(x = gse, y = valor, fill = group)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.35, linewidth = 0.3, width = 0.65) +
  facet_wrap(~ celula, scales = "free_y", nrow = 2) +    # free y: fractions differ by an order of magnitude
  scale_fill_manual(values = c(CNT = "grey85", SCZ = unname(COL["azul"])), name = NULL) +
  labs(x = NULL, y = "Estimated cell fraction (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5), legend.position = "bottom")
guardar_fig(figS4, "Figura_S2_celularidad", 180, 120)

# ---- 8. Final check -------------------------------------------------------------------------------
# Figures that the Results text must reproduce. Figure legends are kept with the manuscript, not here.
cat("Archivos en", DIR_FIG, ":\n"); print(list.files(DIR_FIG, pattern = "\\.tiff$"))
cat("\nCifras de control:\n")
cat("  DMP primario ..................", nrow(sig_prim), "\n")
cat("  DMP análisis principal ........", nrow(sig_main), "\n")
cat("  λ principal antes / después ...", sprintf("%.3f / %.3f", lam_pre, lam_post), "\n")
cat("  concordantes 6/6 (principal) ..", sum(sig_main$conc_aj == 6), "\n")
cat("  LOCO sin GSE147221 ............", loco$n_DMP[loco$escenario == paste0("sin_", DUBLIN)], "DMP,",
    loco$signo_igual_pct[loco$escenario == paste0("sin_", DUBLIN)], "% con el mismo signo\n")
