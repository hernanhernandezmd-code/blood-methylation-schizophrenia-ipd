# 06b — DMRs with DMRcate using the statistics of script 06
#
# Calls DMRs with DMRcate's Gaussian kernel, which smooths the per-CpG statistic along the genome and
# needs no per-CpG entry threshold, so it can find weak but concordant regional signal that dmrff misses
# in blood. DMRcate is not allowed to refit its own simple limma model: the CpGannotated object is built
# by hand from the bacon-corrected statistics of 06, so the regions inherit cohort, cell composition,
# imputed age, ancestry, Rubin's rules and bacon.
# Paper section: Differential Methylation
# Input:  data_dmp_P2_anc[...]/dmp_resultados.csv (from 06, same flags);
#         refs/zhou2017_gene_annotation_450k.tsv.gz
# Output: same folder: dmr_dmrcate.csv (significant regions) and dmr_dmrcate_cpgs.csv (region-to-CpG map)
# Flags must match the run of 06. Folders by arm:
#   data_dmp_P2_anc/            pre-registered primary: EXCLUIR_CPG_TABACO = FALSE, USAR_POSICION = FALSE
#   data_dmp_P2_anc_noSmk_pos/  main analysis (default): EXCLUIR_CPG_TABACO = TRUE, USAR_POSICION = TRUE
#   data_dmp_P2_anc_noAnc/, data_dmp_P2_anc_noAge/  sensitivities: EXCLUIR_CPG_ANCESTRIA or
#                               EXCLUIR_CPG_EDAD = TRUE, with EXCLUIR_CPG_TABACO and USAR_POSICION = FALSE
# EXCLUIR_COHORTE: "" = all six cohorts; a GSE id = that cohort left out (suffix _5coh)

# ---- 1. Packages, paths and statistics from 06 ---------------------------------------------------

suppressMessages({ library(DMRcate); library(data.table); library(GenomicRanges); library(readr) })

# The folder suffix is rebuilt from the same flags used in 06, so each arm reads and writes its own folder.
VARIANTE_MASKING   <- "P2"       # probe-masking variant (02-05)
USAR_ANCESTRIA     <- TRUE       # EPISTRUCTURE ancestry PCs in the 06 model
EXCLUIR_CPG_TABACO <- TRUE   # TRUE = main analysis, no smoking CpGs (_noSmk); FALSE = pre-registered
EXCLUIR_CPG_ANCESTRIA <- FALSE   # TRUE = sensitivity without ancestry CpGs (_noAnc)
EXCLUIR_CPG_EDAD      <- FALSE   # TRUE = sensitivity without clock CpGs (_noAge)
EXCLUIR_COHORTE    <- ""   # empty = all six cohorts
SUF_COH <- if (nchar(EXCLUIR_COHORTE)) "_5coh" else ""   # suffix when a cohort is left out
SUF_DMP <- paste0("_", VARIANTE_MASKING)   # suffix: masking variant
if (USAR_ANCESTRIA)      SUF_DMP <- paste0(SUF_DMP, "_anc")   # suffix: ancestry PCs
SUF_DMP <- paste0(SUF_DMP, SUF_COH)                # suffix: cohort left out
if (EXCLUIR_CPG_TABACO)  SUF_DMP <- paste0(SUF_DMP, "_noSmk")   # suffix: smoking CpGs removed
if (EXCLUIR_CPG_ANCESTRIA) SUF_DMP <- paste0(SUF_DMP, "_noAnc")   # suffix: ancestry CpGs removed
if (EXCLUIR_CPG_EDAD)      SUF_DMP <- paste0(SUF_DMP, "_noAge")   # suffix: clock CpGs removed
USAR_POSICION <- TRUE      # TRUE = chip row and column in the 06 model (_pos, main analysis)
if (USAR_POSICION)        SUF_DMP <- paste0(SUF_DMP, "_pos")

DIR_DMP <- paste0("../data_dmp", SUF_DMP)                # per-CpG statistics of the active arm
cat("Variante activa:", SUF_DMP, "| leyendo/escribiendo en:", DIR_DMP, "\n")
if (!file.exists(file.path(DIR_DMP, "dmp_resultados.csv")))   # stops if 06 was not run with these flags
  stop("No hay dmp_resultados.csv en ", DIR_DMP, " — corre el 06 con estos mismos flags primero.")

dmp <- as.data.frame(fread(file.path(DIR_DMP, "dmp_resultados.csv")))
# 06 columns used: cpg, chr, pos, t (post-bacon z), p, FDR, delta_beta, es_dmp
cat("CpGs (panel P2+anc):", nrow(dmp), "| cols:", paste(head(names(dmp), 9), collapse = ", "), "\n")
stopifnot(all(c("chr", "pos", "t", "p", "FDR", "delta_beta") %in% names(dmp)))

# ---- 2. CpGannotated object from the 06 statistics -----------------------------------------------
# dmrcate() reads the mcols stat, diff, ind.fdr and is.sig; ranges must be sorted by position.

# DMP seeds: column es_dmp (main analysis: FDR < 0.05 and |adjusted Δβ| >= 0.5%); FDR < 0.05 if absent.
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
gr <- GRanges(seqnames = dmp$chr,
              ranges   = IRanges(start = dmp$pos, width = 1),
              stat     = dmp$t,                          # post-bacon z from 06
              diff     = dmp$delta_beta,                 # Δβ (SCZ − CNT), beta scale, read as meandiff
              ind.fdr  = dmp$FDR,                         # per-CpG FDR (BH on post-bacon p)
              is.sig   = sel_dmp(dmp))                    # DMP flag of the active arm
names(gr) <- dmp$cpg
gr <- sort(sortSeqlevels(gr))                            # DMRcate requires genomic order
cat("GRanges construido:", length(gr), "CpGs |", sum(gr$is.sig), "DMP (semillas)\n")

annotated <- new("CpGannotated", ranges = gr)            # object class expected by dmrcate()

# ---- 3. DMR calling ------------------------------------------------------------------------------
# A wide kernel merges concordant neighbouring CpGs into whole regions rather than loose pairs.

LAMBDA    <- 1000       # kernel bandwidth (bp)
C_DEFAULT  <- 2.8      # kernel scaling factor (sigma = lambda / C)
BETACUTOFF <- 0.00108   # minimum mean |Δβ| per region, set at the technical noise floor; averaging
                       # dilutes regional Δβ, and the regions are used for feature preselection
MIN_CPGS   <- 2         # minimum CpGs per region
FDR_REG    <- 0.05      # regional significance (HMFDR, harmonic-mean FDR)

dmrc <- dmrcate(annotated, lambda = LAMBDA, C = C_DEFAULT,
                betacutoff = BETACUTOFF, min.cpgs = MIN_CPGS)
# Region statistics are read from the dmrc slots; extractRanges() would require DMRcatedata for gene
# annotation, which is done with Zhou 2017 in step 4 instead.
regs <- data.frame(
  coord = dmrc@coord, no.cpgs = dmrc@no.cpgs, min_smoothed_fdr = dmrc@min_smoothed_fdr,
  Stouffer = dmrc@Stouffer, HMFDR = dmrc@HMFDR, Fisher = dmrc@Fisher,
  meandiff = dmrc@meandiff, maxdiff = dmrc@maxdiff, stringsAsFactors = FALSE)
cc <- do.call(rbind, strsplit(gsub("-", ":", regs$coord), ":", fixed = TRUE))   # split chrN:start-end
regs$seqnames <- cc[, 1]; regs$start <- as.integer(cc[, 2]); regs$end <- as.integer(cc[, 3])
regs$width <- regs$end - regs$start + 1   # region width (bp)
cat("Regiones DMRcate (λ=", LAMBDA, ", C=", C_DEFAULT, ", βcut=", BETACUTOFF, ", min.cpgs=", MIN_CPGS, "):",
    nrow(regs), "\n")

# ---- 4. Significance filter, region-to-CpG map, gene annotation, output --------------------------

df <- regs[regs$HMFDR < FDR_REG, ]                        # significant regions only
cat("Regiones DMRcate:", nrow(regs), "| significativas (HMFDR<", FDR_REG, "):", nrow(df), "\n")

# Panel CpGs falling inside each significant region (used by script 09).
grReg <- GRanges(df$seqnames, IRanges(df$start, df$end)); grReg$dmr_id <- seq_len(nrow(df))
ov    <- findOverlaps(gr, grReg)                          # gr = all panel CpGs (step 2)
mapa  <- data.frame(cpg = names(gr)[queryHits(ov)],
                    dmr_id = grReg$dmr_id[subjectHits(ov)],
                    delta_beta = gr$diff[queryHits(ov)],
                    p = dmp$p[match(names(gr)[queryHits(ov)], dmp$cpg)],
                    FDR = dmp$FDR[match(names(gr)[queryHits(ov)], dmp$cpg)],
                    stringsAsFactors = FALSE)
df$dmr_id <- grReg$dmr_id   # region identifier shared with the map
df$cpgs_en_region <- vapply(df$dmr_id, function(id) sum(mapa$dmr_id == id), integer(1))   # CpGs per region
cat("CpGs del panel mapeados a las regiones:", nrow(mapa),
    "| media CpGs/región:", round(mean(df$cpgs_en_region), 1), "\n")

# Unique Zhou 2017 (GENCODE) genes of the CpGs in each region, as in 06.
zhou_file <- file.path("..", "refs", "zhou2017_gene_annotation_450k.tsv.gz")
if (file.exists(zhou_file)) {
  z <- data.table::fread(zhou_file, sep = "\t", header = TRUE, na.strings = "NA",
                         colClasses = "character", data.table = FALSE)
  genes_reg <- character(nrow(df))
  for (i in seq_len(nrow(df))) {
    cps <- mapa$cpg[mapa$dmr_id == df$dmr_id[i]]           # CpGs of region i
    gg  <- unique(unlist(strsplit(z$genesUniq[match(cps, z$probeID)], ";", fixed = TRUE)))   # their genes
    gg  <- gg[!is.na(gg) & nzchar(gg) & gg != "NA"]        # drop intergenic CpGs
    genes_reg[i] <- paste(gg, collapse = ";")
  }
  df$genes <- genes_reg
} else { df$genes <- NA_character_; warning("Zhou no encontrado en ", zhou_file) }

cat("\n===== BLOQUE B (DMRcate) — resumen =====\n")
cat("Distribución nº CpGs (no.cpgs):\n"); print(table(df$no.cpgs))
cat("n>=3:", sum(df$no.cpgs >= 3), "| n>=5:", sum(df$no.cpgs >= 5), "\n")
cat("|meandiff| (Δβ) — rango:", paste(round(range(abs(df$meandiff)), 4), collapse = " .. "), "\n")
cat("\nTop 15 por HMFDR:\n")
top <- df[order(df$HMFDR), c("seqnames","start","end","no.cpgs","meandiff","maxdiff","HMFDR","genes")]
print(utils::head(top, 15), row.names = FALSE)

readr::write_csv(df,   file.path(DIR_DMP, "dmr_dmrcate.csv"))          # significant regions
readr::write_csv(mapa, file.path(DIR_DMP, "dmr_dmrcate_cpgs.csv"))     # region-to-CpG map
cat("\nGuardado: dmr_dmrcate.csv (", nrow(df), "regiones) + dmr_dmrcate_cpgs.csv (", nrow(mapa), "CpGs)\n")
