# 08c — Functional enrichment with the CpG-gene pair as the unit
#
# Tests over-representation in the g:Profiler gene sets (release e114_eg62_p19) with duplicated genes
# kept in both lists, so genes with more probes weigh more. The universe is the probes that entered
# limma, not the full 450K array. Hypergeometric test computed locally with Benjamini-Hochberg within
# each source, followed by a probe-level support filter: a term is reported only if it rests on >= 2
# distinct CpGs in >= 2 loci more than 1 Mb apart and no single probe contributes more than half of
# its genes (one probe annotated to several genes, as in the UGT1A locus, inflates the gene count).
# Paper section: Functional Enrichment
# Input:  data_dmp_P2_anc_noSmk_pos/dmp_resultados.csv (from 06); refs/Homo_sapiens.gene_info.gz;
#         refs/gprofiler_e114_eg62_p19/hsapiens.<source>.name.gmt
# Output: data_enrich_P2_anc_noSmk_pos_gp/enriquecimiento_todos_los_terminos.csv (all tested terms),
#         enriquecimiento_crudo.csv (significant terms with probe counts), enriquecimiento_filtrado.csv
#         (terms that pass the probe-level filter)
# EXCLUIR_CPG_TABACO: TRUE (default) reads the main analysis without smoking-associated CpGs (_noSmk);
#   FALSE reads the pre-registered model with all probes.
# USAR_POSICION: TRUE (default) reads the model with chip row/column covariates (_pos, main analysis).

# ---- 1. Packages, paths and parameters ---------------------------------------------------------

suppressMessages({ library(data.table); library(readr) })

# Allows running from the repository root as well as from fase_1/pipeline.
if (!file.exists("../data_dmp_P2_anc") && dir.exists("fase_1/pipeline")) setwd("fase_1/pipeline")

EXCLUIR_CPG_TABACO <- TRUE   # TRUE = main analysis without smoking CpGs (_noSmk)
USAR_POSICION <- TRUE        # TRUE = chip row/column in the model (_pos)
SUF      <- paste0("_P2_anc", if (EXCLUIR_CPG_TABACO) "_noSmk" else "", if (USAR_POSICION) "_pos" else "")
DIR_DMP  <- paste0("../data_dmp", SUF)                   # DMP table of the active analysis
DIR_REFS <- "../refs"                                    # cached gene catalogue and gene sets
DIR_OUT  <- paste0("../data_enrich", SUF, "_gp")         # output folder
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

FDR_LISTA  <- 0.05         # CpG enters the list (selection); the list itself is read from es_dmp
PADJ_SIG   <- 0.05         # term significance (inference)
MIN_CPGS  <- 2             # support filter: at least 2 distinct CpGs
MIN_LOCI  <- 2             # at least 2 loci
SEP_LOCI  <- 1e6           # CpGs more than 1 Mb apart count as separate loci
MAX_SONDA <- 0.50          # no single probe may contribute more than half of the term's genes
# KEGG is not used because its gene sets cannot be redistributed.
FUENTES   <- c("GO:BP", "GO:CC", "GO:MF", "REAC", "WP", "HP")
DIR_GMT   <- file.path(DIR_REFS, "gprofiler_e114_eg62_p19")     # g:Profiler gene sets (see SOURCE.txt)
MIN_GENES <- 2             # minimum term size within the panel

# ---- 2. Gene symbol to Entrez GeneID; build both lists -----------------------------------------
# NCBI gene_info resolves official symbols first and synonyms second (recovers Zhou/GENCODE aliases).
info <- fread(file.path(DIR_REFS, "Homo_sapiens.gene_info.gz"), sep = "\t", quote = "",
              select = c("GeneID", "Symbol", "Synonyms"))
mapa <- new.env(hash = TRUE, parent = emptyenv())        # symbol -> GeneID lookup table
for (i in seq_len(nrow(info))) assign(info$Symbol[i], as.character(info$GeneID[i]), envir = mapa)
for (i in seq_len(nrow(info)))                           # synonyms never overwrite an official symbol
  for (s in strsplit(info$Synonyms[i], "|", fixed = TRUE)[[1]])
    if (s != "-" && !exists(s, envir = mapa, inherits = FALSE))
      assign(s, as.character(info$GeneID[i]), envir = mapa)
a_id <- function(g) if (exists(g, envir = mapa, inherits = FALSE)) get(g, envir = mapa) else NA_character_

dmp <- as.data.frame(fread(file.path(DIR_DMP, "dmp_resultados.csv")))
# DMP of the main analysis = FDR < 0.05 and |adjusted delta-beta| >= 0.5% (column es_dmp from 06);
# tables without that column (pre-registered primary model) fall back to FDR < 0.05.
sel_dmp <- function(d) if ("es_dmp" %in% names(d)) as.logical(d$es_dmp) else d$FDR < 0.05
es_dmp <- sel_dmp(dmp)                                   # TRUE for probes in list 1
cat("sondas que entraron a limma:", nrow(dmp),
    "| chrX/chrY:", sum(dmp$chr %in% c("chrX", "chrY")), "(deben ser 0)\n")

# Genes without an Entrez ID are dropped from both lists, so numerator and denominator describe the
# same set of measured probes.
son_gen <- new.env(hash = TRUE, parent = emptyenv())     # gene -> significant CpGs
gen_cpg <- new.env(hash = TRUE, parent = emptyenv())     # CpG -> genes it is annotated to
pos_cpg <- new.env(hash = TRUE, parent = emptyenv())     # CpG -> "chr pos"
uni_id <- character(0); q_id <- character(0)
for (i in seq_len(nrow(dmp))) {
  if (is.na(dmp$gene[i]) || dmp$gene[i] == "") next      # intergenic probes
  sig <- es_dmp[i]                                       # probe is a DMP
  if (sig) assign(dmp$cpg[i], paste(dmp$chr[i], dmp$pos[i]), envir = pos_cpg)
  for (g in strsplit(dmp$gene[i], ";", fixed = TRUE)[[1]]) {
    id <- a_id(g); if (is.na(id)) next                   # skip genes without Entrez ID
    uni_id <- c(uni_id, id)                              # background record
    if (sig) {
      q_id <- c(q_id, id)                                 # duplicates kept
      k <- paste0("g", id)
      assign(k, c(if (exists(k, envir = son_gen, inherits = FALSE)) get(k, envir = son_gen),
                  dmp$cpg[i]), envir = son_gen)          # add this CpG to the gene's list
      # and the gene to the CpG's list
      assign(dmp$cpg[i], c(if (exists(dmp$cpg[i], envir = gen_cpg, inherits = FALSE))
                             get(dmp$cpg[i], envir = gen_cpg), id), envir = gen_cpg)
    }
  }
}
uni_rec <- uni_id                                         # background: one record per CpG-gene pair analysed
q_rec   <- q_id                                           # list 1: one record per DMP-gene pair
uni_id  <- unique(uni_id)                                 # distinct panel genes
cat(sprintf("LISTA GRANDE: %d registros CpG–gen (%d GeneID únicos) · LISTA 1: %d registros (%d GeneID únicos)\n",
            length(uni_rec), length(uni_id), length(q_rec), length(unique(q_rec))))
stopifnot(all(unique(q_rec) %in% uni_id))                 # list 1 lies within the background

# ---- 3. Hypergeometric test on CpG-gene records ------------------------------------------------
# Per source, the universe is the background records whose gene is annotated in that source.
# N = universe records, K = records in the term, n = list records in the universe, k = list records
# in the term; P = P(X >= k). BH is applied over all terms tested in the source (untested terms with
# k = 0 enter as p = 1).
n_bg <- table(uni_rec); n_q <- table(q_rec)               # records per gene
cuenta <- function(ids, tb) { v <- tb[intersect(ids, names(tb))]; if (length(v)) sum(v) else 0 }

res <- data.frame()
for (src in FUENTES) {
  lineas <- strsplit(readLines(file.path(DIR_GMT, paste0("hsapiens.", src, ".name.gmt"))), "\t", fixed = TRUE)
  conj <- list(); nombre <- character(0)                 # gene sets and their names
  for (l in lineas) {
    ids <- unique(na.omit(vapply(l[-(1:2)], a_id, character(1))))   # GMT symbols to GeneID
    ids <- intersect(ids, uni_id)                          # restrict to panel genes
    if (length(ids) >= MIN_GENES) { conj[[l[1]]] <- ids; nombre[l[1]] <- l[2] }
  }
  anotados <- unique(unlist(conj))                         # panel genes annotated in this source
  N <- cuenta(anotados, n_bg); n <- cuenta(anotados, n_q)
  t_src <- data.frame()
  for (tm in names(conj)) {
    K <- cuenta(conj[[tm]], n_bg); k <- cuenta(conj[[tm]], n_q)
    if (k == 0) next                                       # no list record in the term
    t_src <- rbind(t_src, data.frame(fuente = src, termino = tm, descripcion = nombre[[tm]],
                                     k = k, n = n, K = K, N = N,
                                     p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
                                     tam_termino = length(conj[[tm]]), stringsAsFactors = FALSE))
  }
  m_tests <- length(conj)                                  # BH over all terms tested in the source
  t_src$padj <- pmin(1, p.adjust(c(t_src$p, rep(1, m_tests - nrow(t_src))), method = "BH")[seq_len(nrow(t_src))])
  t_src$genes_termino <- vapply(t_src$termino, function(tm) paste(conj[[tm]], collapse = ";"), character(1))
  res <- rbind(res, t_src)
  cat(sprintf("  %-6s términos probados %5d | universo %7d registros (%5d genes) | lista 1 en el universo %5d | padj<%.2f: %d\n",
              src, m_tests, N, length(anotados), n, PADJ_SIG, sum(t_src$padj < PADJ_SIG)))
}
res$enriquecimiento <- round((res$k / res$n) / (res$K / res$N), 2)   # fold enrichment

# ---- 4. Probe-level support filter -------------------------------------------------------------
# For each significant term: distinct supporting CpGs, loci separated by > 1 Mb, and the share of the
# term's list genes contributed by the single most contributing probe.
q_genes <- unique(q_rec)
tab <- data.frame()
sig <- res[res$padj < PADJ_SIG, ]                         # significant terms
for (j in seq_len(nrow(sig))) {
  # List genes in the term, then the significant CpGs annotated to them.
  ids  <- intersect(strsplit(sig$genes_termino[j], ";", fixed = TRUE)[[1]], q_genes)
  cpgs <- unique(unlist(lapply(ids, function(id) {
    kk <- paste0("g", id); if (exists(kk, envir = son_gen, inherits = FALSE)) get(kk, envir = son_gen) })))
  ub <- unlist(lapply(cpgs, function(c) if (exists(c, envir = pos_cpg, inherits = FALSE)) get(c, envir = pos_cpg)))
  nl <- 0                                                 # number of separate loci
  if (length(ub)) {
    cr <- sub(" .*", "", ub); ps <- as.numeric(sub(".* ", "", ub))   # chromosome and position
    for (c1 in unique(cr)) { v <- sort(ps[cr == c1]); nl <- nl + 1
      if (length(v) > 1) nl <- nl + sum(diff(v) > SEP_LOCI) }   # a gap > 1 Mb starts a new locus
  }
  # Probe that contributes most genes to this term.
  dom_c <- NA_character_; dom_n <- 0
  for (c1 in cpgs) {
    g1 <- if (exists(c1, envir = gen_cpg, inherits = FALSE)) get(c1, envir = gen_cpg) else character(0)
    k1 <- length(intersect(g1, ids))                      # term genes annotated to this probe
    if (k1 > dom_n) { dom_n <- k1; dom_c <- c1 }
  }
  tab <- rbind(tab, data.frame(
    fuente = sig$fuente[j], termino = sig$termino[j], descripcion = sig$descripcion[j],
    padj = sig$padj[j], enriquecimiento = sig$enriquecimiento[j], registros_lista1 = sig$k[j],
    n_genes = length(ids), n_cpgs = length(cpgs), n_loci = nl,
    sonda_dominante = dom_c, genes_de_esa_sonda = dom_n,
    pct_sonda_dominante = round(100 * if (length(ids)) dom_n / length(ids) else NA_real_, 1),
    tam_termino = sig$tam_termino[j], universo_registros = sig$N[j],
    genes = paste(ids, collapse = ";"), cpgs = paste(cpgs, collapse = ";"), stringsAsFactors = FALSE))
}
if (nrow(tab)) tab <- tab[order(tab$padj), ]
write_csv(res[order(res$padj), setdiff(names(res), "genes_termino")], file.path(DIR_OUT, "enriquecimiento_todos_los_terminos.csv"))
write_csv(tab, file.path(DIR_OUT, "enriquecimiento_crudo.csv"))   # significant terms with probe counts

# Terms that pass the probe-level filter.
fin <- if (nrow(tab)) tab[tab$n_cpgs >= MIN_CPGS & tab$n_loci >= MIN_LOCI &
                         tab$pct_sonda_dominante <= 100 * MAX_SONDA, ] else tab
write_csv(fin, file.path(DIR_OUT, "enriquecimiento_filtrado.csv"))

cat(sprintf("\nTérminos significativos (padj<%.2f): %d | superan el control de sondas: %d | descartados por artefacto: %d\n",
            PADJ_SIG, nrow(tab), nrow(fin), nrow(tab) - nrow(fin)))
des <- if (nrow(tab)) tab[!(tab$termino %in% fin$termino), ] else tab   # terms removed by the filter
if (nrow(des)) {
  cat("\nDESCARTADOS por sostenerse en una sola sonda o un solo locus (NO van al manuscrito):\n")
  des$motivo <- ifelse(des$n_cpgs < MIN_CPGS | des$n_loci < MIN_LOCI,
                       "una sola sonda/locus", "una sonda domina el término")
  print(head(des[, c("fuente", "descripcion", "padj", "n_cpgs", "n_loci", "pct_sonda_dominante", "motivo")], 15), row.names = FALSE)
}
if (nrow(fin)) {
  cat("\nTÉRMINOS QUE SE REPORTAN:\n")
  print(fin[, c("fuente", "descripcion", "padj", "enriquecimiento", "n_genes", "n_cpgs", "n_loci", "pct_sonda_dominante")], row.names = FALSE)
} else cat("\nNingún término supera el control: el resultado es NULO y así se reporta.\n")
