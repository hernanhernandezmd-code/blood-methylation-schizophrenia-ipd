# Blood DNA methylation in schizophrenia: one-stage IPD mega-analysis of six cohorts

Analysis code for a pre-registered, one-stage individual participant data (IPD) mega-analysis of
blood DNA methylation in schizophrenia. Six public Illumina HumanMethylation450 cohorts (2,945
participants) are reprocessed uniformly from the methylated and unmethylated signal intensities
deposited in GEO and analyzed in a single model with cohort as a covariate.

- Pre-registration: OSF Registries, [2u56r](https://osf.io/2u56r) (4 June 2026), doi:[10.17605/OSF.IO/2U56R](https://doi.org/10.17605/OSF.IO/2U56R).
- Processed data: GEO [GSE349814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349814)
  (harmonized beta values, one file per source series, and the phenotype table).
- Contact: Hernán Guillermo Hernández Hincapié, Universidad Santo Tomás, Bucaramanga, Colombia.

| Series | Cohort | n deposited |
|---|---|---|
| GSE80417 | UCL | 675 |
| GSE84727 | Aberdeen | 847 |
| GSE116378 | DUTCHSCZ | 64 |
| GSE116379 | CHINFAM | 153 |
| GSE147221 | Dublin | 720 |
| GSE152027 | IoPPN | 800 |

In GSE80417 and GSE84727, `disease status` 1 is control and 2 is case. `01_fenodata.R` stops if the
case and control counts do not match those declared by each GEO series.

## Repository layout

```
fase_1/
  pipeline/   R scripts, numbered in execution order
  refs/       published probe lists and annotation used by the scripts, each with a .SOURCE.txt file
```

All scripts are run from `fase_1/pipeline/`. They read and write sibling folders of `fase_1/`
(`data_input/`, `data_pheno/`, `data_curated/`, `data_normalized_P2/`, `data_harmonized_P2/`,
`data_combat_P2/`, `data_ml_P2/`, `data_dmp_P2_anc*/`, `data_enrich_P2_anc*/`, `figuras_SchizBull/`),
which are created by the pipeline and are not part of the repository. Object and file names, and
console messages, are in Spanish; comments are in English.

## Scripts and paper sections

| Script | Content | Paper section |
|---|---|---|
| `01_fenodata.R` | Phenotype tables from the GEO series matrices, checked against the declared counts | Design, Registration and Data Sources; Table 1 |
| `02_descarga_qc.R` | Download of signal intensities, detection filters, probe masking, sex check | Quality Control and Probe Filtering; Supplementary Figure S1 |
| `03_normalizacion.R` | Joint quantile normalization (dasen) and BMIQ | Normalization and Harmonization |
| `04_armonizacion_covariables.R` | Cross-cohort harmonization (regRCPqn); cell fractions (EpiDISH), smoking score (EpiSmokEr), epigenetic age (methylclock) | Normalization and Harmonization; Covariates |
| `04b_epistructure_ancestria.R` | Ancestry principal components (EPISTRUCTURE) | Covariates; Supplementary Figure S3 |
| `04c_contraste_controles_armonizacion.R` | Between-cohort contrast among controls after harmonization | Normalization and Harmonization |
| `05_combat_particion_edad.R` | Age calibration and multiple imputation (mice), 80/20 partition, ComBat for quality control | Covariates; Data Partition |
| `06_dmp_dmr.R` | DMPs: limma with Rubin's rules over 20 imputations and bacon (adjusted effect and SE on the M scale, crude and adjusted Δβ); DMRs with dmrff | Differential Methylation; Table 2; Figure 1 |
| `06b_dmrcate.R` | DMRs with DMRcate from the statistics of `06` | Differential Methylation |
| `06c_prebacon_qq.R` | Test statistics and inflation before bacon | Differential Methylation; Figure 1 |
| `07_escrutinio_masking.R` | Technical flags of significant DMPs | Robustness and Sensitivity Analyses |
| `07b_vif_colinealidad.R` | Collinearity of model covariates (VIF) | Covariates |
| `07c_tabla1_comparacion_grupos.R` | Case-control comparison of sex, age, smoking score and neutrophils | Table 1 |
| `08c_enriquecimiento_gprofiler.R` | Over-representation with the CpG-gene pair as the unit (g:Profiler gene sets) | Functional Enrichment |
| `010_loco.R` | Leave-one-cohort-out analysis | Robustness and Sensitivity Analyses; Figure 3 |
| `010b_sensibilidad_riesgo_tecnico.R` | Refit without cohorts at high technical risk of bias | Robustness and Sensitivity Analyses |
| `010c_concordancia_direccional.R` | Direction of each DMP within each cohort | Robustness and Sensitivity Analyses; Figure 2 |
| `010d_clozapina.R` | Comparison with CpGs associated with clozapine or treatment-resistant schizophrenia (Hannon et al., 2021) | Robustness and Sensitivity Analyses |
| `010e_solape_hannon2021.R` | Overlap and agreement of effects with the schizophrenia DMPs of Hannon et al. (2021) | Robustness and Sensitivity Analyses; Discussion |
| `010f_sensibilidades_registradas_pendientes.R` | Registered sensitivity analyses: without cell fractions, without ancestry components, after ComBat | Robustness and Sensitivity Analyses; Supplementary Table S1 |
| `010g_sensibilidad_posicion_chip.R` | Chip row and column added to the model | Robustness and Sensitivity Analyses |
| `010h_sensibilidad_emparejamiento.R` | Refit in cases and controls matched 1:1 by sex and age within each cohort | Robustness and Sensitivity Analyses |
| `010i_cohorte_no_es_ruido.R` | Additivity of between-cohort differences in cases and controls; within-cohort concordance against a permutation null | Robustness and Sensitivity Analyses |
| `010j_direccion_cpg_tabaco.R` | Direction of DMPs at smoking-associated CpGs against the effect of smoking | Robustness and Sensitivity Analyses; Discussion |
| `010k_sensibilidades_registradas_restantes.R` | Registered sensitivity analyses: smoking class, epigenetic age acceleration, without regRCPqn | Robustness and Sensitivity Analyses; Supplementary Table S1 |
| `010l_dasen_posicion_chip.R` | dasen with and without its chip-position background refinement | Supplementary Table S1 |
| `11_export_matriz_geo.R` | Export of the processed matrices deposited in GEO | Data and Code Availability |
| `12_figuras_manuscrito.R` | Figures 1-3 and Supplementary Figures S2-S6 | Figures |

The supplementary tables are assembled from the CSV files written by `06`, `06b`, `07`, `08c` and the
`010` series. Each `010` script writes its own folder (`data_dmp_P2_anc_LOCO/`, `_riesgo_tecnico/`,
`_concordancia/`, `_clozapina/`, `_hannon2021/`, `_sens_registradas/`, `_sens_posicion/`, `_sens_emparejado/`,
`_ruido/`, `_tabaco_direccion/`, `_sens_restantes/`, `_dasen_posicion/`).

## Analysis arms

`06_dmp_dmr.R` is run once per arm; later scripts read the folder of the arm selected by their flags.

| Arm | Flags in `06_dmp_dmr.R` | Output folder |
|---|---|---|
| Main analysis: smoking-associated CpGs removed, chip row and column as covariates; DMP = FDR < 0.05 and adjusted \|Δβ\| ≥ 0.5% (column `es_dmp`) | `EXCLUIR_CPG_TABACO <- TRUE`, `USAR_POSICION <- TRUE` (defaults) | `data_dmp_P2_anc_noSmk_pos/` |
| Pre-registered primary model (all probes, FDR < 0.05) | all `EXCLUIR_CPG_*` and `USAR_POSICION` = `FALSE` | `data_dmp_P2_anc/` |
| Main analysis without chip position | `EXCLUIR_CPG_TABACO <- TRUE`, `USAR_POSICION <- FALSE` | `data_dmp_P2_anc_noSmk/` |
| Sensitivity: ancestry CpGs removed | `EXCLUIR_CPG_ANCESTRIA <- TRUE`, `USAR_POSICION <- FALSE` | `data_dmp_P2_anc_noAnc/` |
| Sensitivity: clock CpGs removed | `EXCLUIR_CPG_EDAD <- TRUE`, `USAR_POSICION <- FALSE` | `data_dmp_P2_anc_noAge/` |

`06b`, `07`, `07b` and `08c` default to the main analysis. Chip position and the effect-size threshold were
added after registration and are reported as deviations.

Other flags shared by several scripts:

- `VARIANTE_MASKING`: `"P2"` (default; Chen and Benton cross-reactive probes, Zhou mask, SNPs at the
  CpG or single-base extension with MAF >= 0.05). `"P1"` adds the Naeem list and any probe SNP;
  `"P3"` keeps only Chen and Benton.
- `USAR_ANCESTRIA`: `TRUE` (default) includes four EPISTRUCTURE components in the model.
- `EXCLUIR_COHORTE`: `""` (default, six cohorts). A GEO accession drops that cohort from `03`
  onwards; outputs then go to folders with the suffix `_5coh`.

## Order of execution

`01_fenodata.R` → `02` → `03` → `04` → `04b` → `05` → `06` (each arm) → `04c` →
`06b` → `06c` → `07` → `07b` → `07c` → `08c` → `010` → `010b` → `010c` → `010d` → `010e` → `010f` → `010g` →
`010h` → `010i` → `010j` → `010k` → `010l` → `11` → `12`.

`04c` runs after `06` because it compares the between-cohort contrast with the DMPs of the pre-registered
primary model.

`04` overwrites `data_pheno/covariables_epigeneticas_P2.csv`; `04b` adds the ancestry components to
it, so `04b` is run after every run of `04`.

## Requirements

R 4.3.3 and Bioconductor 3.18. Packages:

- Bioconductor: GEOquery, Biobase, minfi, IlluminaHumanMethylation450kanno.ilmn12.hg19, methylumi,
  illuminaio, wateRmelon, ENmix, lumi, EpiDISH, methylclock, methylclockData, sva, limma, bacon,
  DMRcate, missMethyl, GenomicRanges, BiocParallel.
- CRAN: data.table, readr, readxl, dplyr, tibble, stringr, purrr, glue, jsonlite, httr2, yaml, matrixStats,
  mice, car, ggplot2, ggrastr, ggrepel, patchwork, scales, remotes.
- GitHub: regRCPqn (`regRCPqn/regRCPqn`), EpiSmokEr (`sailalithabollepalli/EpiSmokEr`),
  dmrff (`perishky/dmrff`), maxprobes (`markgene/maxprobes`).

Most scripts install missing packages on first use. `preprocessCore`, used by regRCPqn, may need to
be compiled without threading on some Linux systems:
`BiocManager::install("preprocessCore", configure.args = "--disable-threading", force = TRUE)`.

Scripts `02` to `05` hold all cohorts in memory, because joint quantile normalization and
harmonization are defined over the whole sample set; they were run with 128 GB of RAM. BMIQ in
`03` uses four parallel workers to limit peak memory.

The g:Profiler gene sets used by `08c` are not redistributed; download instructions and the release
are in `fase_1/refs/gprofiler_e114_eg62_p19/SOURCE.txt`. KEGG is not used because its gene sets cannot be
redistributed. The supplementary tables of Joehanes et al. (2016) read by `010j` are not redistributed
either; the source is given in `fase_1/refs/joehanes2016_smoking_cpgs_450k.SOURCE.txt`.

## License

MIT. See `LICENSE`.
