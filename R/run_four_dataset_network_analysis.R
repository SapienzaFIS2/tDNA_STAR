###############################################################################
############### RUN: GBM, CRC, NSCLC OR OVARIAN DIFFERENTIAL NETWORK ANALYSIS ###
###############################################################################
rm(list = ls(all.names = TRUE))
set.seed(123)

###############################################################################
# 1. LIBRARIES
###############################################################################

library(readxl)
library(GEOquery)
library(dplyr)
library(tibble)
library(DESeq2)
library(igraph)
library(readr)

source("functions_four_dataset.R")

###############################################################################
# 2. DATASET SELECTION
###############################################################################

# Choose one dataset:
#dataset <- "GBM"
#dataset <- "CRC"
dataset <- "NSCLC"
#dataset <- "Ovarian"

###############################################################################
# 3. PATHS AND DATASET-SPECIFIC SETTINGS
###############################################################################

gene_annotation_path <- "Human.GRCh38.p13.annot.tsv"

if (dataset == "GBM") {
  
  info_pat_path <- "GSE156902.xlsx"
  counts_path <- "GSE156902_raw_counts_GRCh38.p13_NCBI.tsv"
  output_dir <- "dataset_sol_GBM"
  
  control_group <- "AC"
  disease_group <- "GBM"
  
} else if (dataset == "CRC") {
  
  info_pat_path <- "2022_xu_et_al.csv"
  counts_path <- "2022_xu_et_al.txt"
  output_dir <- "dataset_sol_CRC"
  
  control_group <- "HD"
  disease_group <- "CRC"
  
} else if (dataset == "NSCLC") {
  
  info_pat_path <- "41598_2023_35818_MOESM2_ESM.xlsx"
  counts_path <- "GSE207586_NSCLC_CountMatrix_2ndVersion.txt"
  output_dir <- "dataset_sol_NSCLC"
  
  control_group <- "Control"
  disease_group <- "NSCLC"
  
} else if (dataset == "Ovarian") {
  
  info_pat_path <- "mmc3.xlsx"
  counts_path <- "GSE183635_raw_counts_GRCh38.p13_NCBI.tsv"
  output_dir <- "dataset_sol_Ovarian"
  
  control_group <- "Asymptomatic controls"
  disease_group <- "Ovarian cancer"
  
} else {
  
  stop("dataset must be one of: 'GBM', 'CRC', 'NSCLC', 'Ovarian'")
}

###############################################################################
# 4. LOAD DATA
###############################################################################

gene_annotation <- read.delim(
  gene_annotation_path,
  row.names = 1
)

if (dataset == "GBM") {
  
  info_pat <- read_excel(
    info_pat_path,
    skip = 2
  )
  
  counts <- read.table(
    counts_path,
    header = TRUE,
    row.names = 1
  )
  
} else if (dataset == "CRC") {
  
  info_pat <- read_delim(
    info_pat_path,
    delim = ";",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  counts <- read_delim(
    counts_path,
    delim = "\t"
  )
  
  counts <- counts %>%
    column_to_rownames(var = "gene_id")
  
  counts <- as.data.frame(counts)
  
} else if (dataset == "NSCLC") {
  
  info_pat <- read_excel(
    info_pat_path,
    skip = 2
  )
  
  counts <- read.table(
    counts_path,
    header = TRUE,
    row.names = 1
  )
  
} else if (dataset == "Ovarian") {
  
  info_pat <- read_excel(
    info_pat_path,
    skip = 2
  )
  
  counts <- read.table(
    counts_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
}

###############################################################################
# 5. PREPARE METADATA
###############################################################################

if (dataset == "GBM") {
  
  GBM <- prepare_gbm_metadata(info_pat)
  
  gse <- getGEO(
    "GSE156902",
    GSEMatrix = TRUE
  )
  
  geo_meta <- build_geo_metadata(gse)
  
  meta <- merge_metadata_with_geo(
    GBM = GBM,
    geo_meta = geo_meta
  )
  
  meta <- apply_manual_metadata_fixes(meta)
  
  aligned <- align_counts_metadata(
    counts = counts,
    meta = meta
  )
  
  counts <- aligned$counts
  meta <- aligned$meta
  
} else if (dataset == "CRC") {
  
  meta <- prepare_crc_metadata(
    info_pat = info_pat,
    counts = counts
  )
  
  valid_samples <- intersect(
    rownames(meta),
    colnames(counts)
  )
  
  counts <- counts[
    ,
    valid_samples,
    drop = FALSE
  ]
  
  meta <- meta[
    valid_samples,
    ,
    drop = FALSE
  ]
  
} else if (dataset == "NSCLC") {
  
  prepared <- prepare_nsclc_metadata(
    info_pat = info_pat,
    counts = counts
  )
  
  counts <- prepared$counts
  meta <- prepared$meta
  
  valid_samples <- intersect(
    rownames(meta),
    colnames(counts)
  )
  
  counts <- counts[
    ,
    valid_samples,
    drop = FALSE
  ]
  
  meta <- meta[
    valid_samples,
    ,
    drop = FALSE
  ]
  
} else if (dataset == "Ovarian") {
  
  gse <- getGEO(
    "GSE183635",
    GSEMatrix = TRUE
  )
  
  geo_meta <- build_ovarian_geo_metadata(gse)
  
  meta <- prepare_ovarian_metadata(
    info_pat = info_pat,
    counts = counts,
    geo_meta = geo_meta
  )
  
  valid_samples <- intersect(
    rownames(meta),
    colnames(counts)
  )
  
  counts <- counts[
    ,
    valid_samples,
    drop = FALSE
  ]
  
  meta <- meta[
    valid_samples,
    ,
    drop = FALSE
  ]
}

cat("\nDataset:\n")
print(dataset)

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\nStage table:\n")
print(table(meta$Stage, useNA = "ifany"))

cat("\nGroup x Stage table:\n")
print(table(meta$Group, meta$Stage, useNA = "ifany"))

###############################################################################
# 6. CONTROL VS DISEASE DIFFERENTIAL EXPRESSION
###############################################################################

if (dataset == "GBM") {
  
  contrast_data <- select_ac_vs_t0(
    meta = meta,
    counts = counts
  )
  
  counts_contrast <- contrast_data$counts
  meta_contrast <- contrast_data$meta
  
} else if (dataset %in% c("CRC", "NSCLC", "Ovarian")) {
  
  counts_contrast <- counts
  meta_contrast <- meta
}

cat("\nContrast table:\n")
print(table(meta_contrast$Group, meta_contrast$Stage, useNA = "ifany"))

de_control_vs_disease <- run_deseq2_contrast(
  counts_mat = counts_contrast,
  metadata = meta_contrast,
  group_col = "Group",
  contrast_vec = c("Group", control_group, disease_group),
  min_count = 20,
  min_prop = 0.8,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.33
)

dds <- de_control_vs_disease$dds
res_contr_vs_disease <- de_control_vs_disease$res
deg_contr_vs_disease <- de_control_vs_disease$deg

cat("\nDEGs control vs disease:\n")
print(dim(deg_contr_vs_disease))

###############################################################################
# 7. VST MATRIX FOR DE GENES
###############################################################################

expr_vst_filtr <- get_vst_deg_matrix(
  dds = dds,
  deg_table = deg_contr_vs_disease,
  blind = FALSE
)

expr_log2_filtr <- expr_vst_filtr

###############################################################################
# 8. SELECT DATASET-SPECIFIC STAGES FOR NETWORK ANALYSIS
###############################################################################

# GBM remains unchanged: T0-T1-T2-T3.
# For solid-tumour staging, controls (Stage 0) are excluded from diffcorr:
#   CRC             -> Stage I, II, III
#   NSCLC / Ovarian -> Stage I, II, III, IV
if (dataset == "GBM") {
  network_stages <- 0:3
  transitions <- list(c(0, 1), c(1, 2), c(2, 3))
  meta_T <- select_complete_t0_t3_patients(meta)
  
} else if (dataset == "CRC") {
  network_stages <- 1:3
  transitions <- list(c(1, 2), c(2, 3))
  meta_T <- select_stage_samples(meta, stages = network_stages)
  
} else if (dataset %in% c("NSCLC", "Ovarian")) {
  network_stages <- 1:4
  transitions <- list(c(1, 2), c(2, 3), c(3, 4))
  meta_T <- select_stage_samples(meta, stages = network_stages)
}

cat("\nStages used for network analysis:\n")
print(network_stages)

cat("\nSamples per stage used for network analysis:\n")
print(table(meta_T$Stage, useNA = "ifany"))

cat("\nNumber of samples used for network analysis:\n")
print(nrow(meta_T))

# Fisher's z comparison requires >3 samples in every compared stage.
stage_n <- table(factor(meta_T$Stage, levels = network_stages))
if (any(stage_n <= 3)) {
  bad_stages <- names(stage_n)[stage_n <= 3]
  stop(
    paste0(
      "At least 4 samples are required in every stage used for diffcorr. ",
      "Insufficient stages: ", paste(bad_stages, collapse = ", ")
    )
  )
}

###############################################################################
# 9. NORMALIZED EXPRESSION MATRIX FOR STAGE NETWORK ANALYSIS
###############################################################################

deg_genes <- rownames(deg_contr_vs_disease)
samples_network <- rownames(meta_T)

counts_network <- counts[
  ,
  samples_network,
  drop = FALSE
]

normalized_counts_filtr <- get_normalized_deg_matrix(
  counts_mat = counts_network,
  metadata = meta_T,
  deg_genes = deg_genes
)

cat("\nFiltered normalized matrix dimensions:\n")
print(dim(normalized_counts_filtr))

###############################################################################
# 10. STAGE-SPECIFIC EXPRESSION MATRICES
###############################################################################

stage_expr_list <- make_stage_expr_list(
  expr_mat = normalized_counts_filtr,
  metadata = meta_T,
  stage_col = "Stage"
)

cat("\nStage expression dimensions:\n")
print(lapply(stage_expr_list, dim))

###############################################################################
# 11. CORRELATION MATRICES
###############################################################################

cor_matrices <- compute_stage_correlations(
  stage_expr_list = stage_expr_list,
  method = "pearson"
)

###############################################################################
# 12. DIFFERENTIAL NETWORKS BETWEEN DATASET-SPECIFIC STAGES
###############################################################################

transition_results <- analyze_consecutive_transitions(
  cor_matrices = cor_matrices,
  metadata = meta_T,
  gene_annotation = gene_annotation,
  transitions = transitions,
  z_cutoff = 3,
  percentile = 0.90
)

network_summary <- summarise_transition_results(
  transition_results
)

cat("\nNetwork summary:\n")
print(network_summary)

###############################################################################
# 13. EXTRACT MAIN GENE SETS
###############################################################################

intersection_genes <- lapply(
  transition_results,
  function(x) x$key_genes$genes_intersection
)

hubs_symbols <- lapply(
  transition_results,
  function(x) x$key_genes$hubs_symbols
)

betweenness_symbols <- lapply(
  transition_results,
  function(x) x$key_genes$betweenness_symbols
)

closeness_symbols <- lapply(
  transition_results,
  function(x) x$key_genes$closeness_symbols
)

cat("\nIntersection genes per transition:\n")
for (nm in names(intersection_genes)) {
  cat(nm, ":", length(intersection_genes[[nm]]), "\n")
}

###############################################################################
# 14. SAVE BASIC OUTPUTS
###############################################################################
'''
dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

save_transition_genes(
  transition_results = transition_results,
  output_dir = output_dir
)

write.csv(
  as.data.frame(res_contr_vs_disease),
  file = file.path(output_dir, "deseq2_results_control_vs_disease.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(deg_contr_vs_disease),
  file = file.path(output_dir, "deseq2_deg_control_vs_disease.csv"),
  row.names = TRUE
)

write.csv(
  network_summary,
  file = file.path(output_dir, "network_summary.csv"),
  row.names = FALSE
)

for (nm in names(intersection_genes)) {
  write.csv(
    data.frame(genes = intersection_genes[[nm]]),
    file = file.path(output_dir, paste0("genes_", nm, ".csv")),
    row.names = FALSE
  )
}
'''
###############################################################################
# 15. HUBS / CENTRAL GENES OVERLAP
###############################################################################

# Genes recurrent in every transition available for the selected dataset.
# Examples:
#   GBM / NSCLC / Ovarian: intersection across three transitions
#   CRC: intersection across two transitions
intersect_all <- function(x) {
  if (length(x) == 0) return(character(0))
  if (length(x) == 1) return(unique(x[[1]]))
  Reduce(intersect, x)
}

recurrent_hubs <- intersect_all(hubs_symbols)
recurrent_betweenness <- intersect_all(betweenness_symbols)
recurrent_closeness <- intersect_all(closeness_symbols)

top_top <- unique(
  c(
    recurrent_hubs,
    recurrent_betweenness,
    recurrent_closeness
  )
)

cat("\nGenes recurrent across all available stage transitions:\n")
print(top_top)

# Pairwise overlap rows, followed by overlap across all transitions.
make_pairwise_overlap <- function(gene_lists, metric) {
  transition_names <- names(gene_lists)
  
  if (length(gene_lists) < 2) {
    return(data.frame(
      metric = character(0),
      comparison = character(0),
      n_genes = integer(0),
      genes = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  pairs <- combn(transition_names, 2, simplify = FALSE)
  
  do.call(
    rbind,
    lapply(pairs, function(pair) {
      genes_i <- intersect(gene_lists[[pair[1]]], gene_lists[[pair[2]]])
      
      data.frame(
        metric = metric,
        comparison = paste(pair, collapse = "_AND_"),
        n_genes = length(genes_i),
        genes = paste(genes_i, collapse = ";"),
        stringsAsFactors = FALSE
      )
    })
  )
}

transition_overlap_summary <- rbind(
  make_pairwise_overlap(hubs_symbols, "degree"),
  make_pairwise_overlap(betweenness_symbols, "betweenness"),
  make_pairwise_overlap(closeness_symbols, "closeness"),
  data.frame(
    metric = c("degree", "betweenness", "closeness", "combined"),
    comparison = c(
      "all_available_transitions",
      "all_available_transitions",
      "all_available_transitions",
      "degree_betweenness_closeness_recurrent"
    ),
    n_genes = c(
      length(recurrent_hubs),
      length(recurrent_betweenness),
      length(recurrent_closeness),
      length(top_top)
    ),
    genes = c(
      paste(recurrent_hubs, collapse = ";"),
      paste(recurrent_betweenness, collapse = ";"),
      paste(recurrent_closeness, collapse = ";"),
      paste(top_top, collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
)

cat("\nTransition overlap summary:\n")
print(transition_overlap_summary)



#NSCLC e OVARIAN

hubs_common_degree <- Reduce(
  intersect,
  list(
    transition_results[["1_vs_2"]]$key_genes$hubs_symbols,
    transition_results[["2_vs_3"]]$key_genes$hubs_symbols,
    transition_results[["3_vs_4"]]$key_genes$hubs_symbols
  )
)

print(hubs_common_degree)


#CRC

hubs_common_degree <- Reduce(
  intersect,
  list(
    transition_results[["1_vs_2"]]$key_genes$hubs_symbols,
    transition_results[["2_vs_3"]]$key_genes$hubs_symbols
  )
)

print(hubs_common_degree)



# GBM

hubs_common_degree <- Reduce(
  intersect,
  list(
    transition_results[["0_vs_1"]]$key_genes$hubs_symbols,
    transition_results[["1_vs_2"]]$key_genes$hubs_symbols,
    transition_results[["2_vs_3"]]$key_genes$hubs_symbols
  )
)

print(hubs_common_degree)


##############################################################################àà
##############################################################################àà
##############################################################################àà

# GBM e OVARIAN

deg_with_symbols <- as.data.frame(deg_contr_vs_disease)

deg_with_symbols$gene_id <- rownames(deg_with_symbols)
 
deg_with_symbols$Symbol <- gene_annotation[
         deg_with_symbols$gene_id,
         "Symbol"]

deg_with_symbols <- deg_with_symbols[ ,
       c(
           "gene_id",
             "Symbol",
             "baseMean",
             "log2FoldChange",
             "lfcSE",
             "stat",
             "pvalue",
             "padj"
         )]


print(deg_with_symbols)


##############################################################################àà
##############################################################################àà
##############################################################################àà


# CRC e NSCLC

deg_with_symbols <- as.data.frame(deg_contr_vs_disease)

# ID Ensembl presenti nelle righe dei risultati DESeq2
deg_with_symbols$gene_id <- rownames(deg_with_symbols)

# Rimuove eventuali versioni, ad esempio ENSG000001234.5
deg_with_symbols$gene_id_clean <- sub(
  "\\..*$",
  "",
  deg_with_symbols$gene_id
)

# Pulisce anche gli ID Ensembl dell'annotazione
annotation_ensembl_clean <- sub(
  "\\..*$",
  "",
  gene_annotation$EnsemblGeneID
)

# Associa ogni ID Ensembl al relativo Symbol
deg_with_symbols$Symbol <- gene_annotation$Symbol[
  match(
    deg_with_symbols$gene_id_clean,
    annotation_ensembl_clean
  )
]

# Seleziona e ordina le colonne finali
deg_with_symbols <- deg_with_symbols[
  ,
  c(
    "gene_id",
    "Symbol",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )
]

print(deg_with_symbols)


##############################################################################àà
##############################################################################àà
##############################################################################àà


write.csv(
  deg_with_symbols,
  file = file.path("DEGs_with_symbols.csv"),
  row.names = FALSE
)

















