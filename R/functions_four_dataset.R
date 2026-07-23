###############################################################################
##################### FUNCTIONS: FOUR DATASET NETWORK ANALYSIS ################
###############################################################################

###############################################################################
# 1. Clean GEO titles - GBM / GSE156902
###############################################################################

clean_titles <- function(meta) {
  
  meta$title_clean <- sub("^Blood_Platelets_(.*?) \\[.*$", "\\1", meta$title)
  meta$title_clean <- sub("-TR\\d+$", "", meta$title_clean)
  
  idx_vumc <- grepl("^\\d", meta$title_clean) &
    grepl("VUMC", meta$title, ignore.case = TRUE)
  
  meta$title_clean[idx_vumc] <- paste0("VU", meta$title_clean[idx_vumc])
  
  replace_from_bracket <- function(keyword,
                                   pattern = ".*\\[([^]]+)\\].*",
                                   subpattern = NULL) {
    
    idx <- grepl(keyword, meta$title, ignore.case = TRUE)
    val <- sub(pattern, "\\1", meta$title[idx])
    
    if (!is.null(subpattern)) {
      val <- sub(subpattern, "\\1", val, ignore.case = TRUE)
    }
    
    val <- sub("-TR\\d+$", "", val)
    meta$title_clean[idx] <<- val
  }
  
  replace_from_bracket("Vienna")
  replace_from_bracket("VUMC-HC")
  replace_from_bracket("HD-NKI")
  replace_from_bracket(
    "DOT-HD",
    subpattern = "^VUMC-(.*?)-DOT-HD.*$"
  )
  
  meta
}


###############################################################################
# 2. Prepare GBM metadata - GSE156902
###############################################################################

prepare_gbm_metadata <- function(info_pat) {
  
  GBM <- info_pat[
    info_pat$`Patient group` %in% c("GBM", "AC"),
  ]
  
  GBM <- GBM[
    c(
      "Patient group",
      "Age",
      "Gender",
      "Baseline sample ID",
      "Sample name"
    )
  ]
  
  GBM <- GBM[!duplicated(GBM$`Sample name`), ]
  
  mentioned_in_baseline <- unique(
    GBM$`Sample name`[
      GBM$`Sample name` %in%
        GBM$`Baseline sample ID`[!is.na(GBM$`Baseline sample ID`)]
    ]
  )
  
  GBM$stage <- ifelse(
    GBM$`Baseline sample ID` == "NA" &
      GBM$`Sample name` %in% mentioned_in_baseline,
    0,
    -1
  )
  
  GBM$stage <- as.integer(
    sub(".*[Tt](\\d+)$", "\\1", GBM$`Sample name`)
  )
  
  GBM$stage <- ifelse(is.na(GBM$stage), -1, GBM$stage)
  
  GBM$stage <- ifelse(
    GBM$stage == -1 & GBM$`Patient group` == "GBM",
    0,
    GBM$stage
  )
  
  GBM
}


###############################################################################
# 3. Build GEO metadata table - GBM / GSE156902
###############################################################################

build_geo_metadata <- function(gse) {
  
  meta1 <- pData(gse[[1]])
  meta2 <- pData(gse[[2]])
  
  meta1 <- clean_titles(meta1)
  meta2 <- clean_titles(meta2)
  
  meta1 <- meta1[, intersect(colnames(meta1), colnames(meta2))]
  meta2 <- meta2[, intersect(colnames(meta1), colnames(meta2))]
  
  m <- rbind(meta1, meta2)
  
  m <- m[
    c(
      "title",
      "geo_accession",
      "title_clean"
    )
  ]
  
  if ("GSM4746640" %in% rownames(m)) {
    m["GSM4746640", ]$title_clean <- "VU609 T0"
  }
  
  m
}


###############################################################################
# 4. Merge clinical metadata with GEO accessions - GBM / GSE156902
###############################################################################

merge_metadata_with_geo <- function(GBM, geo_meta) {
  
  meta <- GBM %>%
    left_join(
      geo_meta %>% dplyr::select(title_clean, geo_accession),
      by = c("Sample name" = "title_clean")
    )
  
  meta <- meta %>%
    tibble::column_to_rownames(var = "geo_accession")
  
  meta
}


###############################################################################
# 5. Apply manual metadata corrections - GBM / GSE156902
###############################################################################

apply_manual_metadata_fixes <- function(meta) {
  
  if ("GSM4746493" %in% rownames(meta)) {
    meta["GSM4746493", ]$`Sample name` <- "Vumc-GBM-477t0"
  }
  
  if ("GSM4746495" %in% rownames(meta)) {
    meta["GSM4746495", ]$`Sample name` <- "Vumc-GBM-491t0"
  }
  
  if ("GSM4746484" %in% rownames(meta)) {
    meta["GSM4746484", ]$`Sample name` <- "Vumc-GBM-493t0"
  }
  
  if ("GSM4746637" %in% rownames(meta)) {
    meta["GSM4746637", ]$`Baseline sample ID` <- "Vumc-GBM-493t0"
  }
  
  if ("GSM4746506" %in% rownames(meta)) {
    meta["GSM4746506", ]$`Sample name` <- "Vumc-GBM-follow-up-VU501t0"
  }
  
  if ("GSM4746638" %in% rownames(meta)) {
    meta["GSM4746638", ]$`Baseline sample ID` <- "VU635 T0"
  }
  
  meta
}


###############################################################################
# 6. Prepare CRC metadata - Xu 2022
###############################################################################

prepare_crc_metadata <- function(info_pat, counts) {
  
  meta <- info_pat[
    info_pat$Group %in% c("hd", "stage1", "stage2", "stage3"),
  ]
  
  meta$`Sample ID` <- gsub(
    "-([[:alnum:]]{2,3})$",
    "\\1",
    meta$`Sample ID`
  )
  
  meta[meta$`Sample ID` == "4-2-1922690A", ]$`Sample ID` <- "1922690A"
  meta[meta$`Sample ID` == "4-3-1927708A", ]$`Sample ID` <- "1927708A"
  
  meta <- meta[
    meta$`Sample ID` %in% colnames(counts),
  ]
  
  meta <- meta %>%
    tibble::column_to_rownames(var = "Sample ID")
  
  meta$Stage <- ifelse(
    meta$Group == "hd",
    0,
    ifelse(
      meta$Group == "stage1",
      1,
      ifelse(meta$Group == "stage2", 2, 3)
    )
  )
  
  meta$Stage <- as.integer(meta$Stage)
  
  meta$Group <- ifelse(
    meta$Group == "hd",
    "HD",
    "CRC"
  )
  
  meta
}


###############################################################################
# 7. Prepare NSCLC metadata - Antunes-Ferreira / GSE207586
###############################################################################

convert_nsclc_stage <- function(stage, group) {
  
  stage_chr <- trimws(as.character(stage))
  group_chr <- trimws(as.character(group))
  
  stage_chr <- toupper(stage_chr)
  group_chr <- toupper(group_chr)
  
  stage_chr <- gsub("^STAGE\\s*", "", stage_chr)
  stage_chr <- gsub("[^A-Z0-9]", "", stage_chr)
  
  out <- rep(NA_integer_, length(stage_chr))
  
  out[group_chr == "CONTROL"] <- 0L
  
  out[stage_chr %in% c("0", "CONTROL", "CTRL", "HEALTHY", "HD")] <- 0L
  out[stage_chr %in% c("1", "I", "IA", "IB")] <- 1L
  out[stage_chr %in% c("2", "II", "IIA", "IIB")] <- 2L
  out[stage_chr %in% c("3", "III", "IIIA", "IIIB", "IIIC")] <- 3L
  out[stage_chr %in% c("4", "IV", "IVA", "IVB")] <- 4L
  
  out
}

prepare_nsclc_metadata <- function(info_pat, counts) {
  
  info_pat <- as.data.frame(info_pat)
  
  new_colnames <- gsub("\\.", "-", colnames(counts))
  new_colnames <- sub("-$", "", new_colnames)
  colnames(counts) <- new_colnames
  
  meta <- info_pat %>%
    tibble::column_to_rownames(var = "Sample ID")
  
  meta <- meta[
    ,
    c(
      "Patient Group",
      "Gender",
      "Age",
      "Smoking",
      "Cancer Stage"
    ),
    drop = FALSE
  ]
  
  colnames(meta) <- c(
    "Group",
    "Gender",
    "Age",
    "Smoking",
    "OriginalStage"
  )
  
  meta$Stage <- convert_nsclc_stage(
    stage = meta$OriginalStage,
    group = meta$Group
  )
  
  meta$Group <- as.character(meta$Group)
  meta$Group <- ifelse(meta$Group == "Control", "Control", "NSCLC")
  
  meta$Gender <- as.factor(meta$Gender)
  meta$Smoking <- as.factor(meta$Smoking)
  meta$Stage <- as.integer(meta$Stage)
  
  meta <- meta[
    rownames(meta) %in% colnames(counts),
    ,
    drop = FALSE
  ]
  
  list(
    counts = counts,
    meta = meta
  )
}




###############################################################################
# 8. Build GEO metadata table - Ovarian / GSE183635
###############################################################################

build_ovarian_geo_metadata <- function(gse) {
  
  meta1 <- pData(gse[[1]])
  meta2 <- pData(gse[[2]])
  
  common_cols <- intersect(
    colnames(meta1),
    colnames(meta2)
  )
  
  meta1 <- meta1[, common_cols, drop = FALSE]
  meta2 <- meta2[, common_cols, drop = FALSE]
  
  m <- rbind(meta1, meta2)
  
  m <- m[
    ,
    c(
      "title",
      "geo_accession"
    ),
    drop = FALSE
  ]
  
  m
}


###############################################################################
# 9. Prepare Ovarian metadata - GSE183635
###############################################################################

convert_ovarian_stage <- function(stage) {
  
  stage_chr <- trimws(as.character(stage))
  stage_chr <- toupper(stage_chr)
  stage_chr <- gsub("^STAGE\\s*", "", stage_chr)
  stage_chr <- gsub("[^A-Z0-9]", "", stage_chr)
  
  out <- rep(NA_integer_, length(stage_chr))
  
  out[stage_chr %in% c("NA", "N\\.A", "NA", "NAN", "NAA", "N.A.", "N.A", "N/A", "", "0", "CONTROL", "ASYMPTOMATICCONTROLS")] <- 0L
  out[stage_chr %in% c("I", "1", "IA", "IB", "IC")] <- 1L
  out[stage_chr %in% c("II", "2", "IIA", "IIB", "IIC")] <- 2L
  out[stage_chr %in% c("III", "3", "IIIA", "IIIB", "IIIC")] <- 3L
  out[stage_chr %in% c("IV", "4", "IVA", "IVB", "IVC")] <- 4L
  
  out
}

prepare_ovarian_metadata <- function(info_pat,
                                     counts,
                                     geo_meta) {
  
  info_pat <- as.data.frame(info_pat)
  
  info_pat <- info_pat %>%
    left_join(
      geo_meta %>% dplyr::select(title, geo_accession),
      by = c("Sample ID" = "title")
    )
  
  meta <- info_pat[
    info_pat$Group %in% c("Ovarian cancer", "Asymptomatic controls"),
    ,
    drop = FALSE
  ]
  
  meta <- meta[
    !(meta$Stage == "n.a." & meta$Group == "Ovarian cancer"),
    ,
    drop = FALSE
  ]
  
  meta <- meta[
    !(meta$Stage == "n.a." & meta$Sex != "F"),
    ,
    drop = FALSE
  ]
  
  meta$Stage <- convert_ovarian_stage(meta$Stage)
  
  meta <- meta[
    ,
    c(
      "Group",
      "Stage",
      "Age",
      "Sex",
      "geo_accession"
    ),
    drop = FALSE
  ]
  
  meta <- meta[
    !is.na(meta$geo_accession),
    ,
    drop = FALSE
  ]
  
  meta <- meta[
    !duplicated(meta$geo_accession),
    ,
    drop = FALSE
  ]
  
  meta <- as.data.frame(meta)
  rownames(meta) <- NULL
  
  meta <- meta %>%
    tibble::column_to_rownames(var = "geo_accession")
  
  meta$Age <- ifelse(
    meta$Age == "n.a.",
    NA,
    as.integer(as.character(meta$Age))
  )
  
  meta$Stage <- as.integer(meta$Stage)
  meta$Sex <- as.factor(meta$Sex)
  meta$Group <- as.character(meta$Group)
  
  meta <- meta[
    rownames(meta) %in% colnames(counts),
    ,
    drop = FALSE
  ]
  
  meta
}


###############################################################################
# 10. Align counts and metadata
###############################################################################

align_counts_metadata <- function(counts, meta) {
  
  valid_samples <- intersect(rownames(meta), colnames(counts))
  
  counts <- counts[, valid_samples, drop = FALSE]
  meta <- meta[valid_samples, , drop = FALSE]
  
  colnames(meta)[colnames(meta) == "Patient group"] <- "Group"
  colnames(meta)[colnames(meta) == "stage"] <- "Stage"
  
  meta$Stage <- as.integer(meta$Stage)
  
  list(
    counts = counts,
    meta = meta
  )
}


###############################################################################
# 9. Run DESeq2 contrast
###############################################################################

run_deseq2_contrast <- function(counts_mat,
                                metadata,
                                group_col = "Group",
                                contrast_vec,
                                min_count = 20,
                                min_prop = 0.8,
                                padj_cutoff = 0.05,
                                lfc_cutoff = 0.33) {
  
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "integer"
  
  metadata[[group_col]] <- factor(metadata[[group_col]])
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = metadata[, group_col, drop = FALSE],
    design = as.formula(paste0("~ ", group_col))
  )
  
  print(summary(rowMeans(counts(dds))))
  
  keep <- rowSums(counts(dds) >= min_count) >= min_prop * ncol(dds)
  
  cat("Genes retained:", sum(keep), "\n")
  cat("Fraction retained:", round(sum(keep) / length(keep), 3), "\n")
  
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = contrast_vec
  )
  
  res <- res[order(res$padj), ]
  
  deg <- subset(
    res,
    !is.na(padj) &
      padj < padj_cutoff &
      abs(log2FoldChange) > lfc_cutoff
  )
  
  list(
    dds = dds,
    res = res,
    deg = deg,
    keep = keep
  )
}


###############################################################################
# 10. Get VST matrix for DE genes
###############################################################################

get_vst_deg_matrix <- function(dds,
                               deg_table,
                               blind = FALSE) {
  
  vsd <- vst(dds, blind = blind)
  expr_vst <- assay(vsd)
  
  genes_to_keep <- intersect(
    rownames(deg_table),
    rownames(expr_vst)
  )
  
  expr_vst[
    genes_to_keep,
    ,
    drop = FALSE
  ]
}


###############################################################################
# 11. Select control and disease baseline/stage samples for DESeq2
###############################################################################

select_control_vs_baseline_disease <- function(meta, counts) {
  
  meta_contrast <- meta[
    meta$Stage %in% c(-1, 0),
    ,
    drop = FALSE
  ]
  
  valid_samples <- intersect(
    rownames(meta_contrast),
    colnames(counts)
  )
  
  counts_contrast <- counts[, valid_samples, drop = FALSE]
  meta_contrast <- meta_contrast[valid_samples, , drop = FALSE]
  
  list(
    counts = counts_contrast,
    meta = meta_contrast
  )
}

select_ac_vs_t0 <- function(meta, counts) {
  select_control_vs_baseline_disease(meta = meta, counts = counts)
}


###############################################################################
# 12. Select patients present from T0 to T3 - GBM longitudinal subset
###############################################################################

select_complete_t0_t3_patients <- function(meta) {
  
  stage_0 <- meta[meta$Stage == 0, , drop = FALSE]
  stage_1 <- meta[meta$Stage == 1, , drop = FALSE]
  stage_2 <- meta[meta$Stage == 2, , drop = FALSE]
  stage_3 <- meta[meta$Stage == 3, , drop = FALSE]
  
  meta$t0_t3 <- ifelse(
    meta$`Sample name` %in% stage_1$`Baseline sample ID` &
      meta$`Sample name` %in% stage_2$`Baseline sample ID` &
      meta$`Sample name` %in% stage_3$`Baseline sample ID`,
    1,
    0
  )
  
  meta$t1_t3 <- ifelse(
    meta$Stage == 1 &
      meta$`Baseline sample ID` %in% stage_2$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_3$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_0$`Sample name`,
    1,
    0
  )
  
  meta$t2_t3 <- ifelse(
    meta$Stage == 2 &
      meta$`Baseline sample ID` %in% stage_3$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_1$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_0$`Sample name`,
    1,
    0
  )
  
  meta$t3_t3 <- ifelse(
    meta$Stage == 3 &
      meta$`Baseline sample ID` %in% stage_2$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_1$`Baseline sample ID` &
      meta$`Baseline sample ID` %in% stage_0$`Sample name`,
    1,
    0
  )
  
  meta$in_all_stages <- ifelse(
    meta$t0_t3 + meta$t1_t3 + meta$t2_t3 + meta$t3_t3 > 0,
    1,
    0
  )
  
  meta[
    meta$in_all_stages == 1,
    ,
    drop = FALSE
  ]
}


###############################################################################
# 13. Select stage 0-3 samples for network analysis
###############################################################################

select_stage_samples <- function(meta, stages) {
  
  meta[
    meta$Stage %in% stages,
    ,
    drop = FALSE
  ]
}

# Kept as a compatibility wrapper for existing scripts.
select_stage_0_to_3_samples <- function(meta) {
  select_stage_samples(meta, stages = 0:3)
}


###############################################################################
# 14. Get normalized matrix for network analysis
###############################################################################

get_normalized_deg_matrix <- function(counts_mat,
                                      metadata,
                                      deg_genes) {
  
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "integer"
  
  dds_all <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = metadata,
    design = ~ 1
  )
  
  dds_all <- estimateSizeFactors(dds_all)
  
  norm_counts_all <- counts(
    dds_all,
    normalized = TRUE
  )
  
  genes_to_keep <- intersect(
    deg_genes,
    rownames(norm_counts_all)
  )
  
  norm_counts_all[
    genes_to_keep,
    ,
    drop = FALSE
  ]
}


###############################################################################
# 15. Create stage expression list
###############################################################################

make_stage_expr_list <- function(expr_mat,
                                 metadata,
                                 stage_col = "Stage") {
  
  stages <- sort(unique(metadata[[stage_col]]))
  stages <- stages[!is.na(stages)]
  
  stage_expr_list <- lapply(stages, function(stage_i) {
    
    samples_i <- rownames(metadata)[metadata[[stage_col]] == stage_i]
    
    expr_mat[
      ,
      samples_i,
      drop = FALSE
    ]
  })
  
  names(stage_expr_list) <- as.character(stages)
  
  stage_expr_list
}


###############################################################################
# 16. Compute stage correlations
###############################################################################

compute_stage_correlations <- function(stage_expr_list,
                                       method = "spearman") {
  
  lapply(stage_expr_list, function(mat) {
    cor(
      t(mat),
      method = method
    )
  })
}


###############################################################################
# 17. Fisher z-transform
###############################################################################

z_transform <- function(rho) {
  
  rho[rho >= 1] <- 0.999999
  rho[rho <= -1] <- -0.999999
  
  0.5 * log((1 + rho) / (1 - rho))
}


###############################################################################
# 18. Build differential correlation network
###############################################################################

build_diff_network <- function(cor_mat_a,
                               cor_mat_b,
                               n_a,
                               n_b,
                               z_cutoff = 3) {
  
  if (n_a <= 3 || n_b <= 3) {
    stop("Fisher z differential correlation requires more than 3 samples per stage.")
  }
  
  z_a <- z_transform(cor_mat_a)
  z_b <- z_transform(cor_mat_b)
  
  z_diff <- (z_a - z_b) / sqrt(
    1 / (n_a - 3) + 1 / (n_b - 3)
  )
  
  z_diff_sym <- (z_diff + t(z_diff)) / 2
  
  adj_diff <- ifelse(
    abs(z_diff_sym) >= z_cutoff,
    1,
    0
  )
  
  diag(adj_diff) <- 0
  adj_diff[is.na(adj_diff)] <- 0
  
  igraph::graph_from_adjacency_matrix(
    adj_diff,
    mode = "upper"
  )
}


###############################################################################
# 19. Map graph node IDs to gene symbols
###############################################################################

map_ids_to_symbols <- function(ids,
                               gene_annotation) {
  
  symbols <- gene_annotation$Symbol[
    match(ids, rownames(gene_annotation))
  ]
  
  if ("EnsemblGeneID" %in% colnames(gene_annotation)) {
    idx_na <- is.na(symbols)
    
    symbols[idx_na] <- gene_annotation$Symbol[
      match(ids[idx_na], gene_annotation$EnsemblGeneID)
    ]
  }
  
  symbols <- ifelse(
    is.na(symbols),
    ids,
    symbols
  )
  
  unique(symbols)
}


###############################################################################
# 20. Extract hubs and central genes
###############################################################################

extract_network_key_genes <- function(graph,
                                      gene_annotation,
                                      percentile = 0.90) {
  
  deg <- degree(graph)
  deg_nonzero <- sort(
    deg[deg > 0],
    decreasing = TRUE
  )
  
  if (length(deg_nonzero) > 0) {
    hubs <- deg_nonzero[
      deg_nonzero >= quantile(
        deg_nonzero,
        percentile,
        na.rm = TRUE
      )
    ]
  } else {
    hubs <- deg_nonzero
  }
  
  bet <- betweenness(
    graph,
    directed = FALSE
  )
  
  top_bet <- bet[
    bet >= quantile(
      bet,
      percentile,
      na.rm = TRUE
    )
  ]
  
  cen <- closeness(graph)
  
  top_cen <- cen[
    cen >= quantile(
      cen,
      percentile,
      na.rm = TRUE
    )
  ]
  
  hubs_symbols <- map_ids_to_symbols(
    names(hubs),
    gene_annotation
  )
  
  bet_symbols <- map_ids_to_symbols(
    names(top_bet),
    gene_annotation
  )
  
  cen_symbols <- map_ids_to_symbols(
    names(top_cen),
    gene_annotation
  )
  
  genes_intersection <- intersect(
    intersect(hubs_symbols, bet_symbols),
    cen_symbols
  )
  
  list(
    degree = deg,
    degree_nonzero = deg_nonzero,
    hubs = hubs,
    betweenness = bet,
    top_betweenness = top_bet,
    closeness = cen,
    top_closeness = top_cen,
    hubs_symbols = hubs_symbols,
    betweenness_symbols = bet_symbols,
    closeness_symbols = cen_symbols,
    genes_intersection = genes_intersection
  )
}


###############################################################################
# 21. Analyze one stage transition
###############################################################################

analyze_stage_transition <- function(stage_a,
                                     stage_b,
                                     cor_matrices,
                                     metadata,
                                     gene_annotation,
                                     stage_col = "Stage",
                                     z_cutoff = 3,
                                     percentile = 0.90) {
  
  stage_a_chr <- as.character(stage_a)
  stage_b_chr <- as.character(stage_b)
  
  if (!stage_a_chr %in% names(cor_matrices)) {
    stop(paste0("Missing correlation matrix for stage ", stage_a_chr))
  }
  
  if (!stage_b_chr %in% names(cor_matrices)) {
    stop(paste0("Missing correlation matrix for stage ", stage_b_chr))
  }
  
  n_a <- nrow(
    metadata[
      metadata[[stage_col]] == stage_a,
      ,
      drop = FALSE
    ]
  )
  
  n_b <- nrow(
    metadata[
      metadata[[stage_col]] == stage_b,
      ,
      drop = FALSE
    ]
  )
  
  graph <- build_diff_network(
    cor_mat_a = cor_matrices[[stage_a_chr]],
    cor_mat_b = cor_matrices[[stage_b_chr]],
    n_a = n_a,
    n_b = n_b,
    z_cutoff = z_cutoff
  )
  
  key_genes <- extract_network_key_genes(
    graph = graph,
    gene_annotation = gene_annotation,
    percentile = percentile
  )
  
  comp <- components(graph)
  
  summary <- data.frame(
    transition = paste0(stage_a, "_vs_", stage_b),
    nodes = gorder(graph),
    edges = gsize(graph),
    density = edge_density(graph),
    connected = is_connected(graph),
    n_components = comp$no,
    lcc_size = max(comp$csize),
    n_hubs = length(key_genes$hubs_symbols),
    n_top_betweenness = length(key_genes$betweenness_symbols),
    n_top_closeness = length(key_genes$closeness_symbols),
    n_intersection_genes = length(key_genes$genes_intersection),
    stringsAsFactors = FALSE
  )
  
  list(
    graph = graph,
    key_genes = key_genes,
    summary = summary
  )
}


###############################################################################
# 22. Analyze all consecutive transitions
###############################################################################

analyze_consecutive_transitions <- function(cor_matrices,
                                            metadata,
                                            gene_annotation,
                                            transitions = list(
                                              c(0, 1),
                                              c(1, 2),
                                              c(2, 3)
                                            ),
                                            z_cutoff = 3,
                                            percentile = 0.90) {
  
  results <- lapply(transitions, function(pair) {
    
    analyze_stage_transition(
      stage_a = pair[1],
      stage_b = pair[2],
      cor_matrices = cor_matrices,
      metadata = metadata,
      gene_annotation = gene_annotation,
      z_cutoff = z_cutoff,
      percentile = percentile
    )
  })
  
  names(results) <- sapply(
    transitions,
    function(pair) paste0(pair[1], "_vs_", pair[2])
  )
  
  results
}


###############################################################################
# 23. Extract transition summaries
###############################################################################

summarise_transition_results <- function(transition_results) {
  
  bind_rows(
    lapply(
      transition_results,
      function(x) x$summary
    )
  )
}


###############################################################################
# 24. Save key genes
###############################################################################

save_transition_genes <- function(transition_results,
                                  output_dir) {
  
  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  for (nm in names(transition_results)) {
    
    genes <- transition_results[[nm]]$key_genes$genes_intersection
    
    write.csv(
      data.frame(genes = genes),
      file = file.path(output_dir, paste0("hubs_", nm, ".csv")),
      row.names = FALSE
    )
  }
}