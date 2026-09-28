prepare_ora_inputs <- function(
  de,
  term2gene,
  direction,
  alpha,
  lfc_threshold
) {
  require_columns(
    de,
    c("ensembl_id", "padj", "log2fc_shrunk"),
    "differential-expression results"
  )
  require_columns(term2gene, c("term_id", "ensembl_id"), "term-to-gene mapping")
  if (!direction %in% c("up", "down")) {
    stop("ORA direction must be 'up' or 'down'", call. = FALSE)
  }
  mapped_ids <- unique(as.character(term2gene$ensembl_id))
  tested <- !is.na(de$padj) & de$ensembl_id %in% mapped_ids
  universe <- unique(as.character(de$ensembl_id[tested]))
  significant <- de$padj < alpha
  if (direction == "up") {
    significant <- significant & de$log2fc_shrunk > lfc_threshold
  } else {
    significant <- significant & de$log2fc_shrunk < -lfc_threshold
  }
  genes <- unique(as.character(de$ensembl_id[tested & significant]))
  list(genes = genes, universe = universe)
}

map_tested_genes <- function(de, term2gene, database) {
  require_columns(de, c("ensembl_id", "padj"), "differential-expression results")
  require_columns(term2gene, "ensembl_id", "term-to-gene mapping")
  tested <- unique(as.character(de$ensembl_id[!is.na(de$padj)]))
  mapped <- intersect(tested, unique(as.character(term2gene$ensembl_id)))
  list(
    mapped = mapped,
    audit = data.frame(
      database = database,
      total_de_rows = as.integer(nrow(de)),
      tested_genes = as.integer(length(tested)),
      mapped_tested_genes = as.integer(length(mapped)),
      unmapped_tested_genes = as.integer(length(setdiff(tested, mapped))),
      mapping_rate = if (length(tested) == 0) NA_real_ else length(mapped) / length(tested),
      stringsAsFactors = FALSE
    )
  )
}

empty_ora_result <- function(status, foreground_size, universe_size) {
  result <- data.frame(
    term_id = character(),
    term_name = character(),
    pvalue = numeric(),
    padj = numeric(),
    gene_count = integer(),
    foreground_size = integer(),
    universe_size = integer(),
    direction = character(),
    gene_ids = character(),
    stringsAsFactors = FALSE
  )
  attr(result, "status") <- status
  attr(result, "foreground_size") <- as.integer(foreground_size)
  attr(result, "universe_size") <- as.integer(universe_size)
  result
}

run_ora <- function(
  de,
  term2gene,
  direction,
  alpha,
  lfc_threshold,
  min_size = 5,
  max_size = 500
) {
  require_columns(
    term2gene,
    c("term_id", "term_name", "ensembl_id"),
    "term-to-gene mapping"
  )
  inputs <- prepare_ora_inputs(de, term2gene, direction, alpha, lfc_threshold)
  if (length(inputs$genes) == 0) {
    return(empty_ora_result(
      "no_significant_genes",
      length(inputs$genes),
      length(inputs$universe)
    ))
  }
  mapping <- unique(term2gene[, c("term_id", "ensembl_id", "term_name")])
  enrichment <- clusterProfiler::enricher(
    gene = inputs$genes,
    universe = inputs$universe,
    TERM2GENE = mapping[, c("term_id", "ensembl_id")],
    TERM2NAME = mapping[, c("term_id", "term_name")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = min_size,
    maxGSSize = max_size
  )
  if (is.null(enrichment)) {
    return(empty_ora_result(
      "no_testable_terms",
      length(inputs$genes),
      length(inputs$universe)
    ))
  }
  result <- as.data.frame(enrichment)
  if (nrow(result) == 0) {
    return(empty_ora_result(
      "no_testable_terms",
      length(inputs$genes),
      length(inputs$universe)
    ))
  }
  output <- data.frame(
    term_id = result$ID,
    term_name = result$Description,
    pvalue = result$pvalue,
    padj = result$p.adjust,
    gene_count = as.integer(result$Count),
    foreground_size = as.integer(length(inputs$genes)),
    universe_size = as.integer(length(inputs$universe)),
    direction = direction,
    gene_ids = result$geneID,
    stringsAsFactors = FALSE
  )
  attr(output, "status") <- "completed"
  attr(output, "foreground_size") <- as.integer(length(inputs$genes))
  attr(output, "universe_size") <- as.integer(length(inputs$universe))
  output
}

prepare_ranked_list <- function(de) {
  require_columns(de, c("ensembl_id", "stat"), "differential-expression results")
  keep <- !is.na(de$ensembl_id) & de$ensembl_id != "" & is.finite(de$stat)
  ranked <- de[keep, c("ensembl_id", "stat"), drop = FALSE]
  ranked <- ranked[order(ranked$stat, decreasing = TRUE), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$ensembl_id), , drop = FALSE]
  values <- as.numeric(ranked$stat)
  names(values) <- ranked$ensembl_id
  values
}

run_ranked_enrichment <- function(
  de,
  term2gene,
  seed,
  min_size = 10,
  max_size = 500
) {
  require_columns(term2gene, c("term_id", "ensembl_id"), "term-to-gene mapping")
  ranked <- prepare_ranked_list(de)
  pathways <- split(
    as.character(term2gene$ensembl_id),
    as.character(term2gene$term_id)
  )
  pathways <- lapply(pathways, unique)
  set.seed(seed)
  result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked,
    minSize = min_size,
    maxSize = max_size,
    eps = 0
  )
  result <- as.data.frame(result)
  if (nrow(result) == 0) return(result)
  term_names <- unique(term2gene[, c("term_id", "term_name")])
  result$term_name <- term_names$term_name[match(result$pathway, term_names$term_id)]
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  result[order(result$padj, -abs(result$NES)), , drop = FALSE]
}

build_go_bp_mapping <- function() {
  ensembl_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
  mapping <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_keys,
    columns = c("GOALL", "ONTOLOGYALL"),
    keytype = "ENSEMBL"
  ))
  mapping <- mapping[
    !is.na(mapping$GOALL) & mapping$ONTOLOGYALL == "BP",
    c("ENSEMBL", "GOALL"),
    drop = FALSE
  ]
  term_names <- AnnotationDbi::mapIds(
    GO.db::GO.db,
    keys = unique(mapping$GOALL),
    column = "TERM",
    keytype = "GOID",
    multiVals = "first"
  )
  result <- data.frame(
    term_id = mapping$GOALL,
    term_name = unname(term_names[mapping$GOALL]),
    ensembl_id = strip_ensembl_version(mapping$ENSEMBL),
    database_version = as.character(utils::packageVersion("GO.db")),
    stringsAsFactors = FALSE
  )
  unique(result[!is.na(result$term_name) & result$ensembl_id != "", ])
}

build_kegg_legacy_mapping <- function() {
  mapping <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  )
  result <- data.frame(
    term_id = as.character(mapping$gs_name),
    term_name = as.character(mapping$gs_description),
    ensembl_id = strip_ensembl_version(mapping$ensembl_gene),
    database_version = as.character(mapping$db_version),
    stringsAsFactors = FALSE
  )
  unique(result[!is.na(result$ensembl_id) & result$ensembl_id != "", ])
}
