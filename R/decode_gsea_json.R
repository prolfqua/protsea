#' Reconstruct Native GSEA Results from JSON
#'
#' Decodes a JSON string (including a string read from MuData `uns`) or a file
#' written with [write_gsea_result_json()]. Categories must contain the native
#' `gsea_result` extension. For STRING enrichment tables use [read_gsea_json()].
#'
#' @param json JSON text or a JSON file path.
#' @return Named nested list, indexed by contrast and category, of
#'   [DOSE::gseaResult-class] objects.
#' @export
decode_gsea_json <- function(json) {
  document <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  lapply(document$data, function(contrast) {
    lapply(contrast$categories, function(category) {
      restore_gsea_result(category$gsea_result, contrast$gene_pool,
                          document$rank_lists[[contrast$contrast]])
    })
  })
}

restore_gsea_result <- function(native, gene_pool, rank_list) {
  if (is.null(native)) {
    stop("Native gsea_result data is required to reconstruct clusterProfiler statistics")
  }
  columns <- lapply(names(native$result$columns), function(name) {
    values <- native$result$columns[[name]]
    values <- lapply(values, function(x) if (is.null(x)) NA else x)
    methods::as(unlist(values, use.names = FALSE), native$result$types[[name]])
  })
  names(columns) <- names(native$result$columns)
  result <- as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
  attr(result, "row.names") <- methods::as(
    unlist(native$result$row_names, use.names = FALSE), native$result$row_name_type
  )
  pool <- gene_pool[order(vapply(gene_pool, function(hit) hit$rank, numeric(1)))]
  ids <- vapply(pool, function(hit) hit$input_label, character(1))
  ranks <- vapply(rank_list$entries[ids], as.numeric, numeric(1))
  names(ranks) <- ids
  gene_sets <- lapply(native$gene_sets, function(x) as.character(unlist(x, use.names = FALSE)))
  symbols <- unlist(native$gene2symbol, use.names = TRUE)
  if (!length(symbols)) symbols <- character()
  methods::new(
    "gseaResult", result = result, geneList = ranks, geneSets = gene_sets,
    params = lapply(native$params, function(x) unlist(x, use.names = TRUE)),
    organism = native$organism, setType = native$set_type, keytype = native$key_type,
    readable = native$readable, gene2Symbol = symbols,
    permScores = matrix(nrow = 0, ncol = 0)
  )
}
