#' Build a gseaResult S4 object from one category's JSON data
#'
#' Constructs a [DOSE::gseaResult-class] suitable for `enrichplot::ridgeplot()`
#' and other GSEA-specific visualisations. Native `gsea_result` extensions
#' restore the original statistics. STRING-only results have missing ES/NES
#' and p-values; their ridgeplot is a descriptive view of mapped members.  The ranked gene list is taken from
#' the `input_value` field of `gene_pool`; term genes become `core_enrichment`.
#'
#' @param category_data List from JSON for one category (fields: category,
#'   contrast, terms).
#' @param gene_pool Named list from the JSON gene_pool for this contrast.
#' @param rank_list Named list from rank_lists for this contrast.
#' @return A [DOSE::gseaResult-class] object.
#' @importFrom methods new
#' @export
build_gseaResult <- function(category_data, gene_pool, rank_list) {
  if (!is.null(category_data$gsea_result)) {
    return(restore_gsea_result(category_data$gsea_result, gene_pool, rank_list))
  }
  terms <- category_data[["terms"]]
  if (length(terms) == 0) {
    stop("No terms found in category '", category_data[["category"]], "'")
  }

  # geneList: named numeric vector of input scores, sorted descending
  gene_values <- vapply(gene_pool, function(e) as.numeric(e[["input_value"]]), numeric(1))
  names(gene_values) <- vapply(gene_pool, function(e) e[["label"]], character(1))
  gene_list <- sort(gene_values, decreasing = TRUE)

  rows <- lapply(terms, function(term) {
    gene_ids <- unlist(term[["gene_ids"]])
    labels <- resolve_gene_ids(gene_ids, gene_pool)
    fdr_val <- as.numeric(term[["fdr"]])
    k <- as.integer(term[["genes_mapped"]])

    # STRING does not provide native GSEA ES/NES or unadjusted p-values.
    nes <- NA_real_

    gene_ranks <- match(labels, names(gene_list))
    med_rank <- as.integer(stats::median(gene_ranks[!is.na(gene_ranks)]))

    data.frame(
      ID = term[["term_id"]],
      Description = term[["description"]],
      setSize = k,
      enrichmentScore = nes,
      NES = nes,
      pvalue = NA_real_,
      p.adjust = fdr_val,
      qvalues = NA_real_,
      rank = med_rank,
      leading_edge = paste0("tags=", k, ", list=", length(gene_list)),
      core_enrichment = paste(labels, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, rows)
  rownames(result_df) <- result_df$ID

  gene_sets <- lapply(terms, function(term) {
    resolve_gene_ids(unlist(term[["gene_ids"]]), gene_pool)
  })
  names(gene_sets) <- vapply(terms, function(t) t[["term_id"]], character(1))

  methods::new(
    methods::getClass("gseaResult", where = asNamespace("DOSE")),
    result = result_df,
    organism = "unknown",
    setType = category_data[["category"]],
    geneSets = gene_sets,
    geneList = gene_list,
    keytype = "STRING",
    permScores = matrix(nrow = 0, ncol = 0),
    params = list(pvalueCutoff = 1.0, pAdjustMethod = "BH", minGSSize = 1L, maxGSSize = 500L, exponent = 1L),
    gene2Symbol = character(0),
    readable = TRUE
  )
}
