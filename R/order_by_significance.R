#' Order Enrichment Result Tables as clusterProfiler Does
#'
#' enrichplot's `dotplot()`, `ridgeplot()` and the other `showCategory` views
#' take the first rows of `@result`, relying on the order clusterProfiler
#' returns. Documents written by other programs, such as gseapy-based MEA or
#' STRING, carry no ordering guarantee, so every decoded result is put into
#' clusterProfiler's order here: GSEA by adjusted p-value then descending
#' absolute NES (`DOSE::GSEA`), over-representation by p-value (`enricher`).
#' Missing values sort last.
#'
#' @param result The `@result` data frame of an enrichment object.
#' @return `result` with its rows reordered.
#' @name order_by_significance
#' @keywords internal
NULL

#' @rdname order_by_significance
order_gsea_result <- function(result) {
  result[order(result$p.adjust, -abs(result$NES), na.last = TRUE), , drop = FALSE]
}

#' @rdname order_by_significance
order_enrich_result <- function(result) {
  result[order(result$pvalue, result$p.adjust, na.last = TRUE), , drop = FALSE]
}
