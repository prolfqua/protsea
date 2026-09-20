.gsea_running_trace <- function(ranks, members, exponent) {
  hits <- names(ranks) %in% members
  if (!any(hits) || all(hits)) {
    stop("A running enrichment score requires both hits and misses")
  }
  weights <- abs(ranks)^exponent
  increments <- ifelse(
    hits,
    weights / sum(weights[hits]),
    -1 / sum(!hits)
  )
  list(
    running_scores = as.list(unname(cumsum(increments))),
    hit_indices = as.list(which(hits))
  )
}

#' Convert clusterProfiler GSEA Results to Shared JSON Data
#'
#' Extends the STRING `data` / `rank_lists` structure with a category-level
#' `gsea_result` block containing native statistics, gene sets and parameters.
#' The pool's integer `rank` preserves input order, including tied scores.
#' Every native block also stores the running enrichment score and one-based
#' hit positions for each term. These fields have the same representation for
#' every producer.
#' Identifiers are retained exactly as submitted; no biological ID mapping is
#' performed. Permutation matrices are not stored.
#'
#' @param results Named list of [DOSE::gseaResult-class] objects, one per contrast.
#' @param category Gene-set collection name.
#' @param method Algorithm label, e.g. `"fgsea"` or `"DOSE"`.
#' @return A list with `data` and `rank_lists`, suitable for MuData `uns` or
#'   [write_gsea_result_json()].
#' @export
gsea_result_data <- function(results, category, method = "fgsea") {
  if (is.null(names(results)) || any(!nzchar(names(results))) || anyDuplicated(names(results))) {
    stop("results must have unique, nonempty contrast names")
  }
  data <- lapply(names(results), function(contrast) {
    res <- results[[contrast]]
    stopifnot(methods::is(res, "gseaResult"))
    ranks <- res@geneList
    ids <- names(ranks)
    if (anyDuplicated(ids) || length(ids) != length(ranks) || any(!is.finite(ranks))) {
      stop("geneList must have unique identifiers and finite scores")
    }
    pool <- stats::setNames(lapply(seq_along(ranks), function(i) {
      list(protein_id = ids[[i]], label = ids[[i]], input_label = ids[[i]],
           input_value = unname(ranks[[i]]), rank = i)
    }), ids)
    terms <- lapply(seq_len(nrow(res@result)), function(i) {
      row <- res@result[i, , drop = FALSE]
      members <- intersect(res@geneSets[[row$ID]], ids)
      leading <- strsplit(row$core_enrichment, "/", fixed = TRUE)[[1]]
      if (res@readable) {
        leading <- names(res@gene2Symbol)[res@gene2Symbol %in% leading]
      }
      list(term_id = row$ID, category = category, description = row$Description,
           enrichment_score = row$NES,
           direction = if (row$NES > 0) "top" else if (row$NES < 0) "bottom" else "both ends",
           fdr = row$p.adjust, method = method,
           genes_mapped = length(members), genes_in_set = length(res@geneSets[[row$ID]]),
           gene_ids = as.list(members), leading_edge_ids = as.list(intersect(leading, members)))
    })
    traces <- lapply(res@result$ID, function(id) {
      .gsea_running_trace(ranks, res@geneSets[[id]], res@params$exponent)
    })
    names(traces) <- res@result$ID
    native <- list(
      result = list(columns = lapply(res@result, as.list),
                    types = as.list(vapply(res@result, typeof, character(1))),
                    row_names = as.list(rownames(res@result)),
                    row_name_type = typeof(attr(res@result, "row.names"))),
      gene_sets = lapply(res@geneSets, as.list), params = res@params,
      param_types = as.list(vapply(res@params, typeof, character(1))),
      organism = res@organism, set_type = res@setType, key_type = res@keytype,
      readable = res@readable, gene2symbol = as.list(res@gene2Symbol),
      running_scores = lapply(traces, `[[`, "running_scores"),
      hit_indices = lapply(traces, `[[`, "hit_indices")
    )
    categories <- stats::setNames(list(list(
      category = category, contrast = contrast, terms = terms, gsea_result = native
    )), category)
    list(contrast = contrast, gene_pool = pool, categories = categories)
  })
  names(data) <- names(results)
  rank_lists <- lapply(names(results), function(contrast) {
    list(contrast = contrast, entries = as.list(results[[contrast]]@geneList))
  })
  names(rank_lists) <- names(results)
  list(data = data, rank_lists = rank_lists)
}

#' Write Shared GSEA Data to JSON
#'
#' @param gsea_result Data returned by [gsea_result_data()].
#' @param path Output JSON file.
#' @return `path`, invisibly.
#' @export
write_gsea_result_json <- function(gsea_result, path) {
  jsonlite::write_json(gsea_result, path, auto_unbox = TRUE, digits = NA, na = "null")
  invisible(path)
}
