native_gsea <- function() {
  # Include a tie, both signs and a non-default exponent.
  ranks <- stats::setNames(c(8, 7, 7, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -6, -7, -8),
                          paste0("g", seq_len(16)))
  sets <- data.frame(term = rep(c("positive", "negative"), each = 4),
                     gene = names(ranks)[c(1, 2, 4, 6, 11, 13, 15, 16)])
  set.seed(42)
  suppressWarnings(clusterProfiler::GSEA(
    ranks, exponent = 1.5, minGSSize = 2, maxGSSize = 10,
    pvalueCutoff = 1, TERM2GENE = sets, verbose = FALSE, seed = TRUE,
    nPermSimple = 100
  ))
}

test_that("native GSEA statistics and running scores survive JSON", {
  original <- native_gsea()
  expect_equal(nrow(original@result), 2)
  expect_true(any(original@result$NES > 0) && any(original@result$NES < 0))
  doc <- gsea_result_data(list(A_vs_B = original), "PTMSEA")
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA, na = "null")
  restored <- decode_gsea_json(json)$A_vs_B$PTMSEA
  expect_equal(restored@result, original@result)
  expect_equal(restored@geneList, original@geneList)
  expect_equal(restored@geneSets, original@geneSets)
  expect_equal(restored@params, original@params)
  gs_info <- utils::getFromNamespace("gsInfo", "enrichplot")
  for (id in original@result$ID) {
    expect_equal(gs_info(restored, id), gs_info(original, id))
  }
  expect_silent(enrichplot::gseaplot2(restored, geneSetID = original@result$ID[[1]]))
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  write_gsea_result_json(doc, path)
  expect_equal(decode_gsea_json(path)$A_vs_B$PTMSEA@result, original@result)
  raw <- attr(read_gsea_json(path), "raw_data")$A_vs_B
  expect_equal(build_gseaResult(raw$cats$PTMSEA, raw$gene_pool, raw$rank_list)@result,
               original@result)
})

test_that("order is taken from rank positions, not JSON object key order", {
  original <- native_gsea()
  # Data frames restored from MuData may have automatic integer row names.
  rownames(original@result) <- NULL
  doc <- gsea_result_data(list(A = original), "test")
  doc$data$A$gene_pool <- rev(doc$data$A$gene_pool)
  doc$rank_lists$A$entries <- rev(doc$rank_lists$A$entries)
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA)
  restored <- decode_gsea_json(json)$A$test
  expect_equal(restored@geneList, original@geneList)
  expect_equal(restored@result, original@result)
})

test_that("empty results and readable identifiers round-trip", {
  original <- native_gsea()
  original@result <- original@result[FALSE, ]
  original@gene2Symbol <- stats::setNames(paste0("symbol", seq_along(original@geneList)),
                                         names(original@geneList))
  original@readable <- TRUE
  doc <- gsea_result_data(list(A = original), "test")
  restored <- decode_gsea_json(jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA))$A$test
  expect_equal(restored@result, original@result)
  expect_equal(restored@gene2Symbol, original@gene2Symbol)
  expect_true(restored@readable)
})

test_that("STRING conversions do not invent native statistics", {
  doc <- read_gsea_json(system.file("extdata", "WU2848501_gsea_result.json.gz", package = "protsea"))
  raw <- attr(doc, "raw_data")[[1]]
  gr <- build_gseaResult(raw$cats[[1]], raw$gene_pool, raw$rank_list)
  expect_true(all(is.na(gr@result$NES)))
  expect_true(all(is.na(gr@result$pvalue)))
  expect_error(decode_gsea_json('{"data":{"A":{"contrast":"A","categories":{"X":{}}}}}'),
               "Native gsea_result")
})
