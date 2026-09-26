native_gsea <- function() {
  # Include a tie, both signs and a non-default exponent.
  ranks <- stats::setNames(c(8, 7, 7, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -6, -7, -8), paste0("g", seq_len(16)))
  sets <- data.frame(
    term = rep(c("positive", "negative"), each = 4),
    gene = names(ranks)[c(1, 2, 4, 6, 11, 13, 15, 16)]
  )
  set.seed(42)
  suppressWarnings(clusterProfiler::GSEA(
    ranks,
    exponent = 1.5,
    minGSSize = 2,
    maxGSSize = 10,
    pvalueCutoff = 1,
    TERM2GENE = sets,
    verbose = FALSE,
    seed = TRUE,
    nPermSimple = 100
  ))
}

expected_gsea_trace <- function(ranks, members, exponent) {
  hits <- names(ranks) %in% members
  weights <- abs(ranks)^exponent
  increments <- ifelse(
    hits,
    weights / sum(weights[hits]),
    -1 / sum(!hits)
  )
  data.frame(
    runningScore = cumsum(increments),
    position = as.integer(hits)
  )
}

test_that("native GSEA statistics and running scores survive JSON", {
  original <- native_gsea()
  expect_equal(nrow(original@result), 2)
  expect_true(any(original@result$NES > 0) && any(original@result$NES < 0))
  doc <- gsea_result_data(list(A_vs_B = original), "PTMSEA")
  native <- doc$data$A_vs_B$categories$PTMSEA$gsea_result
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA, na = "null")
  restored <- decode_gsea_json(json)$A_vs_B$PTMSEA
  expect_equal(restored@result, original@result)
  expect_equal(restored@geneList, original@geneList)
  expect_equal(restored@geneSets, original@geneSets)
  expect_equal(restored@params, original@params)
  expect_identical(restored@params, original@params)
  gs_info <- utils::getFromNamespace("gsInfo", "enrichplot")
  for (id in original@result$ID) {
    expect_equal(gs_info(restored, id), gs_info(original, id))
    source_trace <- expected_gsea_trace(
      original@geneList,
      original@geneSets[[id]],
      exponent = original@params$exponent
    )
    expect_equal(
      unlist(native$running_scores[[id]], use.names = FALSE),
      source_trace$runningScore
    )
    expect_equal(
      unlist(native$hit_indices[[id]], use.names = FALSE),
      which(source_trace$position == 1L)
    )
  }
  expect_silent(enrichplot::gseaplot2(restored, geneSetID = original@result$ID[[1]]))
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  write_gsea_result_json(doc, path)
  expect_equal(decode_gsea_json(path)$A_vs_B$PTMSEA@result, original@result)
  raw <- attr(read_gsea_json(path), "raw_data")$A_vs_B
  expect_equal(build_gseaResult(raw$cats$PTMSEA, raw$gene_pool, raw$rank_list)@result, original@result)
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
  expect_equal(as.list(restored@result), as.list(original@result))
  expect_identical(rownames(restored@result), restored@result$ID)
})

test_that("empty results and readable identifiers round-trip", {
  original <- native_gsea()
  original@result <- original@result[FALSE, ]
  original@gene2Symbol <- stats::setNames(paste0("symbol", seq_along(original@geneList)), names(original@geneList))
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
  expect_error(decode_gsea_json('{"data":{"A":{"contrast":"A","categories":{"X":{}}}}}'), "Native gsea_result")
})

test_that("GSEApy MEA uses the same native GSEA JSON structure", {
  path <- test_path("fixtures/gseapy-mea.json")
  document <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  native <- document$data$A_vs_B$categories$MEA$gsea_result
  restored <- decode_gsea_json(path)$A_vs_B$MEA

  expect_s4_class(restored, "gseaResult")
  expect_identical(
    names(restored@result),
    c(
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust",
      "qvalues",
      "rank",
      "leading_edge",
      "core_enrichment"
    )
  )
  expect_equal(restored@params$exponent, 1.5)
  expect_equal(names(restored@geneList), paste0("g", seq_len(16L)))
  expect_setequal(names(restored@geneSets), c("positive", "negative"))
  expect_identical(names(gene_in_category(restored)), restored@result$ID)
  expect_s3_class(enrichplot::ridgeplot(restored, showCategory = 2), "ggplot")

  for (id in restored@result$ID) {
    reproduced <- expected_gsea_trace(
      restored@geneList,
      restored@geneSets[[id]],
      exponent = restored@params$exponent
    )
    expect_equal(
      unlist(native$running_scores[[id]], use.names = FALSE),
      reproduced$runningScore,
      tolerance = 1e-12
    )
    expect_equal(
      unlist(native$hit_indices[[id]], use.names = FALSE),
      which(reproduced$position == 1L)
    )
  }
  expect_silent(enrichplot::gseaplot2(restored, geneSetID = 1L))
})

test_that("decoded native results are in clusterProfiler order whatever the document order", {
  original <- native_gsea()
  doc <- gsea_result_data(list(A_vs_B = original), "PTMSEA")
  native <- doc$data$A_vs_B$categories$PTMSEA$gsea_result
  reversed <- rev(seq_along(native$result$columns$ID))
  native$result$columns <- lapply(native$result$columns, function(x) x[reversed])
  native$result$row_names <- native$result$row_names[reversed]
  doc$data$A_vs_B$categories$PTMSEA$gsea_result <- native
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA, na = "null")
  restored <- decode_gsea_json(json)$A_vs_B$PTMSEA
  expect_equal(restored@result, original@result)
})

test_that("GSEA order is adjusted p-value then absolute NES, missing values last", {
  result <- data.frame(
    ID = c("na", "weak", "strong", "best"),
    p.adjust = c(NA, 0.05, 0.05, 0.01),
    NES = c(2, 1.2, -1.8, 1.1)
  )
  expect_equal(order_gsea_result(result)$ID, c("best", "strong", "weak", "na"))
})

test_that("over-representation order is p-value, missing values last", {
  result <- data.frame(ID = c("c", "a", "b"), pvalue = c(NA, 0.01, 0.001), p.adjust = c(0.2, 0.02, 0.01))
  expect_equal(order_enrich_result(result)$ID, c("b", "a", "c"))
})
