test_that("PTM report fixture contains every enrichment tool", {
  path <- system.file(
    "extdata", "PTM_DPA_gsea_result.json.gz",
    package = "protsea"
  )
  expect_true(nzchar(path))
  restored <- decode_gsea_json(path)

  expect_identical(names(restored), c("a_vs_b", "c_vs_b"))
  for (contrast in restored) {
    expect_identical(names(contrast), c("PTM-SEA", "KinaseLib", "MEA"))
    expect_true(all(vapply(contrast, methods::is, logical(1), "gseaResult")))
  }
})

test_that("PTM report renders shared JSON without rerunning enrichment", {
  path <- test_path("..", "..", "vignettes", "PTM_enrichment.qmd")
  skip_if_not(file.exists(path), "QMD source is not installed with the package")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  top_tabs <- c(
    "# Overview", "# PTM-SEA", "# Kinase GSEA", "# MEA",
    "# About methods", "# Session Info"
  )
  positions <- vapply(top_tabs, function(tab) regexpr(tab, source, fixed = TRUE)[[1]], integer(1))

  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
  expect_match(source, "fig.retina: 1", fixed = TRUE)
  expect_match(source, "protsea::decode_gsea_json", fixed = TRUE)
  expect_match(source, "enrichplot::dotplot", fixed = TRUE)
  expect_match(source, "enrichplot::ridgeplot", fixed = TRUE)
  expect_false(grepl("## Running score", source, fixed = TRUE))
  expect_false(grepl("clusterProfiler::GSEA(", source, fixed = TRUE))
  expect_false(grepl("fgsea::fgsea", source, fixed = TRUE))
})
