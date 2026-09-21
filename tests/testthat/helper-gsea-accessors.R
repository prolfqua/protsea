# geneInCategory() moved from DOSE to enrichit when DOSE 4.6 split the
# enrichment classes out into that package. Resolve it from whichever namespace
# of the installed stack exports it, so the test runs on both.
gene_in_category <- function(x) {
  for (ns in c("enrichit", "DOSE")) {
    if (
      requireNamespace(ns, quietly = TRUE) &&
        "geneInCategory" %in% getNamespaceExports(ns)
    ) {
      return(getExportedValue(ns, "geneInCategory")(x))
    }
  }
  testthat::skip("geneInCategory() is not exported by DOSE or enrichit")
}
