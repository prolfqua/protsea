# Quarto rendering of a completed STRING-GSEA result directory.
#
# The STRING report source is a runtime template under `inst/templates/`. Its
# FGCZ styling comes from a directory-level `_metadata.yml` and its
# Find/Download toolbar is wired via `include-after-body: fgcz-plot-finder.html`.
# `data-raw/sync_quarto_assets.R` keeps vignette build-time assets synchronized
# with `fgczQuartoTemplate`.
#
# Rendering does not use those installed copies: it calls
# `fgczQuartoTemplate::fgcz_render()`, which stages the assets next to the qmd
# from the installed `fgczQuartoTemplate` before rendering, so a runtime render
# picks up the template the way `prolfquapp` does -- one package owns the theme.
# `buttons = FALSE`: the reports already wire the toolbar themselves, so
# `fgcz_render()` must not inject it a second time.

.fgcz_render_one <- function(qmd_source, output_dir, execute_params) {
  if (!file.exists(qmd_source)) {
    stop("Quarto report source not found: ", qmd_source, call. = FALSE)
  }
  staged <- file.path(output_dir, basename(qmd_source))
  if (!isTRUE(file.copy(qmd_source, staged, overwrite = TRUE))) {
    stop("Could not stage ", basename(qmd_source), " into ", output_dir, call. = FALSE)
  }
  # fgcz_render() copies these itself; calling it here names them for cleanup.
  assets <- fgczQuartoTemplate::fgcz_copy_assets(staged)
  on.exit(unlink(c(staged, assets)), add = TRUE)

  fgczQuartoTemplate::fgcz_render(
    input = staged,
    buttons = FALSE,
    execute_params = execute_params
  )
  rendered <- sub("[.]qmd$", ".html", staged)
  if (!file.exists(rendered)) {
    stop("Quarto render produced no HTML for ", basename(qmd_source), call. = FALSE)
  }
  # Quarto leaves an intermediate .rmarkdown beside a report it knits.
  unlink(sub("[.]qmd$", ".rmarkdown", staged))
  invisible(rendered)
}

#' Render the STRING-GSEA reports for one result directory
#'
#' Renders the enrichment report and the result landing page into
#' `output_dir` through [fgczQuartoTemplate::fgcz_render()], which applies the
#' FGCZ theme, banner and toolbar from the installed `fgczQuartoTemplate`. The
#' staged sources and template assets are removed again, leaving only the
#' rendered HTML.
#'
#' @param output_dir Result directory to render in place; holds the GSEA result
#'   JSON the report reads.
#' @param workunit_id B-Fabric workunit identifier, which names the result JSON
#'   and is shown on the landing page.
#' @param doc_dir Directory holding `GSEA_report.qmd`, normally
#'   `system.file("templates", package = "protsea")`.
#' @param templates_dir Directory holding `index.qmd`, normally
#'   `system.file("templates", package = "protsea")`.
#' @return Invisibly, the paths of the rendered HTML files.
#' @export
render_gsea_reports <- function(output_dir, workunit_id, doc_dir, templates_dir) {
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  report <- .fgcz_render_one(
    file.path(doc_dir, "GSEA_report.qmd"),
    output_dir,
    list(json_path = sprintf("WU%s_gsea_result.json", workunit_id))
  )
  index <- .fgcz_render_one(
    file.path(templates_dir, "index.qmd"),
    output_dir,
    list(workunit_id = workunit_id, package_dir = ".")
  )
  invisible(c(report, index))
}
