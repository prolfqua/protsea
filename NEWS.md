# protsea 0.1.0

- Add a separate PTM enrichment QMD for one analysis document containing PTM-SEA, Kinase GSEA, and MEA. The report reconstructs every plot from the shared JSON and uses 96 DPI without retina scaling to keep its self-contained HTML compact.
- Store per-term source running enrichment scores and one-based hit positions
  in the shared native GSEA block, including documents produced by GSEApy MEA,
  while reconstructing the same transient clusterProfiler object for plotting.
- Keep the single round-trip vignette visible in the generated website's Articles navigation.
- Add a Quarto vignette demonstrating an fgsea-backed clusterProfiler analysis, JSON serialization and reconstruction, equality checks, and plots made from the reconstructed result. Build package documentation with altdoc's Quarto website backend.
- Build the compact round-trip article as the package vignette. Ship `GSEA_report.qmd` as an installed runtime report template without executing its 112 MB example during every vignette or documentation build.
- Construct DOSE result classes from their registered definitions, supporting releases where those S4 classes are not namespace exports.
- Promote the STRING JSON readers, plotting conversions, and report templates from stringdbpy/stringGSEAplot into the standalone protsea package.
- Write clusterProfiler GSEA results to the existing JSON structure with native statistics, complete gene sets, ranking order, and parameters; reconstruct running-score plots after serialization.
- Stop deriving artificial NES values from STRING FDR values.
- Verify running-score round trips without relying on the private `DOSE::gseaScores` implementation, which moved to enrichit in Bioconductor 3.23.
