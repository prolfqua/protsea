# protsea 0.1.0

- Add a Quarto vignette demonstrating an fgsea-backed clusterProfiler analysis, JSON serialization and reconstruction, equality checks, and plots made from the reconstructed result. Build package documentation with altdoc's Quarto website backend.
- Build the compact round-trip article as the package vignette. Ship `GSEA_report.qmd` as an installed runtime report template without executing its 112 MB example during every vignette or documentation build.
- Install with current DOSE releases, where `compareClusterResult` is registered but is not an exported class.
- Promote the STRING JSON readers, plotting conversions, and report templates from stringdbpy/stringGSEAplot into the standalone protsea package.
- Write clusterProfiler GSEA results to the existing JSON structure with native statistics, complete gene sets, ranking order, and parameters; reconstruct running-score plots after serialization.
- Stop deriving artificial NES values from STRING FDR values.
