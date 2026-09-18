# protsea development

Standalone R package for shared enrichment JSON and clusterProfiler conversions, extracted from stringdbpy. See README.md for the JSON contract and NEWS.md for changes.

| Task | Command |
|---|---|
| Generate roxygen documentation | `make document` |
| Run unit tests | `make test` |
| Check the package without rebuilding vignettes | `make check-fast` |
| Build and install with vignettes | `make install` |
| Build vignettes | `make build-vignettes` |

- Edit R code in R/ and tests in tests/testthat/. Generate NAMESPACE and man/ with roxygen; do not edit them manually.
- Record behavior changes in NEWS.md under the current DESCRIPTION version.
- Keep the existing STRING JSON fields compatible. Preserve actual native GSEA statistics, ranking order, identifiers, gene sets and exponent; never manufacture NES from FDR.
- Verify native result → JSON → reconstructed result, including numeric running-score curves, when changing the codec. Coordinate persisted format changes with stringdbpy and prophosqua.
- Edit report sources only in vignettes/. Render changed QMDs and provide the generated HTML path. PTM-specific report templates belong in prophosqua.
- The Makefile follows the ecosystem template. Package-specific targets belong in make/package.mk.
