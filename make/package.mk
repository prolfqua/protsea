build build-vignettes test check-fast: sync-quarto-assets

.PHONY: sync-quarto-assets
sync-quarto-assets:
	Rscript data-raw/sync_quarto_assets.R
