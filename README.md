# protsea

[Documentation and vignette](https://prolfqua.github.io/protsea/)

Shared enrichment JSON and clusterProfiler plotting conversions, extracted from `stringdbpy/stringGSEAplot` (source commit `7032a69`). Install with `make install`; run tests with `make test`. Report QMD sources live in `vignettes/` and are built and installed with the package.

## Native GSEA round trip

```r
document <- protsea::gsea_result_data(
  list(treatment_vs_control = result), # clusterProfiler::GSEA() result
  category = "PTMSEA"
)
json <- as.character(jsonlite::toJSON(document, auto_unbox = TRUE, digits = NA, na = "null"))
# Store json directly in MuData uns; no intermediate file is required.
restored <- protsea::decode_gsea_json(json)$treatment_vs_control$PTMSEA
enrichplot::gseaplot2(restored, geneSetID = 1)
```

For standalone file exchange, use `write_gsea_result_json(document, path)` and `decode_gsea_json(path)`. Existing STRING files remain readable through `read_gsea_json(path)`, returning nested `enrichResult` objects; `build_compareClusterResult()` supports cross-contrast plots.

## Existing format and additive native extension

The common document retains `data[contrast]` (a shared `gene_pool` and `categories`) and `rank_lists[contrast].entries`. STRING `metadata` and `links` are optional. Term fields remain compatible with stringdbpy and PTM consumers. `leading_edge_ids` carries the leading subset separately from full mapped `gene_ids`.

For clusterProfiler results, each category additionally contains `gsea_result`:

| Field | Content |
|---|---|
| `result.columns`, `result.types`, `result.row_names` | Original result table, including ES, NES, raw p-value, adjusted p-value, q-value, leading edge and any backend-specific columns; empty columns retain their types. |
| `gene_sets` | Complete original gene sets, including members absent from the ranking. |
| `params`, `param_types` | Original analysis parameters and their R storage types, including the actual weighting exponent. |
| `running_scores` | Per-term running enrichment-score vectors in ranked-list order. |
| `hit_indices` | Per-term one-based positions of gene-set members in the ranked list. |
| `organism`, `set_type`, `key_type`, `readable`, `gene2symbol` | Identifier and collection metadata required by the reconstructed object. |

Every ranked input is in the pool. Its integer `rank` determines the exact original order, including ties; `input_label` addresses the original identifier in `rank_lists.entries`. Object-key ordering and display labels are never used to determine the native ranking. PTM-specific canonical pool identifiers may differ from `input_label` without losing the submitted identifiers.

The original `enrichment_score` field keeps its existing producer-specific meaning: STRING's score for STRING output and NES for native PTM/GSEA output. Exact native statistics are unambiguous in `gsea_result.result.columns.enrichmentScore`, `NES`, `pvalue`, and `p.adjust`. Consumers needing those statistics should decode that block instead of interpreting the legacy summary field.

Reconstruction restores `gseaResult` without re-running enrichment. The common native block records source running scores and hit positions for direct comparison, while `enrichplot` can independently reproduce them from the retained ranks, gene sets and exponent. Permutation matrices are omitted; the exchange format supports reporting, not resuming permutation calculations. JSON uses full numeric precision and `null` for missing table values.

STRING documents without a native block remain valid common documents. STRING does not provide native fgsea NES or raw p-values; its plotting adapter leaves those statistics missing. The STRING ridgeplot remains a descriptive view of mapped members, and `decode_gsea_json()` rejects documents without the native inputs rather than inventing them. Native MEA producers use the same block as fgsea producers and retain the exact GSEApy running scores and hit positions.

The STRING report can be rendered from the installed package through `render_gsea_reports()`. PTM report templates remain in prophosqua.
