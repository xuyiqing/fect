# Triple-difference manuscript figures 2–9

These scripts redraw the supplied three-case manuscript figures following the
3 November 2026 figure review. They do not rerun the saved 200/1,000-draw
bootstrap analyses or alter their confidence intervals.

## Package-wide plotting options

```r
plot(fit, preset = "publication")
plot(fit, type = "heterogeneous", covariate = moderator,
     preset = "publication", show.average = FALSE)
esplot(event_summary, preset = "publication")
```

The preset is opt-in. Existing defaults remain available. Titles and notes
belong in manuscript captions; numerical mean annotations are available using
`average.label = TRUE`. The companion `interflex` `tianzhu` branch supports
`plot(fit, preset = "publication")` as well.

## Exact manuscript redraw

Supply the original `case-materials-bundle-20260928` directory, containing saved
summary tables. No research data or large fitted objects are included here.

```sh
mkdir -p /path/to/figure-revision
Rscript tools/tripleDD/recover_inputs.R /path/to/case-materials-bundle-20260928 /path/to/figure-revision
Rscript tools/tripleDD/redraw_figures.R /path/to/case-materials-bundle-20260928 /path/to/figure-revision
```

Dependencies: `fect`, `ggplot2`, `patchwork`, `ragg`, `haven`, `digest`,
`jsonlite`, and Cairo support in R. Public source files are retrieved from
Harvard Dataverse and verified against original checksums. Fig. 4 uses the
saved GPS coefficients, without refitting the propensity model. Fig. 8 uses
a common FE fit with standard errors disabled, then the original loess formula.
The recovery script stops if the reconstructed inputs fail to match the saved
overlap or subgroup summaries. The validated recovery run used `fect` 2.2.0.

`figures/` contains eight combined vector PDFs and 600 dpi PNGs. `panels/`
contains the 13 separate manuscript assets with their original base names.
`captions.csv` identifies left/right panels and interval definitions outside
the plot. `inputs/recovered_input_validation.json` records source checksums
and numerical agreement.

## Figure-specific choices

- Fig. 2: quarter-width side margins, common y-scale, separate count strip.
- Fig. 3: side-by-side cue panels, a single count strip per panel, compact
  legend, explicit ticks at 0 and 1. The diamond at day 0 is a requested
  zero reference without a CI, not an observation. Former day-7 horizontal
  references are numerical annotations, not mislabeled postperiod means.
- Fig. 4: probability range 0.4–0.6; an assertion prevents hiding observations.
- Fig. 5: black estimates; original native fixed-gap t intervals retained.
- Fig. 6: exactly the clothing/non-clothing profiles, padded axes, uncapped CIs.
- Fig. 7: consistent style and uncapped intervals, unchanged samples and CIs.
- Fig. 8: larger sentence-case axis labels; original fixed-gap loess intervals.
- Fig. 9: no misleading pooled-average line or unsupported period −2 tick;
  common y-scale and separated counts.

The generic preset controls package presentation. The manuscript-specific
renderer combines saved summary tables so that it preserves heterogeneous
interval definitions and layouts across the three original analyses.
