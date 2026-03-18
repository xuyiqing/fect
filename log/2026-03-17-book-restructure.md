# 2026-03-17 — Quarto Book Restructure

> Run: `REQ-book-cleanup` | Profile: r-package | Verdict: Pending — skeptic review follows scribe.

## What Changed

Major restructure of the fect Quarto user manual from 13 chapters to 12. Diagnostics and Model Selection chapters were dissolved into method-specific chapters (Ch2, Ch3, Ch4) to co-locate diagnostic tools with the estimators they test. Effect Heterogeneity was promoted from a section within Ch2 to a standalone chapter. The esplot.R cex multiplier was reverted after discovering it double-multiplied when called from plot.fect. Simulated datasets (simdata1, simdata2) were registered as package data with documentation. All 590 tests pass; all chapters render.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `vignettes/02-fect.Rmd` | modified | Restructured into 6 sections: Simulated data, FEct estimator, Visualization & diagnostics, Other estimands, Additional notes. Brief placebo material absorbed from old Ch5. |
| `vignettes/03-ife-mc.Rmd` | modified | Title changed to "Low-Rank Factor Methods". Full diagnostic suite absorbed from old Ch5. CV material absorbed from old Ch6. |
| `vignettes/04-cfe.Rmd` | modified | fect_mspe material absorbed from old Ch6. Existing diagnostics retained. |
| `vignettes/05-hte.Rmd` | created | New standalone chapter "Effect Heterogeneity" — promoted from Ch2 section. |
| `vignettes/06-plots.Rmd` | modified | Renumbered from 07. Availability table absorbed from old Ch6. |
| `vignettes/07-gsynth.Rmd` | modified | Renumbered from 08. |
| `vignettes/08-panel.Rmd` | modified | Renumbered from 09. |
| `vignettes/09-sens.Rmd` | modified | Renumbered from 10. |
| `vignettes/05-diagnostics.Rmd` | deleted | Dissolved into Ch2 (brief placebo) and Ch3 (full suite). |
| `vignettes/06-model-selection.Rmd` | deleted | Dissolved into Ch3 (CV), Ch4 (fect_mspe), Ch7/06 (availability table). |
| `vignettes/03-plots.Rmd` | deleted | Old superseded file (not in _quarto.yml). |
| `vignettes/04-gsynth.Rmd` | deleted | Old superseded file (not in _quarto.yml). |
| `vignettes/05-panel.Rmd` | deleted | Old superseded file (not in _quarto.yml). |
| `vignettes/_quarto.yml` | modified | Updated chapter list to 12 entries. |
| `vignettes/rscript/05-diagnostics.R` | deleted | Replaced by method-specific scripts. |
| `vignettes/rscript/06-model-selection.R` | deleted | Replaced by method-specific scripts. |
| `vignettes/rscript/05-hte.R` | created | R script for new Ch5. |
| `vignettes/rscript/06-plots.R` | modified | Renumbered from 07-plots.R. |
| `vignettes/rscript/07-gsynth.R` | modified | Renumbered from 08-gsynth.R. |
| `vignettes/rscript/08-panel.R` | modified | Renumbered from 09-panel.R. |
| `vignettes/rscript/09-sens.R` | modified | Renumbered from 10-sens.R. |
| `R/esplot.R` | modified | Reverted cex multiplier — was double-multiplying when called via plot.fect(). |
| `data/simdata1.rda` | created | Registered simulated dataset for Ch2. |
| `data/simdata2.rda` | created | Registered simulated dataset for Ch3. |
| `man/simdata1.Rd` | created | Documentation for simdata1. |
| `man/simdata2.Rd` | created | Documentation for simdata2. |
| `vignettes/05-diagnostics.Rmd` | modified | Fixed duplicate chunk name clash before dissolution. |
| `vignettes/09-panel.Rmd` | modified | Fixed polars/XQuartz rendering issue before renumbering. |
| `architecture.md` | modified | Updated Quarto book section with 12-chapter structure table. |

## Process Record

This section captures the full workflow history: what was proposed, what was tested, what problems arose, and how they were resolved.

### Proposal (from theorist)

**Implementation spec summary** (from `spec.md`):
- Dissolve Ch5 (Diagnostics) and Ch6 (Model Selection) into method-specific chapters to reduce redundancy and improve co-location of estimator + diagnostics content
- Promote Effect Heterogeneity from a section in Ch2 to standalone Ch5
- Renumber downstream chapters to fill gaps
- Register simdata1/simdata2 as proper package datasets with .Rd documentation
- Fix esplot.R cex multiplier (double-multiplying when plot.fect delegates to esplot)

**Test spec summary** (from `test-spec.md`):
- All 590 existing tests must continue to pass
- All 12 chapters must render without errors
- No orphaned chapter files (every .Rmd in vignettes/ must be in _quarto.yml or explicitly auxiliary)
- R scripts in rscript/ must align 1:1 with code chapters

### Implementation Notes (from builder)

- Ch2 restructured into 6 clean sections rather than appending diagnostics material at the end
- Diagnostics dissolution required splitting content: brief placebo intro to Ch2, full diagnostic suite (F-test, TOST, placebo, carryover) to Ch3 where IFE/MC methods need the most validation
- Model Selection dissolution: CV content naturally belongs with Ch3 (factor methods do CV for r), fect_mspe with Ch4 (CFE comparison tool), availability table with Ch6 (visualization chapter)
- esplot.R fix: reverted to passing raw cex values since plot.fect already applies the multiplier before calling esplot — applying it again produced double-sized text
- Fixed 05-diagnostics.Rmd chunk name clash (duplicate `placebo` chunk names) before dissolving content
- Fixed 09-panel.Rmd polars/XQuartz issue (rendering failed on macOS without X11 forwarding)

### Validation Results (from auditor)

| Test | Metric | Expected | Actual | Tolerance | Rel. Error | Verdict |
| --- | --- | --- | --- | --- | --- | --- |
| Full test suite | PASS count | 590 | 590 | 0 | 0% | PASS |
| Full test suite | FAIL count | 0 | 0 | 0 | — | PASS |
| Chapter rendering | Chapters built | 12 | 12 | 0 | 0% | PASS |
| Orphan check | Orphaned .Rmd files | 0 | 0 | 0 | — | PASS |
| Script alignment | Scripts matching chapters | 8 | 8 | 0 | 0% | PASS |

Summary: 590 tests executed, 590 passed, 0 failed. WARN 740 (all upstream: ggplot2 `size` deprecation, HonestDiDFEct NaN).

**Before/After Comparison Table**:

| Metric | Before (old) | After (new) | Change | Interpretation |
| --- | --- | --- | --- | --- |
| Chapter count | 13 | 12 | -1 | Improvement — consolidated redundant chapters |
| R scripts | 9 (02-10) | 8 (02-09) | -1 | Aligned — diagnostics/model-selection scripts removed |
| Orphaned .Rmd files | 5 | 0 | -5 | Improvement — all dead files deleted |
| Test pass count | 590 | 590 | 0 | Neutral — no regressions |

### Problems Encountered and Resolutions

| # | Problem | Signal | Routed To | Resolution |
| --- | --- | --- | --- | --- |
| 1 | 05-diagnostics.Rmd had duplicate `placebo` chunk names causing render failure | BLOCK | builder | Renamed chunks to unique names before dissolving content |
| 2 | 09-panel.Rmd failed to render due to polars/XQuartz dependency on macOS | BLOCK | builder | Added conditional rendering guard for X11-dependent plots |
| 3 | esplot.R cex values were double-multiplied (plot.fect applies multiplier, then esplot applied it again) | BLOCK | builder | Reverted esplot.R to pass raw cex values; multiplier applied only in plot.fect |

### Review Summary (from skeptic, if available)

Pending — skeptic review follows scribe.

- **Pipeline isolation**: pending
- **Convergence**: pending
- **Tolerance integrity**: pending
- **Verdict**: pending

## Design Decisions

1. **Dissolve diagnostics into method chapters rather than keeping a standalone chapter**: Diagnostics are method-specific — placebo tests make sense for FE (Ch2), the full F-test/TOST/carryover suite matters most for IFE/MC (Ch3) where factor misspecification is the primary risk, and CFE (Ch4) already had its own diagnostics. A standalone chapter forced readers to cross-reference constantly. Co-location reduces cognitive overhead.

2. **Dissolve model selection into method chapters rather than keeping standalone**: Cross-validation for r (number of factors) is inherently an IFE/MC concern (Ch3). fect_mspe is a comparison tool most useful when evaluating CFE vs alternatives (Ch4). The availability table is a visualization aid (Ch6). Grouping these under "Model Selection" created an artificial category that broke the estimator-centric narrative flow.

3. **Promote Effect Heterogeneity to standalone Ch5**: HTE content (group-level ATT, cohort effects, calendar effects) applies to all estimators equally and was growing large enough to warrant its own chapter. Keeping it buried in Ch2 made it hard to discover and gave the false impression it was FE-specific.

4. **Revert esplot.R cex multiplier rather than changing plot.fect**: plot.fect is the primary user-facing function and its multiplier pattern (`base_size * user_value`) is the established convention. esplot() is the lower-level function that plot.fect delegates to. The correct fix is for esplot() to accept raw values and let the caller (plot.fect) handle the multiplication. This avoids breaking direct esplot() users who pass absolute pt values.

5. **Register simdata1/simdata2 as package data**: The simulated datasets were previously generated inline in chapter code. Registering them as `.rda` with `.Rd` documentation makes chapters load faster (no simulation step), ensures reproducibility (fixed seed in data-raw), and lets users access them via `data(simdata1)` for their own experiments.

## Handoff Notes

- The 12-chapter structure is the final intended layout. Each method chapter (Ch2-Ch4) now contains its own diagnostics and model selection material — do not re-extract these into standalone chapters.
- R scripts in `rscript/` are numbered 02-09 matching code chapters 2-9. Ch1 (Get Started), Ch A (Cheatsheet), index, and references have no scripts.
- The esplot.R cex behavior is: NULL args use hardcoded defaults (16pt main, 12pt lab, 10pt axis, 10pt text); non-NULL args are treated as raw pt values. plot.fect applies its own multiplier before passing to esplot. Direct esplot() users pass absolute pt values.
- All 590 tests pass. The 740 warnings are all upstream (ggplot2 `size` aesthetic deprecation, HonestDiDFEct NaN in sensitivity analysis). These are not regressions.
- simdata1 uses a simple two-way FE DGP (for Ch2). simdata2 uses an interactive FE DGP with latent factors (for Ch3). Both are documented in `man/simdata1.Rd` and `man/simdata2.Rd`.
