
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fect

<!-- badges: start -->

[![Lifecycle:
stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[<img src="https://www.r-pkg.org/badges/version/fect" alt="CRAN status"/>](https://CRAN.R-project.org/package=fect)
[<img src="https://cranlogs.r-pkg.org/badges/grand-total/fect" alt="CRAN downloads"/>](https://cran.r-project.org/web/packages/fect/index.html)
<!-- badges: end -->

**R** package for implementing counterfactual estimators, also known as
imputation estimators, in panel fixed-effect settings. Suitable for
causal panel analysis with binary treatments under (hypothetically)
baseline randomization. It allows a treatment to switch on and off and
limited carryover effects. It supports two-way fixed effects, linear
factor models, and the matrix completion method.

Starting from v.2.0.0, all **gsynth** functionalities have been merged
into **fect**.

**Source Code:** [GitHub](https://github.com/xuyiqing/fect)

**User Manual:** [Quarto Book](https://yiqingxu.org/packages/fect/)

**Main References:**

Xu, Yiqing (2017). [Generalized Synthetic Control Method: Causal
Inference with Interactive Fixed Effects
Models](https://www.cambridge.org/core/journals/political-analysis/article/generalized-synthetic-control-method-causal-inference-with-interactive-fixed-effects-models/B63A8BD7C239DD4141C67DA10CD0E4F3).
*Political Analysis* 25 (1): 57–76.

Liu, Licheng, Ye Wang, Yiqing Xu (2024). [A Practical Guide to
Counterfactual Estimators for Causal Inference with Time-Series
Cross-Sectional
Data](https://yiqingxu.org/papers/english/2022_fect/LWX2022.pdf).
*American Journal of Political Science*, 68 (1): 160–76.

Chiu, Albert, Xingchen Lan, Ziyi Liu, and Yiqing Xu. (2025). [Causal
Panel Analysis Under Parallel Trends: Lessons from a Large Reanalysis
Study](https://www.cambridge.org/core/journals/american-political-science-review/article/causal-panel-analysis-under-parallel-trends-lessons-from-a-large-reanalysis-study/219275E0CE901F099F2CFFBA07079243).
*American Political Science Review*, First View.

**Report bugs:** Please report any bugs by submitting an issue on
[GitHub](https://github.com/xuyiqing/fect/issues) or emailing me
(yiqingxu \[at\] stanford.edu). We’d really appreciate it if you can
include your minimally replicable code & data file and a **panelView**
treatment status plot. Your feedback is highly valued!
