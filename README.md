
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fect

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[<img src="https://cranlogs.r-pkg.org/badges/grand-total/gsynth" alt="downloads: CRAN"/>](https://cran.r-project.org/web/packages/gsynth/index.html)
[<img src="https://cranlogs.r-pkg.org/badges/grand-total/fect" alt="downloads: CRAN"/>](https://cran.r-project.org/web/packages/fect/index.html)
<!-- badges: end -->

**R** package for implementing counterfactual estimators, also known as
imputation estimators, in panel fixed-effect settings. Suitable for
causal panel analysis with binary treatments under (hypothetically)
baseline randomization. It allows a treatment to switch on and off and
limited carryover effects. It supports two-way fixed effects, linear
factor models, and the matrix completion method.

Starting from v.2.0.0, all **gsynth** functionalities have been merged
into **fect**.

**Github Repo:** [GitHub](https://github.com/xuyiqing/fect) (2.0.0)

**User Manual:** R code used in the [User
Manual](https://yiqingxu.org/packages/fect/) can be downloaded from
[here](fect_examples.R).

**Main References:**

Xu, Yiqing (2017). [Generalized Synthetic Control Method: Causal
Inference with Interactive Fixed Effects
Models](https://www.cambridge.org/core/journals/political-analysis/article/generalized-synthetic-control-method-causal-inference-with-interactive-fixed-effects-models/B63A8BD7C239DD4141C67DA10CD0E4F3).
*Political Analysis* 25 (1): 57–76.

Licheng Liu, Ye Wang, Yiqing Xu (2024). [A Practical Guide to
Counterfactual Estimators for Causal Inference with Time-Series
Cross-Sectional
Data](https://yiqingxu.org/papers/english/2022_fect/LWX2022.pdf).
*American Journal of Political Science*, 68 (1): 160–76.

## Installation

To install **fect** from CRAN, run the code chunk below:

``` r
install.packages("fect")
```

We recommend users to install the most up-to-date version of **fect**
from Github using:

``` r
devtools::install_github("xuyiqing/fect")
```

After installation, check **fect** version to make sure the package is
up-to-date.

``` r
installed.packages()["fect", "Version"]
#> [1] "2.0.0"
```

**fect** depends on the following packages, which should be installed
automatically when **fect** is being installed. You can also install
them manually.

``` r
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}
packages <- c("abind", "doParallel", "doRNG", "fixest", "foreach", "future", 
              "GGally", "ggplot2", "grid", "gridExtra", "Mass", 
              "panelView", "Rcpp")
install_all(packages)
```

### Notes on installation failures

1.  Intel Mac users may encounter compilation problems. See
    [here](http://yiqingxu.org/public/BigSurError.pdf) for a potential
    solution.
2.  Windows users please consider upgrading R to 4.0.0 or higher and
    installing the [latest
    Rtools](https://cran.r-project.org/bin/windows/Rtools/) to avoid
    C++17 complier errors when installing fastplm.
3.  For Rcpp, RcppArmadillo and MacOS “-lgfortran” and “-lquadmath”
    error, click
    [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)
    for details.
4.  Installation failure related to OpenMP on MacOS, click
    [here](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/)
    for a solution.
5.  To fix these issues, try installing gfortran from
    [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS%20clang4%20R%20Binaries%20from%20https://github.com/coatless/r-macos-clang).

## Report bugs

Please report any bugs to me (yiqingxu \[at\] stanford.edu) or submit an
issue on [GitHub](https://github.com/xuyiqing/fect/issues). Please
include your minimally replicable code & data file and a **panelView**
treatment status plot. Your feedback is highly valued!
