
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fect

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**R** package for implementing counterfactual estimators in panel
fixed-effect settings. It is suitable for panel/TSCS analysis with
binary treatments under (hypothetically) baseline randomization. It
allows a treatment to switch on and off and limited carryover effects.
It supports linear factor models—hence, a generalization of
[**gsynth**](https://yiqingxu.org/packages/gsynth/index.html)—and the
matrix completion method.

**Repo:** [GitHub](https://github.com/xuyiqing/fect) (0.4.1)

**Examples:** R code used in the
[tutorial](https://yiqingxu.org/packages/fect/articles/tutorial.html)
can be downloaded from [here](fect_examples.R).

**Reference:** Licheng Liu, Ye Wang, Yiqing Xu (2021). [A Practical
Guide to Counterfactual Estimators for Causal Inference with Time-Series
Cross-Sectional Data](https://papers.ssrn.com/abstract=3555463).
*American Journal of Political Science*, conditionally accepted.

## Installation

<!---
You can install **fect** directly from CRAN by typing the following command in the **R** console: 


```r
install.packages('fect', type = 'source')
```
--->

You can install the development version of **fect** from GitHub by
typing the following commands:

``` r
devtools::install_github('xuyiqing/fect')
```

**panelview** for panel data visualization is also highly recommended:

``` r
devtools::install_github('xuyiqing/panelView')
```

**fect** depends on the following packages, which will be installed
automatically when **fect** is being installed. You can also install
them manually.

``` r
## for processing C++ code
require(Rcpp) 
## for plotting
require(ggplot2)  
require(GGally) 
require(grid)
require(gridExtra)
## for parallel computing 
require(foreach)
require(future)  
require(doParallel) 
require(abind) 
```

### Notes on installation failures

1.  Mac users who have updated to MacOS BigSur or Monterey will likely
    encounter compilation problems. See
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

Please report bugs to **yiqingxu \[at\] stanford.edu** with your sample
code and data file. Much appreciated!
