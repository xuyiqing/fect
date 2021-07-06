# fect

## Fixed Effects Counterfactual Estimators
---

**Authors:** Licheng Liu [<liulch.16@sem.tsinghua.edu.cn>]; Ye Wang[<yw1576@nyu.edu>], Yiqing Xu [<yiqingxu@stanford.edu>]

**Maintainer:** Licheng Liu

**How to Uses:** [Examples](https://yiqingxu.org/packages/fect/fect.html)

**Reference:**  Licheng Liu, Ye Wang, Yiqing Xu (2021). "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data." [Working Paper](https://arxiv.org/abs/2107.00856), Stanford University.

**Note:**

Rcpp, RcppArmadillo and MacOS "-lgfortran" and "-lquadmath" error, see: http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

Installation failture related to OpenMP on MacOS, see:
http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/

To fix these problems, consider installing: 
gfortran 6.1 from https://gcc.gnu.org/wiki/GFortranBinaries#MacOS
clang4 R Binaries from https://github.com/coatless/r-macos-clang
