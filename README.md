pedigreeR
=========
The ```pedigreeR``` package offers a suit of functions for pedigree analyses, including: sorting and editing pedigree data, computing inbreeding, additive relationships and functions of it (e.g., cholesky decomposistions, inverse of a numerator relationship matrix, etc.). The package was originally co-developed by Douglas M. Bates and Ana I. Vazquez and was incorporated in the [lme4](https://cran.r-project.org/web/packages/lme4/index.html) R-package. Recently, pedigreeR evolved to an stand-alone package. We have added new features for analysis of data from self-pollination.

Developers: Ana I. Vazquez & Douglas M. Bates.

Contributors: Paulino Perez-Rodriguez & Gustavo de los Campos.


###Installing pedigreeR from GitHub

```R
#R functions related to pedigrees
install.packages(pkg='devtools',repos='https://cran.r-project.org/')  #1# install devtools
library(devtools)                                                     #2# load the library
install_git('https://github.com/Rpedigree/pedigreeR/')                #3# install pedigreeR from GitHub
```
### Examples
  1. [Installing and loading library](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/install.md).
  2. [Completing and sorting a pedigree](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/example1.md).
  3. [Computing inbreeding, additive relationships and funcitons of it](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/example2.md).
  4. [Fiting pedigree mixed effects models using pedigreeR and lmer](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/example3.md).
  5. [Pedigree analyses for self-polinated species](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/example4.md).
  6. [Cross-validation analyses with pedigreeR and lmer](https://github.com/Rpedigree/pedigreeR/blob/master/inst/examples/example5.md).

