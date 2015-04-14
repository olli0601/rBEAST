# rBEAST
Automated generation of BEAST XML files, with a focus on dating phylogenies of HIV sequences

# Installation

The easiest way to install `rBEAST` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("olli0601/rBEAST")
```
# Content:

* 15-04-14: Generate multi locus XML files for the HKY and GTR model, with fixed dated tree. The two new functions are `beastscript.multilocus.hky` and `beastscript.multilocus.codon.gtr`.
