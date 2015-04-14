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
* 15-04-14: Pool phylogenetic clusters according to several options. The two new function is `beast.pool.clusters`. This can be useful if mixing is poor on a BEAST run comprising all taxa. In this case, sequences of the same cluster can be pooled into smaller taxon sets with this function.
