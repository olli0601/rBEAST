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

* 15-04-14: Generate multi locus XML files for the HKY and GTR model, with fixed dated tree. The two new functions are `beastxml.multilocus.hky.fixed.tree` and `beastxml.multilocus.codon.gtr`.
* 15-04-14: Pool sequences in phylogenetic clusters for a BEAST analysis according to several options. The new function is `beast.pool.clusters`. This can be useful if mixing is poor on a BEAST run comprising all taxa. In this case, sequences of the same cluster can be pooled into smaller taxon sets with this function.
* 15-04-16: Generate relatively simple XML files from a template. The new function is `beastxml.from.template`. This can be useful to generate multiple XML files, each for a pool of sequences, when a (simple) template file is already at hand. This function can add monophyly constraints and logged tmrca statistics for each cluster. In this case, it may be necessary to specify a starting tree with the function `beast.get.startingtree`.

