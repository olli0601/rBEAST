context("beast2out.read.trees")

test_that("beast2out.read.trees: use", {
  trees_file <- "../../vignettes/example.trees"
  testit::assert(file.exists(trees_file))
  posterior <- beast2out.read.trees(trees_file)
  expect_equal(length(posterior), 10)
  expect_equal(class(posterior[[1]]), "phylo")
})

