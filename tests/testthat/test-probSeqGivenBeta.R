testthat::test_that("probSeqGivenBeta", {
  beta = c(0.15, 0.25, 0.30, 0.3)
  kmers = expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"))
  kmers = apply(kmers, 1, paste0, collapse="")
  
  p = do.call(what = "c", args = lapply(kmers, function(x) probSeqGivenBeta(x,beta)))
  
  testthat::expect_true(1 - abs(sum(p)) < 1e-6)
  
})