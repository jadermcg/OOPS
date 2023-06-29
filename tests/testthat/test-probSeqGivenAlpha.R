testthat::test_that("probSeqGivenAlpha", {
  alpha = t(extraDistr::rdirichlet(n = 5, rep(1, 4)))
  kmers = expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T"))
  kmers = apply(kmers, 1, paste0, collapse="")
  p = do.call(what = "c", args = lapply(kmers, function(x) probSeqGivenAlpha(x, alpha)))
  testthat::expect_true(1 - abs(sum(p)) < 1e-6)
})