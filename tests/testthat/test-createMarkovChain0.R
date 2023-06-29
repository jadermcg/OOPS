testthat::test_that("createMarkovChain0", {
  fasta = c("AACAGAC", "ACCACAG", "TTATATA", "ATCATAC")
  beta = createMarkovChain0(fasta)
  testthat::expect_true(1 - abs(sum(beta)) < 1e-6)
  testthat::expect_true(all(beta >= 0))
})