// [[Rcpp::depends(RcppArmadillo)]]
#pragma once
#include <RcppArmadillo.h>
#include <string.h>
#include <vector>
#include <omp.h>

Rcpp::List oops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w);
Rcpp::List logoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w);
