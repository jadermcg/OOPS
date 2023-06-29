// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// oops
Rcpp::List oops(const std::vector<std::string>& fasta, arma::mat alpha, const arma::mat& beta, double cutoff, int niter, double w);
RcppExport SEXP _OOPS_oops(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cutoffSEXP, SEXP niterSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(oops(fasta, alpha, beta, cutoff, niter, w));
    return rcpp_result_gen;
END_RCPP
}
// logoops
Rcpp::List logoops(const std::vector<std::string>& fasta, arma::mat alpha, const arma::mat& beta, double cutoff, int niter, double w);
RcppExport SEXP _OOPS_logoops(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cutoffSEXP, SEXP niterSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(logoops(fasta, alpha, beta, cutoff, niter, w));
    return rcpp_result_gen;
END_RCPP
}
// Q
double Q(const std::vector<std::string>& fasta, const arma::mat& alpha, const arma::mat& beta);
RcppExport SEXP _OOPS_Q(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(Q(fasta, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// LL
double LL(const std::vector<std::string>& fasta, const arma::mat& alpha, const arma::mat& beta);
RcppExport SEXP _OOPS_LL(SEXP fastaSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(LL(fasta, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// createMarkovChain
arma::mat createMarkovChain(const std::vector<std::string>& fasta, const int tau);
RcppExport SEXP _OOPS_createMarkovChain(SEXP fastaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(createMarkovChain(fasta, tau));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenAlpha
double probSeqGivenAlpha(const std::string& kmer, const arma::mat& alpha);
RcppExport SEXP _OOPS_probSeqGivenAlpha(SEXP kmerSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenAlpha(kmer, alpha));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenAlphaLog
double probSeqGivenAlphaLog(const std::string& kmer, const arma::mat& alpha);
RcppExport SEXP _OOPS_probSeqGivenAlphaLog(SEXP kmerSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenAlphaLog(kmer, alpha));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenBeta
double probSeqGivenBeta(const std::string& seq, const arma::mat& beta);
RcppExport SEXP _OOPS_probSeqGivenBeta(SEXP seqSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenBeta(seq, beta));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenBetaLog
double probSeqGivenBetaLog(const std::string& seq, const arma::mat& beta);
RcppExport SEXP _OOPS_probSeqGivenBetaLog(SEXP seqSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenBetaLog(seq, beta));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenPos
double probSeqGivenPos(const std::string& seq, const arma::mat& alpha, const arma::mat& beta, const int pos);
RcppExport SEXP _OOPS_probSeqGivenPos(SEXP seqSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenPos(seq, alpha, beta, pos));
    return rcpp_result_gen;
END_RCPP
}
// probSeqGivenPosLog
double probSeqGivenPosLog(const std::string& seq, const arma::mat& alpha, const arma::mat& beta, const int pos);
RcppExport SEXP _OOPS_probSeqGivenPosLog(SEXP seqSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(probSeqGivenPosLog(seq, alpha, beta, pos));
    return rcpp_result_gen;
END_RCPP
}
// computeIC
double computeIC(const arma::mat alpha, const arma::mat beta);
RcppExport SEXP _OOPS_computeIC(SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeIC(alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// computeICU
double computeICU(const arma::mat alpha);
RcppExport SEXP _OOPS_computeICU(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeICU(alpha));
    return rcpp_result_gen;
END_RCPP
}
// soft_update
void soft_update(arma::mat& alpha, const std::string& seq, const arma::rowvec& posteriori);
RcppExport SEXP _OOPS_soft_update(SEXP alphaSEXP, SEXP seqSEXP, SEXP posterioriSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type posteriori(posterioriSEXP);
    soft_update(alpha, seq, posteriori);
    return R_NilValue;
END_RCPP
}
// hard_update
void hard_update(arma::mat& alpha, const std::string& seq, const arma::rowvec& posteriori);
RcppExport SEXP _OOPS_hard_update(SEXP alphaSEXP, SEXP seqSEXP, SEXP posterioriSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type posteriori(posterioriSEXP);
    hard_update(alpha, seq, posteriori);
    return R_NilValue;
END_RCPP
}
// kmers2alpha
arma::mat kmers2alpha(const std::vector<std::string>& kmers);
RcppExport SEXP _OOPS_kmers2alpha(SEXP kmersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type kmers(kmersSEXP);
    rcpp_result_gen = Rcpp::wrap(kmers2alpha(kmers));
    return rcpp_result_gen;
END_RCPP
}
// alpha2kmers
std::vector<std::string> alpha2kmers(const arma::mat& alpha, const std::vector<std::string>& fasta);
RcppExport SEXP _OOPS_alpha2kmers(SEXP alphaSEXP, SEXP fastaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha2kmers(alpha, fasta));
    return rcpp_result_gen;
END_RCPP
}
// fasta2alpha
arma::mat fasta2alpha(const std::vector<std::string>& fasta, const int k);
RcppExport SEXP _OOPS_fasta2alpha(SEXP fastaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(fasta2alpha(fasta, k));
    return rcpp_result_gen;
END_RCPP
}
// corr
std::string corr(const std::string& a, const std::string b);
RcppExport SEXP _OOPS_corr(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(corr(a, b));
    return rcpp_result_gen;
END_RCPP
}
// corr_freq
double corr_freq(const std::string& correlation_str);
RcppExport SEXP _OOPS_corr_freq(SEXP correlation_strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type correlation_str(correlation_strSEXP);
    rcpp_result_gen = Rcpp::wrap(corr_freq(correlation_str));
    return rcpp_result_gen;
END_RCPP
}
// fast_corr_freq
double fast_corr_freq(const std::string& a, const std::string b);
RcppExport SEXP _OOPS_fast_corr_freq(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_corr_freq(a, b));
    return rcpp_result_gen;
END_RCPP
}
// computeDKL
double computeDKL(const arma::mat& alpha, const arma::mat& beta, const std::string& kmer);
RcppExport SEXP _OOPS_computeDKL(SEXP alphaSEXP, SEXP betaSEXP, SEXP kmerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDKL(alpha, beta, kmer));
    return rcpp_result_gen;
END_RCPP
}
// computeDKLU
double computeDKLU(const arma::mat& alpha, const std::string& kmer);
RcppExport SEXP _OOPS_computeDKLU(SEXP alphaSEXP, SEXP kmerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kmer(kmerSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDKLU(alpha, kmer));
    return rcpp_result_gen;
END_RCPP
}
// computeSCORE
int computeSCORE(const arma::mat& alpha);
RcppExport SEXP _OOPS_computeSCORE(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(computeSCORE(alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OOPS_oops", (DL_FUNC) &_OOPS_oops, 6},
    {"_OOPS_logoops", (DL_FUNC) &_OOPS_logoops, 6},
    {"_OOPS_Q", (DL_FUNC) &_OOPS_Q, 3},
    {"_OOPS_LL", (DL_FUNC) &_OOPS_LL, 3},
    {"_OOPS_createMarkovChain", (DL_FUNC) &_OOPS_createMarkovChain, 2},
    {"_OOPS_probSeqGivenAlpha", (DL_FUNC) &_OOPS_probSeqGivenAlpha, 2},
    {"_OOPS_probSeqGivenAlphaLog", (DL_FUNC) &_OOPS_probSeqGivenAlphaLog, 2},
    {"_OOPS_probSeqGivenBeta", (DL_FUNC) &_OOPS_probSeqGivenBeta, 2},
    {"_OOPS_probSeqGivenBetaLog", (DL_FUNC) &_OOPS_probSeqGivenBetaLog, 2},
    {"_OOPS_probSeqGivenPos", (DL_FUNC) &_OOPS_probSeqGivenPos, 4},
    {"_OOPS_probSeqGivenPosLog", (DL_FUNC) &_OOPS_probSeqGivenPosLog, 4},
    {"_OOPS_computeIC", (DL_FUNC) &_OOPS_computeIC, 2},
    {"_OOPS_computeICU", (DL_FUNC) &_OOPS_computeICU, 1},
    {"_OOPS_soft_update", (DL_FUNC) &_OOPS_soft_update, 3},
    {"_OOPS_hard_update", (DL_FUNC) &_OOPS_hard_update, 3},
    {"_OOPS_kmers2alpha", (DL_FUNC) &_OOPS_kmers2alpha, 1},
    {"_OOPS_alpha2kmers", (DL_FUNC) &_OOPS_alpha2kmers, 2},
    {"_OOPS_fasta2alpha", (DL_FUNC) &_OOPS_fasta2alpha, 2},
    {"_OOPS_corr", (DL_FUNC) &_OOPS_corr, 2},
    {"_OOPS_corr_freq", (DL_FUNC) &_OOPS_corr_freq, 1},
    {"_OOPS_fast_corr_freq", (DL_FUNC) &_OOPS_fast_corr_freq, 2},
    {"_OOPS_computeDKL", (DL_FUNC) &_OOPS_computeDKL, 3},
    {"_OOPS_computeDKLU", (DL_FUNC) &_OOPS_computeDKLU, 2},
    {"_OOPS_computeSCORE", (DL_FUNC) &_OOPS_computeSCORE, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_OOPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
