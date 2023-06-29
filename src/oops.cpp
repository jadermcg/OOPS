#include "utils.h"
#include "prob_utils.h"
#include "oops.h"

//'Runs Expectation Maximization OOPS and reestimates model parameters.
//'@name oops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model to be reestimated.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability for motif belongs to position w1, w2, w3, ..., wm.
//'@param beta Markov model thats represents the control senquences.
//[[Rcpp::export]]
Rcpp::List oops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w = 1) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  
  
  /**
   * Posteriori
   */
  arma::vec z(m);
  
  
  /**
   * Convergence control
   */  
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    new_alpha.fill(1e-100);
    for (const auto &seq : fasta) { // Foreach sequence
      double marginal_prob = 0.0;
      
      /**
       * E-STEP
       */
      
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        z[j] = w * probSeqGivenPos(seq, alpha, beta, j);
        marginal_prob += z[j];
      }
      
      /**
       * M-STEP
       */
      z = z / marginal_prob;
      soft_update(new_alpha, seq, z.t());
    }
    
    alpha = new_alpha / arma::accu(new_alpha.col(0));
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}

//'Runs Expectation Maximization OOPS and reestimates model parameters.
//'@name logoops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model to be reestimated.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability for motif belongs to position w1, w2, w3, ..., wm.
//'@param beta Markov model thats represents the control senquences.
//[[Rcpp::export]]
Rcpp::List logoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w = 1) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  arma::vec new_w(m);
  /**
   * Posteriori
   */
  arma::vec z(m);
  
  
  /**
   * Convergence control
   */  
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    new_alpha.fill(1e-100);
    for (const auto &seq : fasta) { // Foreach sequence
      double marginal_prob = 0.0;
      
      /**
       * E-STEP
       */
      
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        z[j] = std::log(w) + probSeqGivenPosLog(seq, alpha, beta, j);
        marginal_prob += std::exp(z[j]);
      }
      
      /**
       * M-STEP
       */
      z = arma::exp(z - std::log(marginal_prob));
      soft_update(new_alpha, seq, z.t());
    }
    
    alpha = new_alpha / arma::accu(new_alpha.col(0));
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}




