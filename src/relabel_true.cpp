#include <RcppArmadillo.h>
#include "utils.h"

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// add openmp support if available
#ifdef _OPENMP

#include <omp.h>
//[[Rcpp::plugins(openmp)]]

#endif


//' Relabel membership vector by minimizing KL-algorithm to true labels
//' 
//' Relabels the membership vectors of a mixed membership model or 
//' mixed membership stochastic blockmodel, by minimizing the KL-divergence
//' to a priori known true labels, where probabilities are entered on the
//' log-scale.
//' 
//' @param lphi cube of length \code{S}, each slice of which is a matrix of 
//'   dimension \code{N} times \code{K}, and where \code{N} is the number of 
//'   individuals and K is the number of extreme types (or classes). 
//'   Elements of \code{lphi} should be log-probabilities.
//' @param lphi_true a \code{N} times \code{K} matrix containing the "true"
//'   values of latent-class/mixed-membership probabilities on the log-scale
//' @param renormalize whether to renormalize the rows of \code{lphi}
//' @param nthreads number of threads to use if OPENMP is enabled
//' @param verbose if true, KL-divergence from true labels is calculated
//'   and printed before and after relabeling
//' @return Returns a Rcpp::List with two elements: 1) A arma::cube, 
//'   \code{permuted}, of the same dimensions as \code{phi} but 
//'   with the labels permuted; 2) \code{perms} is a \code{S} times \code{K} 
//'   matrix containing the permutations necessary to produce \code{permuted} 
//'   from \code{phi} (i.e., the mapping from \code{phi} to \code{permuted})
//' @details OpenMP is enabled if available. Due to overhead, the inner-most
//'     loop is parallelized only if the number of latent classes/types
//'     is larger than 3.
//[[Rcpp::export]]
Rcpp::List relabel_true_log(
        const arma::cube & lphi, 
        const arma::mat & lphi_true,
        bool renormalize = false,
        arma::uword nthreads = 0L,
        bool verbose = false
) {
    
    
#ifdef _OPENMP
    
    if (nthreads == 0L)
        nthreads = omp_get_num_threads();
    
    omp_set_num_threads(nthreads);
    
#endif
    
    arma::uword S = lphi.n_slices;   // no. sims
    arma::uword K = lphi.n_cols;     // no. labels
    
    arma::cube res(lphi);
    
    // renormalize if requested
    if (renormalize) {
        
        for (arma::uword s = 0; s < S; ++s) 
            res.slice(s).each_row( [](arma::rowvec& v){ v = v - log_sum_exp(v); } );
        
    }
    
    // generate matrix of possible permutations
    arma::umat perms = gen_permute(K);    
    
    // no. permutations
    arma::uword n_perms = perms.n_rows;     
    
    // matrix to store permutation history (initialize to 0,1,..K; each row)
    arma::umat perm_hist(S, K);
    perm_hist.each_row() = perms.row(0L);
    
    arma::uword s;
    
    // initial avg. KL-dist to post.mean
    
    if (verbose) {
        
        double meankl(0.0);
        
        
#ifdef _OPENMP
        
#pragma omp parallel for reduction(+: meankl)
        for (s = 0; s < S; ++s)
            meankl += kl_dist_log(res.slice(s), lphi_true);
        
#else
        
        for (s = 0; s < S; ++s)
            meankl += kl_dist_log(res.slice(s), lphi_true);
        
#endif
        
        meankl /= S;
        
        Rcpp::Rcout <<
            std::setprecision(3) <<
                std::fixed <<
                    "KL-divergence from true labels before relabeling : " <<
                        meankl <<
                            std::endl;
        
    }
    
    // START RELABELING //
    
    for (s = 0; s < S; ++s) {
        
        arma::mat lP_hat = res.slice(s);
        arma::colvec kl_q(n_perms);
        
        // calc permutation that minimizes KL-dist to p
        arma::uword n;
        
#ifdef _OPENMP
        
        #pragma omp parallel for if (K > 3)
        for (n = 0; n < n_perms; ++n) {
            
            kl_q(n) = kl_dist_log(permute_mat(lP_hat, perms.row(n).t(), 0L), lphi_true);
            
        }
        
#else
        
        for (n = 0; n < n_perms; ++n)
            kl_q(n) = kl_dist_log(permute_mat(lP_hat, perms.row(n).t(), 0L), lphi_true);
        
        
#endif
        
        // index of min-KL permutation
        arma::uword choose_perm = kl_q.index_min();
        
        if (choose_perm != 0) {
            
            arma::uvec best_perm = perms.row(choose_perm).t();
            res.slice(s) = permute_mat(lP_hat, best_perm, 0L);
            perm_hist.row(s) = permute_urowvec(perm_hist.row(s), best_perm);
            
        }
        
    }
    
    if (verbose) {
        
        // recalculate meankl
        double meankl_new(0.0);
        
#ifdef _OPENMP
        
#pragma omp parallel for reduction(+: meankl_new)
        for (s = 0; s < S; ++s) {
            meankl_new += kl_dist_log(res.slice(s), lphi_true);
        }
        
#else
        
        for (s = 0; s < S; ++s) 
            meankl_new += kl_dist_log(res.slice(s), lphi_true);
        
#endif
        
        meankl_new /= S;
        
        Rcpp::Rcout << 
            std::setprecision(3) <<
                std::fixed <<
                    "KL-divergence from true labels after relabeling : " <<
                        meankl_new <<
                            std::endl;
        
    }
    
    return Rcpp::List::create(
        Named("permuted") = res,
        Named("perms") = perm_hist);
    
}
