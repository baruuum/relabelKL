#include <RcppArmadillo.h>
#include "utils.h"

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// add openmp support if available
#ifdef _OPENMP

  #include <omp.h>
  //[[Rcpp::plugins(openmp)]]
  
#endif

//' Relabel membership vector by minimizing KL-algorithm (log-scale)
//' 
//' Relabels the membership vectors of a mixed membership model or 
//' mixed membership stochastic blockmodel, using the KL-algorithm
//' proposed by Stephens (2000).
//' 
//' @param lphi cube of length \code{S}, each element of which is a 
//'     matrix of dimension N times K, where N is the number of individuals 
//'     and K is the number of extreme types (or classes).
//' @param renormalize if \code{true}, renormalizes the rows of each slice of \code{lphi}
//' @param maxit the number of maximum iterations to run, defaults to 100.
//' @param nthreads number of threads to use in parallel processing
//' @param verbose if TRUE, number of iterations and corresponding KL-distance 
//'     values are printed. 
//' @return The function returns a Rcpp::List of three elements. 
//'     1) A cube, \code{permuted}, of the same dimensions as \code{phi} but 
//'     with the labels permuted; 2) \code{perms} is a \code{S} times \code{K} 
//'     matrix containing the permutations necessary to produce \code{permuted} 
//'     from \code{phi} (i.e., the mapping from \code{phi} to \code{permuted});
//'     and 3) the number of iterations run.
//' @details OpenMP is enabled if available. Due to overhead, the inner-most
//'     loop is parallelized only if the number of latent classes/types
//'     is larger than 3.
//[[Rcpp::export]]
Rcpp::List relabel_kl_log(
    const arma::cube & lphi,
    bool renormalize = false,
    arma::uword maxit = 100L,
    arma::uword nthreads = 0L,
    bool verbose = true
) {
  
#ifdef _OPENMP
  
  if (nthreads == 0L)
     nthreads = omp_get_num_threads();
  
  omp_set_num_threads(nthreads);
  
#endif

    arma::uword S = lphi.n_slices;   // no. sims
    arma::uword K = lphi.n_cols;     // no. labels
    arma::uword N = lphi.n_rows;     // no. inds
    
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
        
    // initialize objects
    arma::mat lQ_hat(N, K), lQ_hat_new(N, K);
    double meankl(0.0), meankl_new(0.0);

#ifdef _OPENMP

    #pragma omp parallel for
    for (arma::uword i = 0; i < N; ++i) {

      for (arma::uword j = 0; j < K; ++j) {

        lQ_hat(i,j) = log_sum_exp(arma::vec(res(arma::span(i), arma::span(j), arma::span::all)));

      }

    }

    lQ_hat -= std::log(S);

    #pragma omp parallel for reduction(+: meankl)
    for (arma::uword s = 0; s < S; ++s) {

        meankl += kl_dist_log(res.slice(s), lQ_hat);

    }

#else
    
    for (arma::uword i = 0; i < N; ++i)
      for (arma::uword j = 0; j < K; ++j)
        lQ_hat(i,j) = log_sum_exp(arma::vec(res(arma::span(i), arma::span(j), arma::span::all)));
    
    lQ_hat -= std::log(S);
    
    for (arma::uword s = 0; s < S; ++s)
      meankl += kl_dist_log(res.slice(s), lQ_hat);
    
#endif
    
    meankl /= S;
    
    if (verbose) 
        Rcpp::Rcout << "\nStarting Fixed-Point Iterations ... \n" << std::endl;
    
    // START FIXED POINT INTERATIONS //
    
    for (arma::uword r = 0; r < maxit; ++r) {

        if (verbose) {
            Rcpp::Rcout <<
                std::setprecision(3) <<
                std::fixed <<
                "Iteration : " <<
                r << ",   " <<
                "Mean KL-divergence : " <<
                meankl <<
                std::endl;
        }
        
        for (arma::uword s = 0; s < S; ++s) {
            
            arma::mat lP_hat = res.slice(s);
            arma::colvec kl_q(n_perms);
        
            // calc permutation that minimizes KL-dist to p


#ifdef _OPENMP

            

            #pragma omp parallel for if (K > 3)
            for (arma::uword n = 0; n < n_perms; ++n) {

                kl_q(n) = kl_dist_log(permute_mat(lP_hat, perms.row(n).t(), 0L), lQ_hat);

            }

#else
            
            for (arma::uword n = 0; n < n_perms; ++n)
              kl_q(n) = kl_dist_log(permute_mat(lP_hat, perms.row(n).t(), 0L), lQ_hat);
              
#endif

            // index of min-KL permutation
            arma::uword choose_perm = kl_q.index_min();

            if (choose_perm != 0) {
                
                arma::uvec best_perm = perms.row(choose_perm).t();
                res.slice(s) = permute_mat(lP_hat, best_perm, 0L);
                perm_hist.row(s) = permute_urowvec(perm_hist.row(s), best_perm);
                
            }
            
        }
        
        // update Q_hat and meankl
        meankl_new = 0.0;
        
#ifdef _OPENMP

        #pragma omp parallel for
        for (arma::uword i = 0; i < N; ++i) {

            for (arma::uword j = 0; j < K; ++j) {

                lQ_hat_new(i,j) = log_sum_exp(arma::vec(res(arma::span(i), arma::span(j), arma::span::all)));

            }

        }

        lQ_hat_new -= std::log(S);

        #pragma omp parallel for reduction(+: meankl_new)
        for (arma::uword s = 0; s < S; ++s) {

            meankl_new += kl_dist_log(res.slice(s), lQ_hat);

        }
        
#else
        
        for (arma::uword i = 0; i < N; ++i)
          for (arma::uword j = 0; j < K; ++j)
            lQ_hat_new(i,j) = log_sum_exp(arma::vec(res(arma::span(i), arma::span(j), arma::span::all)));
        
        lQ_hat_new -= std::log(S);
        
        for (arma::uword s = 0; s < S; ++s)
          meankl_new += kl_dist_log(res.slice(s), lQ_hat_new);
        
#endif

        meankl_new /= S;
        
        // if there is no change, return
        if (meankl_new == meankl) {
            
            if (verbose) {
                Rcpp::Rcout << 
                    std::setprecision(3) <<
                    std::fixed <<
                    "Converged!   ( Final KL-divergence : " <<
                    meankl <<
                    " )\n" <<
                    std::endl;
            }
            
            return Rcpp::List::create(
                Named("permuted") = res,
                Named("perms") = perm_hist,
                Named("iterations") = r);
        }
        
        // otherwise, update
        lQ_hat = lQ_hat_new;
        meankl = meankl_new;
        
    }
    
    return Rcpp::List::create(
                Named("permuted") = res,
                Named("perms") = perm_hist,
                Named("iterations") = maxit);

}
