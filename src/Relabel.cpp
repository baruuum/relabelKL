#include <RcppArmadillo.h>
#include "utils.h"
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' Relabel membership vector by minimizing KL-distance to (unknown) optimal labels
//' 
//' Relabels the membership vectors of a mixed membership model or 
//' mixed membership stochastic blockmodel, using the KL-algorithm
//' proposed by Stephens (2000).
//' 
//' @param phi cube of length \code{S}, each element of which is a 
//'        matrix of dimension N times K, where N is the number of individuals 
//'        and K is the number of extreme types (or classes).
//' @param max_iter the number of maximum iterations to run, defaults to 100.
//' @param verbose if TRUE, number of iterations and corresponding KL-distance 
//'        values are printed. 
//' @return A Rcpp::List of three elements. 1) A cube, \code{relabeled}, of 
//'         the same dimensions as \code{phi} but with the labels permuted. 
//'         2) \code{perms} is a \code{S} times \code{K} matrix containing 
//'         the permutations necessary to produce \code{relabeled} from 
//'         \code{phi} (i.e., the mapping from \code{phi} to \code{relabeled}).
//'         3) the number of iterations run.

//[[Rcpp::export]]
Rcpp::List relabel_kl(const arma::cube & phi, 
                      arma::uword maxit = 100,
                      bool verbose = true) {

    arma::uword S = phi.n_slices;   // no. sims
    arma::uword K = phi.n_cols;     // no. labels
    
    arma::cube res(phi);

    // generate matrix of possible permutations
    arma::umat perms = gen_permute(K);    

        // no. permutations
    arma::uword n_perms = perms.n_rows;     
    
    // matrix to store permutation history (initialize to 0,1,..K; each row)
    arma::umat perm_hist(S, K);
    perm_hist.each_row() = perms.row(0L);
        
    // post.mean phi 
    arma::mat Q_hat = mean(phi, 2L);
    
    // initial avg. KL-dist to post.mean
    arma::uword s;
    double meankl(0.0);
    for (s = 0; s < S; ++s)
        meankl += kl_dist(res.slice(s), Q_hat);
    meankl = meankl / S;
    
    if (verbose) 
        Rcpp::Rcout << "Starting Fixed-Point Iterations ... \n" << std::endl;
    
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
                
        for (s = 0; s < S; ++s) {
            
            arma::mat P_hat = res.slice(s);
            arma::colvec kl_q(n_perms);
            
            // calc permutation that minimizes KL-dist to p
            for (arma::uword n = 0; n < n_perms; ++n) {
                
                arma::mat q_perm = permute_mat(P_hat, perms.row(n).t(), 0L);
                kl_q(n) = kl_dist(q_perm, Q_hat);
                
            }
            // index of min-KL permutation
            arma::uword choose_perm = kl_q.index_min();

            if (choose_perm != 0) {
                
                arma::uvec best_perm = perms.row(choose_perm).t();
                res.slice(s) = permute_mat(P_hat, best_perm, 0L);
                perm_hist.row(s) = permute_urowvec(perm_hist.row(s), best_perm);
                
            }
            
        }
        
        // update Q_hat and meankl
        arma::mat Q_hat_new = arma::mean(res, 2);
        
        double meankl_new(0.0);
        for (s = 0; s < S; ++s) 
            meankl_new += kl_dist(res.slice(s), Q_hat_new);
        meankl_new = meankl_new / S;

        // if there is no change, return
        if (meankl_new == meankl) {
            
            if (verbose) {
                Rcpp::Rcout << 
                    std::setprecision(3) <<
                    std::fixed <<
                    "Converged!   ( Final KL-divergence : " <<
                    meankl <<
                    " )" <<
                    std::endl;
            }
            
            return Rcpp::List::create(
                Named("relabeled") = res,
                Named("perms") = perm_hist,
                Named("iterations") = r);
        }
        
        // otherwise, update
        Q_hat = Q_hat_new;
        meankl = meankl_new;
        
    }
    
    return Rcpp::List::create(
                Named("relabeled") = res,
                Named("perms") = perm_hist,
                Named("iterations") = maxit);

}
