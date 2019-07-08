#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' KL Divergence between two distributions
//'
//' Calculates the KL-divergence of the target to the true distribution.
//' Notice that the KL-divergence is defined only if the target distribution is
//' absolutely continuous with respect to the true distribution. Hence, if
//' \code{q(i)} = 0 but \code{p(i)} != 0, for some i, the function will return
//' R_Inf. On the other hand, for \code{p(i)} = 0, \code{p(i)} * log(\code{p(i)}/\code{q(i)}
//' is defined to be equal zero.
//'
//' @param p the "true" distribution
//' @param q the "target" distribution
//' @return Returns the KL-divergence of \code{q} from \code{p}
template <typename T>
inline double kl_dist(
        const T & p,
        const T & q)
{
    if (p.n_elem != q.n_elem)
        Rcpp::stop("size mismatch (kl_dist)");

    if (arma::any(arma::vectorise(p) <= 0.0) || arma::any(arma::vectorise(q) <= 0.0))
        Rcpp::stop("negative probabilities in input vectors (kl_dist)");

    double kl(0.0);
    for (arma::uword i = 0; i < p.n_elem; ++i) {
        if (p(i) != 0)
            kl += p(i) * (std::log(p(i)) - std::log(q(i)));
    }

    return kl;

};

//' Generate permutations
//'
//' Generates a matrix all possible permutations of N objects
//'
//' @param N integer of number of elements
//' @return matrix of all possible permutations of length \code{N}.
//'         Returns a matrix of dimension n! times \code{N}
inline arma::umat gen_permute(
        arma::uword N)
{

    // check argument
    if (N < 0)
        Rcpp::stop("N has to be a positive integer (gen_permute)");

    // no. of possible permutations
    arma::uword K = R::gammafn(N + 1);

    // original sequence in ascending order
    arma::urowvec seq = arma::linspace<arma::urowvec>(0, N - 1, N);
    // matrix to store results
    arma::umat pmat(K, N, arma::fill::zeros);

    // add org. seq. into the first row
    pmat.row(0) = seq;
    // iterator
    arma::uword it(0L);

    while(std::next_permutation(seq.begin(), seq.end())) {

        pmat.row(++it) = seq;

    }

    return pmat;

};

//' Permute vector
//'
//' Permutes a row or column vector in a specified order
//' 
//' @param x vector to permute
//' @param order new order (indices have to start from 0)
//' @return \code{x} permuted in order \code{order}
inline arma::urowvec permute_urowvec(
        const arma::urowvec & x,
        const arma::uvec & order)
{

    if (order.min() != 0)
        Rcpp::stop("order vector has to start from 0 (permute_vec)");


    if (x.n_elem != order.n_elem)
        Rcpp::stop("size mismatch (permute_vec)");

    arma::urowvec res(arma::size(x));

    for (arma::uword k = 0; k < x.n_elem; ++k)
        res(k) = x(order(k));
    
    return res;

};

//' Permute rows or columns of matrix
//'
//' Permutes the rows or the columns of a matrix in a specified order
//' 
//' @param x matrix to permute
//' @param order new order (indices have to start from 0)
//' @param dim if 0, columns are permuted; otherwise, rows are permuted
//' @return \code{x} permuted in order \code{order}
template <typename T>
inline T permute_mat(
    const T & x,
    const arma::uvec & order,
    arma::uword dim)
{

    if (order.min() != 0)
        Rcpp::stop("order vector has to start from 0 (permute_mat)");

    T res(x);
    
    if (dim == 0L) {
        
        if (x.n_cols != order.n_elem)
            Rcpp::stop("dimension mismatch (permute_mat)");
        
        for (arma::uword i = 0; i < x.n_cols; ++i)
            res.col(i) = x.col(order(i));
        
        return res;
        
    } 
        
    if (x.n_rows != order.n_elem)
       Rcpp::stop("dimension mismatch (permute_mat)");
    
    for (arma::uword i = 0; i < x.n_rows; ++i) 
        res.row(i) =  x.row(order(i));
    
    return res;

};

#endif
