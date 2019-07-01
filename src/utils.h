#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Logarithm of sum of exponentials 
//'
//' @param x any object that allows iterators
//' @return returns the log of the sum of exponentiated elements of \code{x}
template <typename T>
inline double log_sum_exp(
        const T & x
) {

  double max_exp = x(0);
  double esum(0.0);
  
  typename T::const_iterator i; 

  for (i = x.begin() + 1 ; i != x.end() ; ++i)
    if (*i > max_exp)
      max_exp = *i;

  for (i = x.begin(); i != x.end() ; ++i)
    esum += std::exp(*i - max_exp);

  return std::log(esum) + max_exp;

};


//' Logarithm of the sum of the exponential of two real numbers
//'
//' @param x,y two real numbers
//' @return returns the log of the sum of the exponential of \code{x} and \code{y}
inline double log_add_exp(const double x, const double y) {

    double d = x - y;
    if (d > 0.0)
        return x + std::log1p(std::exp(-d));
    if (d <= 0.0)
        return y + std::log1p(std::exp(d));

}

//' Logarithm of accumulative sum of exponentials
//'
//' @param x any object that allows iterators
//' @return returns the log of the accumulative sum of exponentiated elements of \code{x}
template <typename T>
inline T log_accu_exp (const T& x) {

    T y(x);
    typename T::iterator i = y.begin() + 1;

    for (; i < y.end(); i++) {
        *i = log_add_exp(*(i - 1), *i);
    }

    return y;
};

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
        Rcpp::stop("size mismatch (kl_dist)!");

    if (arma::any(vectorise(p) <= 0.0) || arma::any(vectorise(q) <= 0.0))
        Rcpp::stop("negative probabilities input vectors (kl_dist)");

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
        Rcpp::stop("N has to be positive (gen_permute)");

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

        // note: std::next_permutation will generate the next
        //       permutation in lexicographical order and return
        //       "false" if the generated sequence is the first seq.

        ++it;

        pmat.row(it) = seq;

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
template <typename T>
inline T permute_vec(
        const T & x,
        const arma::uvec & order)
{

    if (order.min() != 0)
        Rcpp::stop("order vector has to start from 0 (permute_vec)!");


    if (x.n_elem != order.n_elem)
        Rcpp::stop("size mismatch (permute_vec)!");

    T res(x);

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
        Rcpp::stop("order vector has to start from 0 (permute_mat)!");

    T res(x);
    
    if (dim == 0L) {
        
        if (x.n_cols != order.n_elem)
            Rcpp::stop("dimension mismatch (permute_mat)!");
        
        for (arma::uword i = 0; i < x.n_cols; ++i)
            res.col(i) = x.col(order(i));
        
        return res;
        
    } 
        
    if (x.n_rows != order.n_elem)
       Rcpp::stop("dimension mismatch (permute_mat)!");
    
    for (arma::uword i = 0; i < x.n_rows; ++i) 
        res.row(i) =  x.row(order(i));
    
    return res;

};

#endif
