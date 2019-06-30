// #include <RcppArmadillo.h>
// //[[Rcpp::depends(RcppArmadillo)]]
// 
// using namespace Rcpp;
// 
// //' Permute vector
// //'
// //' Permutes a row or column vector in a specified order
// //' 
// //' @param x vector to permute
// //' @param order new order (indices have to start from 0)
// //' @return \code{x} permuted in order \code{order}
// template <typename T>
// inline T permute_vec(
//         const T & x,
//         const arma::uvec & order)
// {
// 
//     if (order.min() != 0)
//         Rcpp::stop("order vector has to start from 0 (permute_vec)!");
// 
// 
//     if (x.n_elem != order.n_elem)
//         Rcpp::stop("size mismatch (permute_vec)!");
// 
//     T res(x);
// 
//     for (arma::uword k = 0; k < x.n_elem; ++k)
//         res(k) = x(order(k));
//     
//     return res;
// 
// };
// 
// //' Permute rows or columns of matrix
// //'
// //' Permutes the rows or the columns of a matrix in a specified order
// //' 
// //' @param x matrix to permute
// //' @param order new order (indices have to start from 0)
// //' @param dim if 0, columns are permuted; otherwise, rows are permuted
// //' @return \code{x} permuted in order \code{order}
// template <typename T>
// inline T permute_mat(
//     const T & x,
//     const arma::uvec & order,
//     arma::uword dim)
// {
// 
//     if (order.min() != 0)
//         Rcpp::stop("order vector has to start from 0 (permute_mat)!");
// 
//     T res(x);
//     
//     if (dim == 0L) {
//         
//         if (x.n_cols != order.n_elem)
//             Rcpp::stop("dimension mismatch (permute_mat)!");
//         
//         for (arma::uword i = 0; i < x.n_cols; ++i)
//             res.col(i) = x.col(order(i));
//         
//         return res;
//         
//     } 
//         
//     if (x.n_rows != order.n_elem)
//        Rcpp::stop("dimension mismatch (permute_mat)!");
//     
//     for (arma::uword i = 0; i < x.n_rows; ++i) 
//         res.row(i) =  x.row(order(i));
//     
//     return res;
// 
// };
// 
// //[[Rcpp::export]]
// arma::vec tmppermvec(arma::vec x, arma::uvec o) {
//     return permute_vec(x,o);
// }
// 
// //[[Rcpp::export]]
// arma::mat tmppermmat(arma::mat x, arma::uvec o, arma::uword d) {
//     return permute_mat(x,o,d);
// }
