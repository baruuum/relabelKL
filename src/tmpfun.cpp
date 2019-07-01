// #include <RcppArmadillo.h>
// //[[Rcpp::depends(RcppArmadillo)]]
// 
// using namespace Rcpp;
// 
// //[[Rcpp::export]]
// arma::umat gen_permute(
//         arma::uword N)
// {
// 
//     // check argument
//     if (N < 0)
//         Rcpp::stop("N has to be positive (gen_permute)");
// 
//     // no. of possible permutations
//     arma::uword K = R::gammafn(N + 1);
// 
//     // original sequence in ascending order
//     arma::urowvec seq = arma::linspace<arma::urowvec>(0, N - 1, N);
//     // matrix to store results
//     arma::umat pmat(K, N, arma::fill::zeros);
// 
//     // add org. seq. into the first row
//     pmat.row(0) = seq;
//     // iterator
//     arma::uword it(0L);
// 
//     while(std::next_permutation(seq.begin(), seq.end())) {
// 
//         // note: std::next_permutation will generate the next
//         //       permutation in lexicographical order and return
//         //       "false" if the generated sequence is the first seq.
// 
//         ++it;
// 
//         pmat.row(it) = seq;
// 
//     }
// 
//     return pmat;
// 
// }
// 
// 
// /*** R
// 
// gen_permute(3)
// */