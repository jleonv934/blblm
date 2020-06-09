//' @title faster linear regression with RcppArmadillo
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//' Below is a simple example of exporting a C++ function to R. You can
//' source this function into an R session using the Rcpp::sourceCpp
//' function or via the Source button on the editor toolbar

//' For more on using Rcpp click the Help button on the editor toolbar
//' @name fast_lm
//' @param x A single integer.
//' @export fast_lm

// [[Rcpp::export]]
List fast_lm(const arma::mat& X, const arma::colvec& y) {
  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals

  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);

  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

  return List::create(Named("coefficients") = coef,
                      Named("stderr")       = std_err,
                      Named("res")  = res,
                      Named("df.residual")  = n - k);
}
