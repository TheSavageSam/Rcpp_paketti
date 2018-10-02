#include <Rcpp.h>
#include "k.h"

Rcpp::IntegerVector rk(Rcpp::IntegerVector smokes_hav, Rcpp::NumericVector age, Rcpp::IntegerVector ymis_index, Rcpp::NumericVector theta) {
  std::vector<int> smokes_hav_ = k(Rcpp::as<std::vector<int> >(smokes_hav), Rcpp::as<std::vector<double> >(age), Rcpp::as<std::vector<int> >(ymis_index), Rcpp::as<std::vector<double> >(theta));
  return( Rcpp::wrap(smokes_hav_) );
}
