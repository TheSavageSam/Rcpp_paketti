#include <Rcpp.h>
#include "cppSample.h"

/* n_ = count of sampled values, values_ given values to sample from, probabilities_ = probabilities to use with sampling */
Rcpp::NumericVector rcppSample(int n_, Rcpp::NumericVector values_, Rcpp::NumericVector probabilities_) {
  std::vector<double> sampled_values = cppSample(n_, as<std::vector<double> >(values_), as<std::vector<double> >(probabilities_));
	return( Rcpp::wrap(sample_values) );
}