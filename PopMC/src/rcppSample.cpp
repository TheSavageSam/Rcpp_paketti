
#include <Rcpp.h>
/* n_ = count of sampled values, values_ given values to sample from, probabilities_ = probabilities to use with sampling */
NumericVector cppSample(int n_, NumericVector values_, NumericVector probabilities_) {
  std::vector<double> sampled_values = cppSample(n_, as<std::vector<double> >(values_), as<std::vector<double> >(probabilities_));
	return( wrap(sample_values) );
}