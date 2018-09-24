
#include <RcppArmadillo.h>
using namespace std;

/* a function to calculate an average over a NumericMatrix */
double average(Rcpp::NumericMatrix nm) {
  double sm = std::accumulate(nm.begin(), nm.end(), 0.0);
  double avg = sm / nm.size();
  return (avg);
}

/* a function to calculate an average over a NumericVector */
double average2(Rcpp::NumericVector nm) {
  double sm = std::accumulate(nm.begin(), nm.end(), 0.0);
  double avg = sm / nm.size();
  return (avg);
}

/* n_ = count of sampled values, values_ given values to sample from, probabilities_ = probabilities to use with sampling */
std::vector<double> cppSample(int n, std::vector<double> values, std::vector<double> probabilities) {

	//Rcpp::RNGScope scope;             // do we really need this?

	std::vector<double> cumsum(probabilities.size());
	std::partial_sum(probabilities.begin(),probabilities.end(),cumsum.begin(),std::plus<double>());

	//Rcpp::NumericVector rand_unif = Rcpp::runif(n, Rcpp::Named("min") = 0,Rcpp::Named("max") = maxint);
	std::vector<double> rand(n);
	for(int i = 0; i < n; i++ ) {
		rand[i] = ::Rf_runif(0,1);
	}

	std::vector<double> sampled_values(n);
	// tehdään eka luku erikseen ja loput loopilla
	for(int k = 0; k < rand.size(); k++) {
		if (rand[k] <= cumsum[0]) {
			sampled_values[k] = values[0];
		}
		for(int i = 1; i < cumsum.size(); i++) {
			if(rand[k] > cumsum[i-1] & rand[k] <= cumsum[i]) {
				sampled_values[k] = values[i];
			}
		}
	}
	return( sampled_values );
}