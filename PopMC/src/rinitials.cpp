
#include <Rcpp.h>
/* double a, double b */
double rinitials(double a, double b) {
  double res = ::Rf_rgamma(a,b);
  return( res );
}