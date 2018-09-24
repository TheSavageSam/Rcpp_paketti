#include <Rcpp.h>
#include "k.h"

IntegerVector rk(IntegerVector smokes_hav, NumericVector age, IntegerVector ymis_index, NumericVector theta) {
  std::vector<int> smokes_hav = k(as<std::vector<int> >(smokes_hav), as<std::vector<double> >(age), as<std::vector<int> >(ymis_index), as<std::vector<double> >(theta));
  return( wrap(smokes_hav) );
}
