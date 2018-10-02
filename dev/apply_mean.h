// Originally by Kevin Ushey: https://gist.github.com/kevinushey/4561281
// the column-wise implementation is just as fast as colMeans,
// but the row-wise operation is not quite as fast as rowMeans.

// can I do better?

#include <Rcpp.h>
using namespace Rcpp;

template <class T>

inline double do_mean( T& x );
NumericVector row_means( NumericMatrix& X );
NumericVector col_means( NumericMatrix& X );
NumericVector Mean( NumericMatrix X, int dim );
