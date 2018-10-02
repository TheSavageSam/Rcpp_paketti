// Originally by Kevin Ushey: https://gist.github.com/kevinushey/4561281
// the column-wise implementation is just as fast as colMeans,
// but the row-wise operation is not quite as fast as rowMeans.

// can I do better?

#include <Rcpp.h>
using namespace Rcpp;

template <class T>
inline double do_mean( T& x ) {
  return mean(x);
}

NumericVector row_means( NumericMatrix& X ) {
  
  int nRows = X.nrow();
  NumericVector out = no_init(nRows);
  
  for( int i=0; i < nRows; i++ ) {
    NumericMatrix::Row tmp = X(i, _);
    out[i] = do_mean( tmp );
  }
  
  return out;
  
}

NumericVector col_means( NumericMatrix& X ) {
  
  int nCols = X.ncol();
  NumericVector out = no_init(nCols);
  
  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_mean( tmp );
  }
  
  return out;
  
}

NumericVector Mean( NumericMatrix X, int dim ) {
  
  if( dim == 1 ) {
    return row_means(X);
  } else if( dim == 2 ) {
    return col_means(X);
  } 
  
}

/*** R
library(rbenchmark)
x <- matrix( rnorm(1E6), ncol=1E3 )
benchmark( replications=5, order=NULL,
           apply(x, 1, mean),
           Mean(x, 1),
           rowMeans(x),
           apply(x, 2, mean),
           Mean(x, 2),
           colMeans(x)
)
*/