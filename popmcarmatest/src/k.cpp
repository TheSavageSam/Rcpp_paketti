
#include <Rcpp.h>
double expit(double x) { return(exp(x)/(1+exp(x))); }

std::vector<int> k(std::vector<int> smokes_hav, std::vector<double> age, std::vector<int> ymis_index, std::vector<double> theta) {
	// df: data
  // theta: theta-parametrivektori (2 alkiota)
  int i; //int s=0;
  //std::vector<double> p_tmp;
  double p_tmp;
  for(int i_ymis = 0; i_ymis < ymis_index.size(); i_ymis++) {
    i = ymis_index[i_ymis]-1;
    p_tmp = expit(theta[0] + age[i] * theta[1]); // p_tmp on [Nhav,M]-matriisi (jokaiselle i = 1,...,M lasketaan imputointi)
    smokes_hav[i] = ::Rf_rbinom(1, p_tmp);
  }
  //smokes_imp[sel] = rbinom(n=sum(sel),size=1,prob=p_tmp) #TODO: Tarvitaan M imputointia
  return(smokes_hav);
}
