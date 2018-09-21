#include <Rcpp.h>
#include "cppSample.h"
#include "k.h"
#include "rinitials.h"

using namespace Rcpp; // Tämä käyttöön, niin ei tarvi aina kirjoittaa Rcpp::
using namespace std;

RcppExport SEXP popmc(SEXP df_, SEXP ymis_index_, SEXP parlist_, SEXP theta_dens, SEXP theta_prop, SEXP d_k_log, SEXP d_g_log) {

try {          // or use BEGIN_RCPP macro
//take parameters in
DataFrame df = DataFrame(df_);
NumericVector ymis_index = NumericVector(ymis_index_);
List parlist = List(parlist_);
// proposal functions and likelihood functions
Function tdens("theta_dens"); // parametrien tiheysfunktio
Function tprop("theta_prop"); // parametrien ehdotusjakauma
Function dklog("d_k_log"); // 
Function dglog("d_g_log"); // 

RNGScope scp; //Rcpp::RNGScope

//parameters
int nsim = as<int>(parlist["nsim"]);
NumericVector x_par = parlist["x"];
NumericVector mu = parlist["mu"];
NumericMatrix sdmat = parlist["sdmat"];

//data
vector<int> smokes_hav = df["smokes_hav"];
NumericVector smokes_hav_r = df["smokes_hav"];
vector<double> age = df["age"];
vector<double> x = df["x"];
vector<int> osal = df["osal"];
vector<int> ymis_index2 = as<vector<int> >(ymis_index); // Eri mittainen, lyhyempi, indeksivektori, joka kertoo mitkä ovat puuttuvia

int ncol = smokes_hav.size()+1;
int nrow = nsim; // partikkelien määrä TAI populaation koko

// kutsutaan funktiota theta_dens
double td_out = as<double>(tdens(x_par,mu,sdmat));
print(wrap(td_out));

// Mallin parametrit, nimeltään thetam
vector<double> thetam;
thetam.push_back(0.125);
thetam.push_back(0.85);

// TODO: Iteroi puuttuvien z ylitse ja imputoi ne kukin nsim=ncol kertaa.
smokes_hav = k(smokes_hav, age, ymis_index2, thetam);
smokes_hav_r = wrap(smokes_hav);
// kutsutaan funktiota theta_prop 
//<NumericVector>

List tpout = List(tprop(smokes_hav,df["age"])); //uutta

// partikkelisuotimelle: 
/*
*NumericMatrix zline_mat(nrow,ncol);   // ehdotetut partikkelit
*NumericMatrix w_tilde_log(nrow,ncol); // painot (logaritmi-asteikko)
*NumericMatrix w_tilde(nrow,ncol);     // painot (tavallinen asteikko)
*NumericMatrix w_line_mat(nrow,ncol);  // normalisoidut painot
*NumericMatrix z_mat(nrow,ncol);       // hyväksytyt partikkelit
*/
// population monte carlolle:
NumericMatrix zline_mat(nrow,ncol);   // ehdotetut imputoinnit puuttuville z (erona on se, että vain osa havainnoista puuttuu, joten ne on käsiteltävä)
NumericMatrix theta(nrow,ncol);       // parametrit (M kertaa parametrien määrä)
NumericMatrix r_mat_log(nrow,ncol);   // painot (logaritmi-asteikko)
NumericMatrix r_mat(nrow,ncol);       // painot (tavallinen asteikko)
NumericMatrix w_mat(nrow,ncol);       // normalisoidut painot
NumericMatrix z_mat(nrow,ncol);       // hyväksytyt partikkelit

// TODO: Generate M number of realisations (population of realisations) for parameters called theta. Theta has two or more parameters, so this must be a matrix.
for (int i=0;i<nrow;i++) {
  // Lets generate a single realisation for theta to get initial values.
  List tpout = List(tprop(smokes_hav,df["age"])); //uutta
  NumericVector tp = tpout["out"];
  theta(i, _) = tp; // Rcpp::_
  // Simulate imputations from k-function.
  smokes_hav = k(smokes_hav, age, ymis_index2, as<vector<double> >(tp)); // argumentit on std-muotoisia! TODO: Mieti muut argumentit kuntoon!
  smokes_hav_r = wrap(smokes_hav); // Rcpp::wrap()
  
  
}
vector<double> w_sum(nrow);


//for(int i_ymis = 0; i_ymis < ymis_index.size(); i_ymis++) {
//smokes_hav[ymis_index2[i_ymis]-1] = ::Rf_rbinom(1, 0.5); // imputoidaan kolikonheitolla!
//	Rprintf("Puuttuvat tiedot: %d",smokes_hav[ymis_index2[i_ymis]-1]);
//	R_FlushConsole();
//	R_ProcessEvents();
//}
// TODO: Simuloi nsim=ncol määrän verran parametreja theta.

// harjoitellaan R:llä kirjoitetun funktion kutsumista
//Funcktion k("k");
//k()
//k<-function(i,df,theta_df)




NumericVector zz(smokes_hav.begin(), smokes_hav.end());


return(List::create(Named("Zline.mat") = zline_mat,
Named("theta.mat") = theta,
Named("Rmat.mat") = r_mat,
Named("W.mat") = w_mat,
Named("Z.mat") = z_mat,
Named("smokes_hav") = zz));

} catch( std::exception &ex ) {    // or use END_RCPP macro
forward_exception_to_r( ex );
} catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}