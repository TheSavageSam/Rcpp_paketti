#include <Rcpp.h>
#include "cppSample.h"
#include "k.h"
#include "rinitials.h"

RcppExport SEXP popmc(SEXP df_, SEXP ymis_index_, SEXP parlist_, SEXP theta_dens, SEXP theta_prop, SEXP d_k_log, SEXP d_g_log) {

try {          // or use BEGIN_RCPP macro
//take parameters in
Rcpp::DataFrame df = Rcpp::DataFrame(df_);
Rcpp::NumericVector ymis_index = Rcpp::NumericVector(ymis_index_);
Rcpp::List parlist = Rcpp::List(parlist_);
// proposal functions and likelihood functions
Rcpp::Function tdens("theta_dens"); // parametrien tiheysfunktio
Rcpp::Function tprop("theta_prop"); // parametrien ehdotusjakauma
Rcpp::Function dklog("d_k_log"); // 
Rcpp::Function dglog("d_g_log"); // 

Rcpp::RNGScope scp;

//parameters
int nsim = Rcpp::as<int>(parlist["nsim"]);
Rcpp::NumericVector x_par = parlist["x"];
Rcpp::NumericVector mu = parlist["mu"];
Rcpp::NumericMatrix sdmat = parlist["sdmat"];

//data
std::vector<int> smokes_hav = df["smokes_hav"];
Rcpp::NumericVector smokes_hav_r = df["smokes_hav"];
std::vector<double> age = df["age"];
std::vector<double> x = df["x"];
std::vector<int> osal = df["osal"];
std::vector<int> ymis_index2 = Rcpp::as<std::vector<int> >(ymis_index); // Eri mittainen, lyhyempi, indeksivektori, joka kertoo mitkä ovat puuttuvia

int ncol = smokes_hav.size()+1;
int nrow = nsim; // partikkelien määrä TAI populaation koko

// kutsutaan funktiota theta_dens
double td_out = Rcpp::as<double>(tdens(x_par,mu,sdmat));
Rcpp::print(Rcpp::wrap(td_out));

// Mallin parametrit, nimeltään thetam
std::vector<double> thetam;
thetam.push_back(0.125);
thetam.push_back(0.85);

// TODO: Iteroi puuttuvien z ylitse ja imputoi ne kukin nsim=ncol kertaa.
smokes_hav = k(smokes_hav, age, ymis_index2, thetam);
smokes_hav_r = Rcpp::wrap(smokes_hav);
// kutsutaan funktiota theta_prop 
//<NumericVector>

Rcpp::List tpout = Rcpp::List(tprop(smokes_hav,df["age"])); //uutta

// partikkelisuotimelle: 
/*
*Rcpp::NumericMatrix zline_mat(nrow,ncol);   // ehdotetut partikkelit
*Rcpp::NumericMatrix w_tilde_log(nrow,ncol); // painot (logaritmi-asteikko)
*Rcpp::NumericMatrix w_tilde(nrow,ncol);     // painot (tavallinen asteikko)
*Rcpp::NumericMatrix w_line_mat(nrow,ncol);  // normalisoidut painot
*Rcpp::NumericMatrix z_mat(nrow,ncol);       // hyväksytyt partikkelit
*/
// population monte carlolle:
Rcpp::NumericMatrix zline_mat(nrow,ncol);   // ehdotetut imputoinnit puuttuville z (erona on se, että vain osa havainnoista puuttuu, joten ne on käsiteltävä)
Rcpp::NumericMatrix theta(nrow,ncol);       // parametrit
Rcpp::NumericMatrix r_mat_log(nrow,ncol);   // painot (logaritmi-asteikko)
Rcpp::NumericMatrix r_mat(nrow,ncol);       // painot (tavallinen asteikko)
Rcpp::NumericMatrix w_mat(nrow,ncol);       // normalisoidut painot
Rcpp::NumericMatrix z_mat(nrow,ncol);       // hyväksytyt partikkelit


std::vector<double> w_sum(nrow);


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




Rcpp::NumericVector zz(smokes_hav.begin(), smokes_hav.end());


return(Rcpp::List::create(Rcpp::Named("Zline.mat") = zline_mat,
Rcpp::Named("theta.mat") = theta,
Rcpp::Named("Rmat.mat") = r_mat,
Rcpp::Named("W.mat") = w_mat,
Rcpp::Named("Z.mat") = z_mat,
Rcpp::Named("smokes_hav") = zz));

} catch( std::exception &ex ) {    // or use END_RCPP macro
forward_exception_to_r( ex );
} catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}