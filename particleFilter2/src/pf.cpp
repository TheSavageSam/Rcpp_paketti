#include <Rcpp.h>
#include "cppSample.h"

RcppExport SEXP pf(SEXP df_, SEXP parlist_) {

try {          // or use BEGIN_RCPP macro
//take parameters in
Rcpp::DataFrame df = Rcpp::DataFrame(df_);
Rcpp::List parlist = Rcpp::List(parlist_);

//parameters
int nsim = Rcpp::as<int>(parlist["nsim"]);
double sigmax_init = Rcpp::as<double>(parlist["sigmax_init"]);
double initial_value = Rcpp::as<double>(parlist["init_sim"]);
double home_eff_init = Rcpp::as<double>(parlist["home_eff_init"]);
double scaler_par = 0.95; // Parametrien arvoja pienennetään tällä kertoimella jokaista päivää kohti
// Käytä: double R_pow(double x, double y); => ::R_pow(x,y);

//data
std::vector<int> goals = df["goals"];
std::vector<double> day = df["day"];
std::vector<double> home = df["home"];

// lasketaan montako päivää kahden ottelun välillä on:
std::vector<double> daysdifference(day.size());
std::vector<double> scale(day.size());
std::vector<double> rate(day.size());
daysdifference[0] = 1;
scale[0] = scaler_par;
rate[0] = sigmax_init * scale[0];
for(int i = 1; i < day.size(); i++) {
daysdifference[i] = day[i]-day[i-1];
scale[i] = ::R_pow(scaler_par,daysdifference[i]);
rate[i] = sigmax_init * scale[i];
}

int ncol = goals.size()+1;
int nrow = nsim;
double tarkkuus = 1/sigmax_init;
Rcpp::NumericMatrix zline_mat(nrow,ncol);
Rcpp::NumericMatrix w_tilde_log(nrow,ncol);
Rcpp::NumericMatrix w_tilde(nrow,ncol);
Rcpp::NumericMatrix w_line_mat(nrow,ncol);
Rcpp::NumericMatrix z_mat(nrow,ncol);

std::vector<double> w_sum(nrow);
//Asetetaan nollia täyteen: Tarviiko tehdä erikseen vai onko valmiiksi?

for(int i = 0; i < nrow; i++) {
double tmp1 = ::Rf_rgamma(initial_value*rate[0], 1/rate[0]);
zline_mat(i,0) = tmp1;
z_mat(i,0) = tmp1; //Tämä on sämplättyjen matriisi, mutta käytetään silti tätä näin syntaktisista syistä.

double tmp_w = ::Rf_dgamma(tmp1,3.0*rate[0], 1/rate[0], true);
w_tilde_log(i,0) = tmp_w;
w_tilde(i,0) = 1;
w_sum[i] = 0;
}


for(int j = 1; j < ncol; j++) {
for(int i = 0; i < nrow; i++) {
double tmp = ::Rf_rgamma(0.5*(z_mat(i,j-1)+goals[j-1])*rate[j-1], 1/rate[j-1]);
double proposal = tmp; //Rcpp::as<double>(tmp);
zline_mat(i,j) = proposal;

// probability weights of proposal distribution
double tmp_w1 = ::Rf_dgamma(tmp,0.5*(z_mat(i,j-1)+goals[j-1])*rate[j-1], 1/rate[j-1], true);

// TODO: probability weights of parameters change in time -distribution
double tmp_w2 = ::Rf_dgamma(proposal, z_mat(i,j-1)*rate[j-1], 1/rate[j-1], 1);

double heff = 1.0;
if (home[j-1] == 1) {
heff = 1.11;
} else {
heff = 2 - 1.11;
}
// probability weights of p(observations | parameters) -distribution (poisson)
double tmp_w3 = ::Rf_dpois(goals[j-1], heff * proposal, 1); // lisää data tähän goals[j-1]

w_tilde_log(i,j) = tmp_w2 + tmp_w3 - tmp_w1;
w_sum[i] = w_tilde_log(i,j); // TODO: Viimeksi tässä ei ollut summaa?
w_tilde(i,j) = ::exp(w_tilde_log(i,j));

}

std::vector<double> w_tilde_exp(nrow);
std::transform(w_sum.begin(),w_sum.end(),w_tilde_exp.begin(),exp); // lasketaan log-painoista tavan painot
std::vector<double> w_line(nrow);
// Suoritetaan normeeraus käyttäen STL-funktiota. Normeerattu todennäköisyys on nyt vektori z.
double summa = std::accumulate(w_tilde_exp.begin(),w_tilde_exp.end(),0.0);
//if (summa > 0) {
std::transform(w_tilde_exp.begin(), w_tilde_exp.end(), w_line.begin(), std::bind1st(std::multiplies<double>(),1/summa));
//} else {
//  z[0] <- 1; // TODO: Toimiikohan?
//}
//TODO: tarkasta, että summa ei ole 0

//sämpläys
std::vector<double> sample_from(nrow);  // = Rcpp::as<std::vector<double> >(zline_mat(Rcpp::_,j));
for(int i = 0; i < nrow; i++) {
sample_from[i] = zline_mat(i,j);
w_line_mat(i,j) = w_line[i];
}
std::vector<double> sampled_values = cppSample(nrow, sample_from, w_line);
for(int i = 0; i < nrow; i++) {
z_mat(i,j) = sampled_values[i];
}
if (j%5==0) {
  Rprintf("\v%d ajanhetkeä simuloitu...\n",j);
  R_FlushConsole();
  R_ProcessEvents();
}

}

return(Rcpp::List::create(Rcpp::Named("Zline.mat") = zline_mat,
Rcpp::Named("log(Wtilde.mat)") = w_tilde_log,
Rcpp::Named("Wtilde.mat") = w_tilde,
Rcpp::Named("Wline.mat") = w_line_mat,
Rcpp::Named("Z.mat") = z_mat));

} catch( std::exception &ex ) {    // or use END_RCPP macro
forward_exception_to_r( ex );
} catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}