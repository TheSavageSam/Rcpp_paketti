rm(list=ls())
library(Rcpp)
package_name <- "PopMC"
setwd("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/")


unlink(paste("./",package_name,"/",sep=""),recursive=TRUE)

# malliin liittyvät funktiot:

rinitials_txt <- '
#include <Rcpp.h>
/* double a, double b */
double rinitials(double a, double b) {
  double res = ::Rf_rgamma(a,b);
  return( res );
}'
rinitialsh_txt <- '
// rinitials.h
double rinitials(double a, double b);
'


#---------------------------

cppsample_txt <- '
#include <Rcpp.h>
/* n_ = count of sampled values, values_ given values to sample from, probabilities_ = probabilities to use with sampling */
std::vector<double> cppSample(int n, std::vector<double> values, std::vector<double> probabilities) {

\t//Rcpp::RNGScope scope;             // do we really need this?

\tstd::vector<double> cumsum(probabilities.size());
\tstd::partial_sum(probabilities.begin(),probabilities.end(),cumsum.begin(),std::plus<double>());

\t//Rcpp::NumericVector rand_unif = Rcpp::runif(n, Rcpp::Named("min") = 0,Rcpp::Named("max") = maxint);
\tstd::vector<double> rand(n);
\tfor(int i = 0; i < n; i++ ) {
\t\trand[i] = ::Rf_runif(0,1);
\t}

\tstd::vector<double> sampled_values(n);
\t// tehdään eka luku erikseen ja loput loopilla
\tfor(int k = 0; k < rand.size(); k++) {
\t\tif (rand[k] <= cumsum[0]) {
\t\t\tsampled_values[k] = values[0];
\t\t}
\t\tfor(int i = 1; i < cumsum.size(); i++) {
\t\t\tif(rand[k] > cumsum[i-1] & rand[k] <= cumsum[i]) {
\t\t\t\tsampled_values[k] = values[i];
\t\t\t}
\t\t}
\t}
\treturn( sampled_values );
}'
cppsampleh_txt <- '
// cppSample.h
std::vector<double> cppSample(int n, std::vector<double> values, std::vector<double> probabilities);
'
k_txt <- '
#include <Rcpp.h>
double expit(double x) { return(exp(x)/(1+exp(x))); }

std::vector<int> k(std::vector<int> smokes_hav, std::vector<double> age, std::vector<int> ymis_index, std::vector<double> theta) {
\t// df: data
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
'
kh_txt <- '
double expit(double x);
std::vector<int> k(std::vector<int> smokes_hav, std::vector<double> age, std::vector<int> ymis_index, std::vector<double> theta);'


txt2 <- '#include <Rcpp.h>
#include "cppSample.h"
#include "k.h"
#include "rinitials.h"

RcppExport SEXP popmc(SEXP df_, SEXP ymis_index_, SEXP parlist_, SEXP theta_dens, SEXP theta_prop) {

try {          // or use BEGIN_RCPP macro
//take parameters in
Rcpp::DataFrame df = Rcpp::DataFrame(df_);
Rcpp::NumericVector ymis_index = Rcpp::NumericVector(ymis_index_);
Rcpp::List parlist = Rcpp::List(parlist_);
//Rcpp::List functionlist = Rcpp::List(functions_list_);
Rcpp::Function tdens("theta_dens"); // parametrien tiheysfunktio
Rcpp::Function tprop("theta_prop"); // parametrien ehdotusjakauma

Rcpp::RNGScope scp;

//parameters
int nsim = Rcpp::as<int>(parlist["nsim"]);
Rcpp::NumericVector x_par = parlist["x"];
Rcpp::NumericVector mu = parlist["mu"];
Rcpp::NumericMatrix sdmat = parlist["sdmat"];

double out = Rcpp::as<double>(tdens(x_par,mu,sdmat));

//data
std::vector<int> smokes_hav = df["smokes_hav"];
std::vector<double> age = df["age"];
std::vector<double> x = df["x"];
std::vector<int> osal = df["osal"];
std::vector<int> ymis_index2 = Rcpp::as<std::vector<int> >(ymis_index); // Eri mittainen, lyhyempi, indeksivektori, joka kertoo mitkä ovat puuttuvia


int ncol = smokes_hav.size()+1;
int nrow = nsim; // partikkelien määrä TAI populaation koko

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

// Mallin parametrit, nimeltään thetam
std::vector<double> thetam;
thetam.push_back(0.125);
thetam.push_back(0.85);

// TODO: Iteroi puuttuvien z ylitse ja imputoi ne kukin nsim=ncol kertaa.
smokes_hav = k(smokes_hav, age, ymis_index2, thetam);
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
}'

tmp<-'
for(int i = 0; i < nrow; i++) {
\t// TODO: Askel (a) Simuloi sekä imputoinnit että parametrien alkuarvot
\tdouble initvalues = rinitials(initial_value*rate[0], 1/rate[0]);
\tzline_mat(i,0) = initvalues;
\tz_mat(i,0) = initvalues; //Tämä on sämplättyjen matriisi, mutta käytetään silti tätä näin syntaktisista syistä.

\tdouble tmp_w = ::Rf_dgamma(initvalues,3.0*rate[0], 1/rate[0], true);
\tw_tilde_log(i,0) = tmp_w;
\tw_tilde(i,0) = 1;
\tw_sum[i] = 0;
}


for(int j = 1; j < ncol; j++) {
\tfor(int i = 0; i < nrow; i++) {
\t\tdouble tmp = ::Rf_rgamma(0.5*(z_mat(i,j-1)+smokes_hav[j-1])*rate[j-1], 1/rate[j-1]);
\t\tdouble proposal = tmp; //Rcpp::as<double>(tmp);
\t\tzline_mat(i,j) = proposal;
\t\t
\t\t// probability weights of proposal distribution
\t\tdouble loglik_proposal = ::Rf_dgamma(tmp,0.5*(z_mat(i,j-1)+smokes_hav[j-1])*rate[j-1], 1/rate[j-1], true);
\t\t
\t\t// TODO: probability weights of parameters change in time -distribution
\t\tdouble loglik_model = ::Rf_dgamma(proposal, z_mat(i,j-1)*rate[j-1], 1/rate[j-1], 1);
\t\t
\t\t// probability weights of p(observations | parameters) -distribution (poisson)
\t\tdouble loglik_obsmodel = ::Rf_dpois(smokes_hav[j-1], proposal, 1);
\t\t
\t\tw_tilde_log(i,j) = loglik_model + loglik_obsmodel - loglik_proposal;
\t\tw_sum[i] = w_tilde_log(i,j); // TODO: Viimeksi tässä ei ollut summaa?
\t\tw_tilde(i,j) = ::exp(w_tilde_log(i,j));
\t\t
\t}
\t
\tstd::vector<double> w_tilde_exp(nrow);
\tstd::transform(w_sum.begin(),w_sum.end(),w_tilde_exp.begin(),exp); // lasketaan log-painoista tavan painot
\tstd::vector<double> w_line(nrow);
\t// Suoritetaan normeeraus käyttäen STL-funktiota. Normeerattu todennäköisyys on nyt vektori z.
\tdouble summa = std::accumulate(w_tilde_exp.begin(),w_tilde_exp.end(),0.0);
\t//if (summa > 0) {
\tstd::transform(w_tilde_exp.begin(), w_tilde_exp.end(), w_line.begin(), std::bind1st(std::multiplies<double>(),1/summa));
\t//} else {
\t//  z[0] <- 1; // TODO: Toimiikohan?
\t//}
\t//TODO: tarkasta, että summa ei ole 0
\t
\t//sämpläys
\tstd::vector<double> sample_from(nrow);  // = Rcpp::as<std::vector<double> >(zline_mat(Rcpp::_,j));
\tfor(int i = 0; i < nrow; i++) {
\t\tsample_from[i] = zline_mat(i,j);
\t\tw_line_mat(i,j) = w_line[i];
\t}
\tstd::vector<double> sampled_values = cppSample(nrow, sample_from, w_line);
\tfor(int i = 0; i < nrow; i++) {
\t\tz_mat(i,j) = sampled_values[i];
\t}
\tif (j%5==0) {
\t\tRprintf("%d ajanhetkeä simuloitu...",j);
\t\tR_FlushConsole();
\t\tR_ProcessEvents();
\t}
\t
}
'
rfunc_txt <- '
popmc <- function(df_,par_,theta_dens) {
#	mis<-apply(df_,2,function(x) any(is.na(x)))
#	newname<-paste(names(mis),"_mis",sep="")[mis]
#	newvalues<-as.integer(is.na(df_[,mis]))
#	df_[,dim(df_)[2]+1]<-newvalues
#	names(df_)<-c(names(df_),newname)
\ty <- .Call("popmc",df_=df_,ymis_index_=which(is.na(df_$smokes_hav)),parlist_=par_,theta_dens=theta_dens,PACKAGE="PopMC")
#\tfor(i in 1:(length(y)-1)) {
#\t\ty[[i]] <- y[[i]][,-1]
#\t}
\ty
}'

mavfunc_txt <- 'mean_and_var <- function(Wline.mat,Zline.mat) {
x_1 <- colSums(Wline.mat * Zline.mat)
x_2 <- colSums(Wline.mat * Zline.mat^2)
v <- x_2 - (x_1)^2
list("means"=x_1,"vars"=v)
}'

plotfunc_txt <- '
plot_pf_results <- function(df,pf_mav,probs,...) {

# TODO: Voisiko käyttäen gamma-jakauman oletuksia laskea luottovälit?
#lasketaan cv, joka on gammalle vakio (voi toki vaihdella ajassa parametrien kehittyessä)
cv <- sqrt(pf_mav$vars)/pf_mav$means
b <- 1/cv
a <- pf_mav$means * b
l1 <- vector("list",length(probs))
for(i in 1:length(probs)) {
l1[[i]] <- qgamma(probs[i],a,b)
}
par(mai=c("bottom"=2.2, "left"=0.7, "top"=0.7, "right"=0.7))
plot.default(x=df$day[df$home==1],y=df$goals[df$home==1],col="blue",ylim=range(l1,df$goals),xlim=range(df$day),xaxt = "n",xlab="",...)
points(x=df$day[df$home==0],y=df$goals[df$home==0],col="red")
axis(side=1, at=df$day, labels=rownames(df),las=2)
lines(x=df$day,y=pf_mav$means,col="green")
for(i in 1:length(probs)) {
lines(x=df$day,y=l1[[i]],lty=2,col="green")
}
}'

hello_example <- '
#include "rcpp_hello_world.h"

SEXP rcpp_hello_world(){
    using namespace Rcpp ;
    
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    //Numericvec y   = Numericvec::create( 0.0, 1.0 ) ;
    NumericVector vec(Dimension(2, 2, 2));
    for(int i = 0; i < 8; i++) {
      vec(i) = i+1;
    }
    List z            = List::create( x, vec ) ;
    
    return z ;
}'
hello_header <- '
#ifndef _PopMC_RCPP_HELLO_WORLD_H
#define _PopMC_RCPP_HELLO_WORLD_H

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call cant find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
 RcppExport SEXP rcpp_hello_world() ;

#endif
'

Rcpp.package.skeleton(name="PopMC")

unlink(paste("./",package_name,"/src/rcpp_hello_world.cpp",sep=""))
unlink(paste("./",package_name,"/src/rcpp_hello_world.h",sep=""))
unlink(paste("./",package_name,"/R/rcpp_hello_world.R",sep=""))
unlink(paste("./",package_name,"/man/rcpp_hello_world.Rd",sep=""))

#cat(hello_example,file=paste("./",package_name,"/src/rcpp_hello_world.cpp",sep=""))
#cat(hello_header,file=paste("./",package_name,"/src/rcpp_hello_world.h",sep=""))
cat(cppsample_txt,file=paste("./",package_name,"/src/cppSample.cpp",sep=""))
cat(cppsampleh_txt,file=paste("./",package_name,"/src/cppSample.h",sep=""))
cat(k_txt,file=paste("./",package_name,"/src/k.cpp",sep=""))
cat(kh_txt,file=paste("./",package_name,"/src/k.h",sep=""))

cat(rinitials_txt,file=paste("./",package_name,"/src/rinitials.cpp",sep=""))
cat(rinitialsh_txt,file=paste("./",package_name,"/src/rinitials.h",sep=""))
cat(txt2,file=paste("./",package_name,"/src/popmc.cpp",sep=""))
cat(rfunc_txt,file=paste("./",package_name,"/R/popmc.R",sep=""))
cat(mavfunc_txt,file=paste("./",package_name,"/R/plot_pf_results.R",sep=""))
cat(plotfunc_txt,file=paste("./",package_name,"/R/mean_and_var.R",sep=""))

compileAttributes(pkgdir = "./PopMC/",verbose=TRUE)

