// Lisätään kommenttiin, miten vastaava rivi alustettaisiin Rcpp::Armadillossa!
#include "RcppArmadillo.h"
#include "cppSample.h"
#include "rcppSample.h"
#include "rk.h"
#include "k.h"
#include "rinitials.h"
#include "apply_mean.h"

using namespace Rcpp; // Tämä käyttöön, niin ei tarvi aina kirjoittaa Rcpp::
using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List popmc2(Rcpp::DataFrame df_, Rcpp::NumericVector ymis_index_, Rcpp::List parlist_, Rcpp::Function theta_dens, Rcpp::Function theta_prop, Rcpp::Function d_k_log, Rcpp::Function d_g_log) {

try {          // or use BEGIN_RCPP macro
//take parameters in
DataFrame df = DataFrame(df_);
NumericVector ymis_index = NumericVector(ymis_index_);
List parlist = List(parlist_);
// proposal functions and likelihood functions
Function thetadens("theta_dens"); // parametrien tiheysfunktio
Function thetaprop("theta_prop"); // parametrien ehdotusjakauma
Function dklog("d_k_log"); // 
Function dglog("d_g_log"); // 

RNGScope scp; //Rcpp::RNGScope

//parameters
NumericVector x_par = parlist["x"];
NumericVector mu = parlist["mu"];
NumericMatrix sdmat = parlist["sdmat"];

//data
//vector<int> smokes_hav = df["smokes_hav"];
IntegerVector smokes_hav = df["smokes_hav"];
NumericVector age = df["age"];
NumericVector x = df["x"];
IntegerVector osal = df["osal"];
IntegerVector ymis_index2 = as<IntegerVector>(ymis_index); // Eri mittainen, lyhyempi, indeksivektori, joka kertoo mitkä ovat puuttuvia

int ncol = smokes_hav.size()+1; // ilmeisesti initial imputations on se eka??

int Tt = 10; // otetaan argumentteina (montako ajanhetki-iteraatiota)
int M = as<int>(parlist["nsim"]);  // miten iso populaatio
int n = smokes_hav.size(); // datasta montako havaintoa
NumericMatrix theta_df(2, M); theta_df.fill(0); //theta_df<-matrix(c(-0.8,-0.1)+rnorm(2*M,sd=0.000001),ncol=2,nrow=M,byrow=T)
arma::cube res_theta_array = arma::zeros<arma::cube>(2, M, Tt); // Tt muutettu viimeiseksi indeksiksi, jotta voi käyttä .slice()-operaatiota //res_theta_array<-array(dim=c(Tt,dim(theta_df))) #
NumericMatrix imputed_data(M, n); //imputed_data<-matrix(NA,ncol=M,nrow=length(df$smokes_hav)) # 
arma::mat theta_all(M, 2); //theta_all<-matrix(NA,nrow=M,ncol=2) # 
NumericVector dens(M); //dens<-vector("double",M) # 
NumericMatrix r(Tt, M); //r<-matrix(NA,nrow=Tt,ncol=M) # 
NumericMatrix w(Tt, M); //w<-matrix(NA,nrow=Tt,ncol=M) # 
arma::mat mu_mat(M, 2); //mu_mat<-matrix(NA,nrow=M,ncol=2) # 
mu_mat.zeros();
arma::cube vcov_arr = arma::zeros<arma::cube>(2,2,M); //oli (Tt,2,2); //vcov_arr<-array(NA,dim=c(M,2,2)) # 

NumericMatrix nn(M, M); //nn<-matrix(NA,nrow=M,ncol=M) # 
NumericMatrix dd(M, M); //dd<-matrix(NA,nrow=M,ncol=M) # 

NumericVector rr(M);

for (int t=0; t < Tt; t++) {
    // Askel (a) (M kpl toistoja)
    for (int i=0; i < M; i++) {
      NumericVector apu = theta_df(_,i);
      imputed_data(_,i) = rk(smokes_hav, age, ymis_index2, apu );        //TODO! imputed_data[,i]<-k(i=i,df=df,theta_df=theta_df) # 
      //IntegerVector rk(IntegerVector smokes_hav, NumericVector age, IntegerVector ymis_index, NumericVector theta);
      List tmp = thetaprop(smokes_hav, age);        //TODO! tmp<-theta_prop(1,smokes_imp=imputed_data[,i],df) # theta_prop<-function(smokes_imp,age)
      
      arma::vec tmp2 = as<arma::vec>(tmp[0]); // tmp[0] on vektori, pituus 2
      theta_all.row(i) = tmp2.t();
      //theta_all(i,0) = tmp2(0,0);        //theta_all[i,1:2]<-tmp[[1]] # 
      //theta_all(i,1) = tmp2(0,1);        //theta_all[i,1:2]<-tmp[[1]] # 
      dens(i) = tmp[1]; //tmp[1] on double // dens[i]<-tmp[[2]] #
      arma::vec tmp3 = as<arma::vec>(tmp[2]); // tmp[2] on vektori, pituus 2 Rcpp::Numeric
      mu_mat.row(i) = tmp3.t();   //TODO: Voiko virhe olla täällä?        //mu_mat[i,]<-tmp[[3]] # 
//      mu_mat(i,0) = tmp3(0);
//      mu_mat(i,1) = tmp3(1);
      
      vcov_arr.slice(i) = as<arma::mat>(tmp[3]);//apu;
    }

//    return(List::create(Named("mu_mat") = mu_mat,
//                        Named("vcov") = vcov_arr,
//                        Named("theta_all") = theta_all));
    // Askel (b)
    
    // Tarvitaan MxM matriisit full data likelihoodista, ehdotusjakaumasta kdim
    // dimensiot: imputointi, theta-parametri
    double log_dens_pi = 1, log_dens_g = 1, log_dens_k_ll = 1, log_dens_k_li = 1;
    for (int l=0; l < M; l++) {
      for (int i=0; i < M; i++) {
        
        // lasketaan pi-jakauman tiheys
        arma::rowvec rv = theta_all.row(i);
        arma::vec rv_t = rv.t();
        arma::vec rv_t2 = theta_all.row(i).t(); // oli .t() vielä lisäksi
        arma::mat A = vcov_arr.slice(l);
        arma::vec muu = mu_mat.row(l); // 0x0 vectori
        return(List::create(Named("rv") = rv,
                            Named("rv_t") = rv_t,
                            Named("rv_t2") = rv_t2,
                            Named("A") = A,
                            Named("muu") = mu_mat.row(l)));
        
        log_dens_pi = as<double>(thetadens(Rcpp::as<NumericVector>(wrap(theta_all.row(i).t())), Rcpp::as<NumericVector>(wrap(muu)), Rcpp::as<NumericMatrix>(wrap(A)) )); //log_dens_pi<-theta_dens(theta_all[i,1:2],mu_mat[l,],vcov_arr[l,,]) # 
        log_dens_g = as<double>(dglog(Rcpp::as<NumericVector>(wrap(imputed_data(_,l))),Rcpp::as<NumericVector>(wrap(theta_all.row(i).t())),age)); //log_dens_g <-d_g_log(smokes=imputed_data[,l],theta_all[i,],df) # 
        //TODO: dklog tarvitsee indeksivektorin siitä, mitkä ovat puuttuvaa tietoa!
        log_dens_k_ll = as<double>(dklog(Rcpp::as<NumericVector>(wrap(imputed_data(_,l))), Rcpp::as<NumericVector>(wrap(theta_df(i,_))), age, ymis_index)); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta //log_dens_k_ll <-d_k_log(smokes=imputed_data[,l],theta_df[i,],df) #
        log_dens_k_li = as<double>(dklog(Rcpp::as<NumericVector>(wrap(imputed_data(_,l))), Rcpp::as<NumericVector>(wrap(theta_df(l,_))), age, ymis_index)); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta //log_dens_k_li <-d_k_log(smokes=imputed_data[,l],theta_df[l,],df) #
        
        nn(i,l) = exp(log_dens_g-log_dens_k_ll-log(M)); // log()-funktio tarvisee ehkä alkuun #include <bits/stdc++.h> //nn[i,l]<-log_dens_g-log_dens_k_ll-log(M) # 
        dd(i,l) = exp(log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M)); //dd[i,l]<-log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M) # 
      }
      Rprintf("Iteraatio %d:%d suoritettu.",t,l); R_FlushConsole(); R_ProcessEvents(); //print(paste("Iteraatio ",t,":",l," suoritettu.",sep="")) # 
    }
    
    double avg_nn = average(nn); // TODO: Miten tästä saisi luokkien NumericVector ja NumericMatrix metodin?
    double avg_dd = average(dd);
    //NumericVector unscaled = row_means(exp(nn - mean(nn))); // jos ei toimi, niin korvaa mean(nn) tällä mean(rowMeans(nn)) //unscaled<-rowMeans(exp(nn-mean(nn))) # 
    NumericVector unscaled = row_means(nn); // jos ei toimi, niin korvaa mean(nn) tällä mean(rowMeans(nn)) //unscaled<-rowMeans(exp(nn-mean(nn))) # 
    double sm = std::accumulate(nn.begin(), nn.end(), 0.0);
    double avg = sm / nn.size();
    double avg3 = average2(unscaled);
    
    // alkuperäiset luvut saadaan unscaled*exp(mean(nn))
    NumericVector jakaja = row_means(dd); //jakaja<-rowMeans(exp(dd)) # 
    for (int i=0;i<rr.size();i++) {
      rr(i)= unscaled(i) / jakaja(i); //rr<-unscaled/jakaja # 
    }
    NumericVector ww = rr / sum(rr); //ww<-rr/sum(rr) # 
    w(t,_) = ww; //w[t,]<-ww # 
    // Askel (c)
    Function seq = Environment("package:base")["seq"];
    NumericVector sekv = seq(1,M,1);
    // rcppSample on versio cppSample-funktiosta, joka käytää NumericVector-muotoa!
    NumericVector sampled(M);
    sampled = rcppSample(M, sekv, ww); //sampled<-sample.int(M,M,replace=T,prob=ww) # 
    // theta_df = theta_all[sampled,_]; //theta_df<-theta_all[sampled,] # //TODO: Tämä R:n tapa samplata ei toimi!
    for(int i=0;i<sampled.size();i++) {
      theta_df(1,i) = theta_all(sampled(i),1);
      theta_df(2,i) = theta_all(sampled(i),2);
//      theta_df(1,i) = theta_all(i,1);
//      theta_df(2,i) = theta_all(i,2);
    }
    
    // tallennetaan theta:t taulukkoon
    //res_theta_array.slice(t) = theta_df; //res_theta_array[t,,]<-theta_all[sampled,] # 
    for(int i=0;i<sampled.size();i++) { //i=1,...,M
      res_theta_array(1,i,t) = theta_df(1,i);
      res_theta_array(2,i,t) = theta_df(2,i);
    }
    Rprintf("ITERAATIO % suoritettu.",t); //print(paste("ITERAATIO ",t," suoritettu.",sep="")) # 
    
}
return(List::create(Named("imputed_data") = imputed_data,
                    Named("theta_all") = theta_all, //TODO: Tarkasta mikä vastaa mitäkin nimeä ja muunna nimet sopiviksi vasta, kun tiedät niiden olevan varmasti oikein!
                    Named("res_theta_array") = res_theta_array));
 
} catch( std::exception &ex ) {    // or use END_RCPP macro
  forward_exception_to_r( ex );
} catch(...) {
  ::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
} // End of function
