// Lisätään kommenttiin, miten vastaava rivi alustettaisiin Rcpp::Armadillossa!
#include <Rcpp.h>
#include <armadillo>
#include "cppSample.h"
#include "rk.h"
#include "k.h"
#include "rinitials.h"
#include "apply_mean.cpp"

using namespace Rcpp; // Tämä käyttöön, niin ei tarvi aina kirjoittaa Rcpp::
using namespace std;

RcppExport SEXP popmc(SEXP df_, SEXP ymis_index_, SEXP parlist_, SEXP theta_dens, SEXP theta_prop, SEXP d_k_log, SEXP d_g_log) {
  
//try {          // or use BEGIN_RCPP macro
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
vector<int> smokes_hav = df["smokes_hav"];
NumericVector smokes_hav_r = df["smokes_hav"];
vector<double> age = df["age"];
vector<double> x = df["x"];
vector<int> osal = df["osal"];
vector<int> ymis_index2 = as<vector<int> >(ymis_index); // Eri mittainen, lyhyempi, indeksivektori, joka kertoo mitkä ovat puuttuvia

int ncol = smokes_hav.size()+1; // ilmeisesti initial imputations on se eka??

int Tt(1); // otetaan argumentteina (montako ajanhetki-iteraatiota)
int M = as<int>(parlist["nsim"]);  // miten iso populaatio
int n = smokes_hav.size(); // datasta montako havaintoa
NumericMatrix theta_df(2, M, fill::zeros); //theta_df.fill(0); //theta_df<-matrix(c(-0.8,-0.1)+rnorm(2*M,sd=0.000001),ncol=2,nrow=M,byrow=T)
arma::dcube res_theta_array(2, M, Tt, fill::zeros); // Tt muutettu viimeiseksi indeksiksi, jotta voi käyttä .slice()-operaatiota //res_theta_array<-array(dim=c(Tt,dim(theta_df))) #
NumericMatrix imputed_data(M, n, fill::zeros); //imputed_data<-matrix(NA,ncol=M,nrow=length(df$smokes_hav)) # 
NumericMatrix theta_all(M, 2, fill::zeros); //theta_all<-matrix(NA,nrow=M,ncol=2) # 
NumericVector dens(M); //dens<-vector("double",M) # 
NumericMatrix r(Tt, M); //r<-matrix(NA,nrow=Tt,ncol=M) # 
NumericMatrix w(Tt, M); //w<-matrix(NA,nrow=Tt,ncol=M) # 
NumericMatrix mu_mat(M, 2); //mu_mat<-matrix(NA,nrow=M,ncol=2) # 
arma::dcube vcov_arr(2,2,Tt); //oli (Tt,2,2); //vcov_arr<-array(NA,dim=c(M,2,2)) # 

NumericMatrix nn(M, M); //nn<-matrix(NA,nrow=M,ncol=M) # 
NumericMatrix dd(M, M); //dd<-matrix(NA,nrow=M,ncol=M) # 


for (int t=0; t < Tt; t++) {
    // Askel (a) (M kpl toistoja)
    for (int i=0; i < M; i++) {
      imputed_data(_,i) = rk(smokes_hav_r, age, ymis_index, theta_df(_,i));        //TODO! imputed_data[,i]<-k(i=i,df=df,theta_df=theta_df) # 
      //IntegerVector rk(IntegerVector smokes_hav, NumericVector age, IntegerVector ymis_index, NumericVector theta);
      List tmp = thetaprop(smokes_hav_r, age);        //TODO! tmp<-theta_prop(1,smokes_imp=imputed_data[,i],df) # theta_prop<-function(smokes_imp,age)
      theta_all(i,_) = tmp[0];        //theta_all[i,1:2]<-tmp[[1]] # 
      dens(i) = tmp[1];               //dens[i]<-tmp[[2]] # 
      mu_mat(i,_) = tmp[2];           //mu_mat[i,]<-tmp[[3]] # 
      vcov_arr.slice(i) = tmp[3];     //vcov_arr[i,,]<-tmp[[4]] #  // TAI vcov_arr.slice(i) = tmp[3]; mutta tällöin pitää määrittää aluksi arma::cube vcov_arr(2,2,Tt); eikä (Tt,2,2)
    }
    
    // Askel (b)
    
    // Tarvitaan MxM matriisit full data likelihoodista, ehdotusjakaumasta kdim
    // dimensiot: imputointi, theta-parametri
    double log_dens_pi = 1, log_dens_g = 1, log_dens_k_ll = 1, log_dens_k_li = 1;
    for (int l=0; l < M; l++) {
      for (int i=0; i < M; i++) {
        
        // lasketaan pi-jakauman tiheys
        log_dens_pi = as<double>(thetadens(theta_all(i,_), mu_mat(l,_), vcov_arr.slice(l))); //log_dens_pi<-theta_dens(theta_all[i,1:2],mu_mat[l,],vcov_arr[l,,]) # 
        log_dens_g = as<double>(dglog(imputed_data(_,l),theta_all(i,_),df)); //log_dens_g <-d_g_log(smokes=imputed_data[,l],theta_all[i,],df) # 
        log_dens_k_ll = as<double>(dklog(imputed_data(_,l),theta_df(i,_),df)); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta //log_dens_k_ll <-d_k_log(smokes=imputed_data[,l],theta_df[i,],df) #
        log_dens_k_li = as<double>(dklog(imputed_data(_,l),theta_df(l,_),df)); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta //log_dens_k_li <-d_k_log(smokes=imputed_data[,l],theta_df[l,],df) #
        
        nn(i,l) = log_dens_g-log_dens_k_ll-log(M); // log()-funktio tarvisee ehkä alkuun #include <bits/stdc++.h> //nn[i,l]<-log_dens_g-log_dens_k_ll-log(M) # 
        dd(i,l) = log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M); //dd[i,l]<-log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M) # 
      }
      Rprintf("Iteraatio %d:%d suoritettu.",t,l); R_FlushConsole(); R_ProcessEvents(); //print(paste("Iteraatio ",t,":",l," suoritettu.",sep="")) # 
    }
    
    NumericVector unscaled = row_means(exp(nn - mean(nn))); // jos ei toimi, niin korvaa mean(nn) tällä mean(rowMeans(nn)) //unscaled<-rowMeans(exp(nn-mean(nn))) # 
    // alkuperäiset luvut saadaan unscaled*exp(mean(nn))
    NumericVector jakaja = row_means(exp(dd)); //jakaja<-rowMeans(exp(dd)) # 
    Numericvector rr = unscaled / jakaja; //rr<-unscaled/jakaja # 
    NumericVector ww = rr / sum(rr); //ww<-rr/sum(rr) # 
    w(t,_) = ww; //w[t,]<-ww # 
    // Askel (c)
    // rcppSample on versio cppSample-funktiosta, joka käytää NumericVector-muotoa!
    NumericVector sampled = rcppSample(M, seq(1,M), ww); //sampled<-sample.int(M,M,replace=T,prob=ww) # 
    theta_df = theta_all(sampled,_); //theta_df<-theta_all[sampled,] # 
    // tallennetaan theta:t taulukkoon
    res_theta_array.slice(t) = theta_df; //res_theta_array[t,,]<-theta_all[sampled,] # 
    
    Rprintf("ITERAATIO % suoritettu.",t); //print(paste("ITERAATIO ",t," suoritettu.",sep="")) # 
    
  }
