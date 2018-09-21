# Lisätään kommenttiin, miten vastaava rivi alustettaisiin Rcpp::Armadillossa!
#using namespace Rcpp;
#integer Tt(1); // otetaan argumentteina (montako ajanhetki-iteraatiota)
#integer M(1);  // miten iso populaatio
#integer n(1); //data_y.size(); // datasta montako havaintoa
theta_df<-matrix(c(-0.8,-0.1)+rnorm(2*M,sd=0.000001),ncol=2,nrow=M,byrow=T) # Rcpp::NumericMatrix theta_df(2, M, fill::zeros); //theta_df.fill(0);
res_theta_array<-array(dim=c(Tt,dim(theta_df))) # arma::dcube res_theta_array(2, M, Tt, fill::zeros); // Tt muutettu viimeiseksi indeksiksi, jotta voi käyttä .slice()-operaatiota
imputed_data<-matrix(NA,ncol=M,nrow=length(df$smokes_hav)) # NumericMatrix imputed_data(M, n, fill::zeros);
theta_all<-matrix(NA,nrow=M,ncol=2) # NumericMatrix theta_all(M, 2, fill::zeros);
dens<-vector("double",M) # NumericVector dens(M);
r<-matrix(NA,nrow=Tt,ncol=M) # NumericMatrix r(Tt, M);
w<-matrix(NA,nrow=Tt,ncol=M) # NumericMatrix w(Tt, M);
mu_mat<-matrix(NA,nrow=M,ncol=2) # NumericMatrix mu_mat(M, 2);
vcov_arr<-array(NA,dim=c(M,2,2)) # arma::cube vcov_arr(2,2,Tt) //oli (Tt,2,2);

nn<-matrix(NA,nrow=M,ncol=M) # NumericMatrix nn(M, M);
dd<-matrix(NA,nrow=M,ncol=M) # NumericMatrix dd(M, M);


aika<-system.time({
  for(t in 1:Tt) { # for (int t=0; t < Tt; t++) {
    # // Askel (a) (M kpl toistoja)
    for(i in 1:M) { # for (int i=0; i < M; i++) {
      imputed_data[,i]<-k(i=i,df=df,theta_df=theta_df) # imputed_data(_,i) = k();
      tmp<-theta_prop(1,smokes_imp=imputed_data[,i],df) # List tmp = theta_prop();
      theta_all[i,1:2]<-tmp[[1]] # theta_all(i,_) = tmp[0];
      dens[i]<-tmp[[2]] # dens(i) = tmp[1];
      mu_mat[i,]<-tmp[[3]] # mu_mat(i,_) = tmp[2];
      vcov_arr[i,,]<-tmp[[4]] # vcov_arr.slice(i) = tmp[3]; TAI vcov_arr.slice(i) = tmp[3]; mutta tällöin pitää määrittää aluksi arma::cube vcov_arr(2,2,Tt); eikä (Tt,2,2)
    } #}
    
    # // Askel (b)
    
    #// Tarvitaan MxM matriisit full data likelihoodista, ehdotusjakaumasta kdim
    #// dimensiot: imputointi, theta-parametri
    # double log_dens_pi = 1, log_dens_g = 1, log_dens_k_ll = 1, log_dens_k_li = 1;
    for(l in 1:M) { # for (int l=0; l < M, l++) {
      for(i in 1:M) { # for (int i=0; i < M, i++) {
        
        #// lasketaan pi-jakauman tiheys
        log_dens_pi<-theta_dens(theta_all[i,1:2],mu_mat[l,],vcov_arr[l,,]) # log_dens_pi = theta_dens(theta_all(i,_), mu_mat(l,_), vcov_arr.slice(l));
        log_dens_g <-d_g_log(smokes=imputed_data[,l],theta_all[i,],df) # log_dens_g = d_g_log(imputed_data(_,l),theta_all(i,_),df);
        log_dens_k_ll <-d_k_log(smokes=imputed_data[,l],theta_df[i,],df) #log_dens_k_ll = d_k_log(imputed_data(_,l),theta_df(i,_),df); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta
        log_dens_k_li <-d_k_log(smokes=imputed_data[,l],theta_df[l,],df) #log_dens_k_li = d_k_log(imputed_data(_,l),theta_df(l,_),df); // tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta
        
        nn[i,l]<-log_dens_g-log_dens_k_ll-log(M) # nn(i,l) = log_dens_g-log_dens_k_ll-log(M); // log()-funktio tarvisee ehkä alkuun #include <bits/stdc++.h> 
        dd[i,l]<-log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M) # dd(i,l) = log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M);
      } # }
      print(paste("Iteraatio ",t,":",l," suoritettu.",sep="")) # Rprintf("Iteraatio %d:%d suoritettu.",t,l); R_FlushConsole(); R_ProcessEvents();
    } # }
    
    unscaled<-rowMeans(exp(nn-mean(nn))) # NumericVector unscaled = sugar::rowMeans(exp(nn - mean(nn))); // jos ei toimi, niin korvaa mean(nn) tällä mean(rowMeans(nn))
    #// alkuperäiset luvut saadaan unscaled*exp(mean(nn))
    jakaja<-rowMeans(exp(dd)) # NumericVector jakara = rowMeans(exp(dd));
    rr<-unscaled/jakaja # Numericvector rr = unscaled / jakaja;
    ww<-rr/sum(rr) # NumericVector ww = rr / sum(rr);
    w[t,]<-ww # w(t,_) = ww;
    #// Askel (c)
    #// rcppSample on versio cppSample-funktiosta, joka käytää NumericVector-muotoa!
    sampled<-sample.int(M,M,replace=T,prob=ww) # NumericVector sampled =  rcppSample(M, seq(1,M), ww);
    theta_df<-theta_all[sampled,] # theta_df = theta_all(sampled,_);
    #// tallennetaan theta:t taulukkoon
    res_theta_array[t,,]<-theta_all[sampled,] # res_theta_array.slice(t) = theta_df;
    
    print(paste("ITERAATIO ",t," suoritettu.",sep="")) # Rprint("ITERAATIO % suoritettu.",t)
    
  } # }
})