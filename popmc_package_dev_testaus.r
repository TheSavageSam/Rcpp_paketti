## compile in between
# cd /home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/
# sudo R CMD build PopMC && sudo R CMD INSTALL PopMC

library(PopMC)

setwd("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/testialue")
load(file="testidata.Rdata",verbose=T)

head(df,20)
names(df)
mu<-c(0.12,0.85)
M<-matrix(c(1,0.12,0.12,1),2,2)

theta_prop<-function(smokes_imp,age) {
  # smokes_imp: imputoitu tupakointivektori
  # age: ikävektori aineistosta
  age2<-df$age-37.5;
  
  fit<-glm(smokes_imp ~ age2,family=binomial)
  d<-coef(fit)
  require(mvtnorm)
  out<-rmvnorm(n=1,mean=d,sigma=vcov(fit))
  dens<-dmvnorm(out,mean=d,sigma=vcov(fit),log=T)
  colnames(out)<-c("beta0","beta1")
  
  return(list("out"=out,"dens"=dens,"d"=d,"vcov"=vcov(fit)))
  # out: simuloidut estimaatit (pituus 2)
  # dens: tiheysfunktion arvo
  # d: parametriestimaatit (pituus 2)
  # vcov: kovarianssimatriisin estimaatti (2x2)
}

theta_dens<-function(theta,d,vcov) {
  require(mvtnorm)
  dmvnorm(theta,mean=d,sigma=vcov,log=T)
}
# uutta:
# k-funktion log-uskottavuus
d_k_log<-function(smokes,theta_par,age,df) { #df on siksi että saadaan indeksi sille, mitkä ovat puuttuvia
  # smokes: tupakointitiedon vektori, johon on jo sijoitettu imputoinnit
  sel<-is.na(df$smokes_hav)
  # posterior for smokes_hav
  tmp<-cbind(1,(age-37.5))
  p_tmp<-expit(tmp %*% theta_par)
  sum(log(c(p_tmp[smokes==1],1-p_tmp[smokes==0])[sel]))
}

# g-funktio eli full data log-uskottavuus
d_g_log<-function(smokes,theta_par,age) { #theta_df on tällä funktiolla aina iteraatiolla t, eikä t-1
  # smokes: tupakointitiedon vektori, johon on jo sijoitettu imputoinnit
  # posterior for smokes_hav
  tmp<-cbind(1,(age-37.5))
  p_tmp<-expit(tmp %*% theta_par)
  sum(log(c(p_tmp[smokes==1],1-p_tmp[smokes==0])))
}
#--

theta_dens(c(0,0.1),mu,M)
#system.time({
  y <- popmc(df_=df,par_=list(nsim = 1000,x=c(0.0,0.1), mu = mu, sdmat = M),theta_dens_=theta_dens,theta_prop_=theta_prop,d_k_log_=d_k_log,d_g_log_=d_g_log)
#})
names(y) #$smokes_hav
summary(y$smokes_hav)
summary(df)
