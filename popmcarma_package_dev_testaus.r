## compile in between
# cd /home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/
# sudo R CMD build PopMC && sudo R CMD INSTALL PopMC
# sudo R CMD build popmcarmatest && sudo R CMD INSTALL popmcarmatest
rm(list=ls())
library(popmcarmatest)

setwd("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/testialue")
load(file="testidata.Rdata",verbose=T)

head(df,20)
names(df)
mu<-c(0.12,0.85)
M<-matrix(c(1,0.12,0.12,1),2,2)

expit <- function(x) { exp(x)/(1+exp(x)); }
theta_prop<-function(smokes_imp,age) {
  # smokes_imp: imputoitu tupakointivektori
  # age: ikävektori aineistosta
  age2<-age-37.5;
  
  fit<-glm(smokes_imp ~ age2,family=binomial)
  d<-coef(fit)
  require(mvtnorm)
  out<-rmvnorm(n=1,mean=d,sigma=vcov(fit))
  dens<-dmvnorm(out,mean=d,sigma=vcov(fit),log=T)
  colnames(out)<-c("beta0","beta1")
  
  return(list("out"=out,"dens"=dens,"d"=as.vector(d),"vcov"=as.matrix(vcov(fit))))
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
# (a) Imputoi puuttuvat arvot jakaumasta k(z|y,theta) ja pi(theta|y,z_t)
k<-function(theta_par,age,smokes_imp,which_missing) {
  # Tämä täytyy ajaa jokaisella imputaatiolla/parametrirealisaatiolla i.
  # posterior for smokes_hav
  sel<-which_missing
  tmp<-cbind(1,(age-37.5))
  p_tmp<-expit(tmp %*% theta_par) # p_tmp on [Nhav,M]-matriisi (jokaiselle i = 1,...,M lasketaan imputointi)
  smokes_imp[sel]<-rbinom(n=sum(sel),size=1,prob=p_tmp) #TODO: Tarvitaan M imputointia
  smokes_imp
}

# k-funktion log-uskottavuus
d_k_log<-function(smokes,theta_par,age,ymis_index) { #ymis_index on siksi että saadaan indeksi sille, mitkä ovat puuttuvia
  # smokes: tupakointitiedon vektori, johon on jo sijoitettu imputoinnit
  sel<-ymis_index
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
tmp=theta_prop(c(0,0,1),c(34,45,56))
class(tmp[[1]])
as.vector(tmp[[1]])
class(tmp[[3]])
tmp[[3]]
library(popmcarmatest)

head(df)
y <- popmc2(df_=df,ymis_index_=which(is.na(df$smokes_hav)),parlist_=list(nsim = 10, x=c(0.0,0.1), mu = mu, sdmat = as.matrix(M)),
              theta_dens=theta_dens,theta_prop=theta_prop,d_k_log=d_k_log,d_g_log=d_g_log)
y
1000*y$vcov
class(as.matrix(y$vcov))
#names(y) #$smokes_hav
#summary(y$smokes_hav)
#summary(df)

# gdb debuggaus:
# library(popmcarmatest)
# source("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/popmcarma_package_dev_testaus.r")
