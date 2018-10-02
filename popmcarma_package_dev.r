rm(list=ls())
library(RcppArmadillo)
package_name <- "popmcarmatest"
setwd("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/")


unlink(paste("./",package_name,"/",sep=""),recursive=TRUE)

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
popmc2 <- function(df_,par_,theta_dens_,theta_prop_,d_k_log_,d_g_log_) {
  .Call(`_popmcarmatest_popmc2`, df_,par_,theta_dens_,theta_prop_,d_k_log_,d_g_log_)
}

#?RcppArmadillo.package.skeleton
#"theta_prop","theta_dens","k","d_k_log","d_g_log"
#RcppArmadillo.package.skeleton(name="hello_world")
RcppArmadillo.package.skeleton(name=package_name,example_code = F)
?RcppArmadillo.package.skeleton

unlink(paste("./",package_name,"/src/rcpparma_hello_world.cpp",sep=""))
unlink(paste("./",package_name,"/src/rcpparma_hello_world.h",sep=""))
#unlink(paste("./",package_name,"/R/rcpparma_hello_world.R",sep=""))
unlink(paste("./",package_name,"/man/rcpparma_hello_world.Rd",sep=""))

#file.copy(from, to, overwrite = recursive, recursive = FALSE,
#          copy.mode = TRUE, copy.date = FALSE)
filut<-c("apply_mean.cpp","popmc2.cpp","rk.cpp","cppSample.cpp","rcppSample.cpp","k.cpp","rinitials.cpp","apply_mean.h","cppSample.h","k.h","rcppSample.h","rinitials.h","rk.h") #,"popmc.cpp"
manut<-c("popmc2.Rd",paste(package_name,"-package.Rd",sep=""))
file.copy(from=paste("dev",filut,sep="/"), paste(package_name,"src",sep="/"))
file.copy(from=paste("man",manut,sep="/"), paste(package_name,"man",sep="/"),overwrite=T)

Rcpp::compileAttributes(pkgdir = paste("./",package_name,"/",collapse="",sep=""),verbose=TRUE)

## compile in between
# cd /home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/
# sudo R CMD build popmcarmatest && sudo R CMD INSTALL popmcarmatest



