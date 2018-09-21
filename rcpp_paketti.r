# Kehitetään Rcpp:n pohjalta R-paketti, jolla voi tehdä Population Monte Carlo -mallin sovitusta
# Juho Kopra
# 10.9.2018

rm(list=ls())
setwd("/home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/")
library(Rcpp)

Rcpp.package.skeleton(name="PopMC",author="Juho Kopra",email="juho.kopra@uef.fi")

rkoodi<-'
# Alustus
M<-50
Tt<-10
# d_g_log
# d_k_log
# theta_dens
# theta_prop
# k


theta_df<-matrix(c(-0.8,-0.1)+rnorm(2*M,sd=0.000001),ncol=2,nrow=M,byrow=T)
res_theta_array<-array(dim=c(Tt,dim(theta_df)))
imputed_data<-matrix(NA,ncol=M,nrow=length(df$smokes_hav))
theta_all<-matrix(NA,nrow=M,ncol=2)
dens<-vector("double",M)
r<-matrix(NA,nrow=Tt,ncol=M)
w<-matrix(NA,nrow=Tt,ncol=M)
mu_mat<-matrix(NA,nrow=M,ncol=2)
vcov_arr<-array(NA,dim=c(M,2,2))

nn<-matrix(NA,nrow=M,ncol=M)
dd<-matrix(NA,nrow=M,ncol=M)


aika<-system.time({
for(t in 1:Tt) {
# Askel (a) (M kpl toistoja)
for(i in 1:M) {
imputed_data[,i]<-k(i=i,df=df,theta_df=theta_df)
tmp<-theta_prop(1,smokes_imp=imputed_data[,i],df)
theta_all[i,1:2]<-tmp[[1]]
dens[i]<-tmp[[2]]
mu_mat[i,]<-tmp[[3]]
vcov_arr[i,,]<-tmp[[4]]
}

# Askel (b)

#Tarvitaan MxM matriisit full data likelihoodista, ehdotusjakaumasta kdim
# dimensiot: imputointi, theta-parametri
for(l in 1:M) {
for(i in 1:M) {

#lasketaan pi-jakauman tiheys
log_dens_pi<-theta_dens(theta_all[i,1:2],mu_mat[l,],vcov_arr[l,,])
log_dens_g <-d_g_log(smokes=imputed_data[,l],theta_all[i,],df)
log_dens_k_ll <-d_k_log(smokes=imputed_data[,l],theta_df[i,],df) # tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta
log_dens_k_li <-d_k_log(smokes=imputed_data[,l],theta_df[l,],df) # tässä nimenomaan theta_{t-1}^(i) eli edelliskerran theta

nn[i,l]<-log_dens_g-log_dens_k_ll-log(M)
dd[i,l]<-log_dens_k_li+log_dens_pi-log_dens_k_ll-log(M)
}
print(paste("Iteraatio ",t,":",l," suoritettu.",sep=""))
}

unscaled<-rowMeans(exp(nn-mean(nn)))
#alkuperäiset luvut saadaan unscaled*exp(mean(nn))
jakaja<-rowMeans(exp(dd))
rr<-unscaled/jakaja
ww<-rr/sum(rr) 
w[t,]<-ww
# Askel (c)
sampled<-sample.int(M,M,replace=T,prob=ww)
theta_df<-theta_all[sampled,]
# tallennetaan theta:t taulukkoon
res_theta_array[t,,]<-theta_all[sampled,]

print(paste("ITERAATIO ",t," suoritettu.",sep=""))

}'

