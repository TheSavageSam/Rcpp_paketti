# Kehitellään R-koodia Laskennallisen tilastotieteen seminaaria varten.
# Population Monte Carlo aka. Iterated Importance Sampling (Celeux et al., 2006, Computational statistics and data analysis)
# 25.10.2014 Juho Kopra
# koodi2: kehitetään koodi s.e. se on selkeämpi, tehdään vain yksi imputointi aluksi (tajutaan koodia)
# koodi3: Rao-Blackwellized Iterated Importance Sampling
rm(list=ls())

# Tarvitaan data, jossa on puuttuvaa (simuloidaan sellainen)
expit<-function(x) { exp(x)/(1+exp(x)); }

# montako hlöä
N<-3000

# Ikä
age<-runif(n=N,25,50)
hist(age,100)

# Tupakointi
mu<-log(0.3/(1-0.3))
age_coef<- -0.125
p<-expit(mu+age_coef*(age-37.5))
smokes<-rbinom(n=N,size=1,prob=p)
plot(table(round(age),smokes))
plot(round(age),p,ylim=c(0,1))

# puuttuminen: satunnaista (MAR), riippuu vain iästä (ikä on havaittu kaikille)
mu_mis<-log(0.94/(1-0.94))
mis_coef<-0.09
rand_disturb<-rnorm(N,sd=0.1)
p_mis<-1-expit(mu_mis+mis_coef*(age-37.5)+rand_disturb) # p(s_i = NA|age,x)
plot(round(age),p_mis,ylim=c(0,1))
# simuloidaan osallistuminen
m_part<-rbinom(n=N,size=1,prob=1-p_mis) # 1 = osallistuminen
summary(m_part)
smokes_hav<-smokes
smokes_hav[m_part==0]<-NA
summary(smokes_hav)
plot(table(round(age),smokes_hav))

df<-data.frame(smokes_hav=smokes_hav,age=age,x=rand_disturb,osal=m_part)
head(df,10)
save.image(file="rao_blackwellized_PMC_28102014_T10_M500_data_only.Rdata")
load(file="rao_blackwellized_PMC_28102014_T10_M500_data_only.Rdata")
#load(file="lask_semma_data_05112014.Rdata")

# Koetetaan ratkaista ongelma Tärkeysotannalla
M<-sum(is.na(df$smokes_hav))
M
# Askel 0: Valitse alkuarvot parametreille

# Askel 1-T:
# (a) Imputoi puuttuvat arvot jakaumasta k(z|y,theta) ja pi(theta|y,z_t)
k<-function(i,df,theta_df) {
  smokes_imp<-df$smokes_hav
  # posterior for smokes_hav
  sel<-is.na(df$smokes_hav)
  tmp<-cbind(1,(df$age-37.5))
  p_tmp<-expit(tmp %*% theta_df[i,]) # p_tmp on [Nhav,M]-matriisi (jokaiselle i = 1,...,M lasketaan imputointi)
  smokes_imp[sel]<-rbinom(n=sum(sel),size=1,prob=p_tmp) #TODO: Tarvitaan M imputointia
  smokes_imp
}


# tarvitaan IS:iä myös ihan theta-parametrien eli regressiokertoimien estimointiin
theta_prop<-function(i,smokes_imp,df) {
  age2<-df$age-37.5;
  
  fit<-glm(smokes_imp ~ age2,family=binomial)
  d<-coef(fit)
  require(mvtnorm)
  out<-rmvnorm(n=1,mean=d,sigma=vcov(fit))
  dens<-dmvnorm(out,mean=d,sigma=vcov(fit),log=T)
  colnames(out)<-c("beta0","beta1")
  
  return(list("out"=out,"dens"=dens,"d"=d,"vcov"=vcov(fit)))
}
theta_dens<-function(theta,d,vcov) {
  dmvnorm(theta,mean=d,sigma=vcov,log=T)
}

# k-funktion log-uskottavuus
d_k_log<-function(smokes,theta_df,df) {
  # smokes: tupakointitiedon vektori, johon on jo sijoitettu imputoinnit
  sel<-is.na(df$smokes_hav)
  # posterior for smokes_hav
  df$age[sel]
  tmp<-cbind(1,(df$age-37.5))
  p_tmp<-expit(tmp %*% theta_df)
  sum(log(c(p_tmp[smokes==1],1-p_tmp[smokes==0])[sel]))
}

# g-funktio eli full data log-uskottavuus
d_g_log<-function(smokes,theta_df,df) { #theta_df on tällä funktiolla aina iteraatiolla t, eikä t-1
  # smokes: tupakointitiedon vektori, johon on jo sijoitettu imputoinnit
  sel<-is.na(df$smokes_hav)
  # posterior for smokes_hav
  df$age[sel]
  tmp<-cbind(1,(df$age-37.5))
  p_tmp<-expit(tmp %*% theta_df)
  sum(log(c(p_tmp[smokes==1],1-p_tmp[smokes==0])))
}


#------------------------------------------------
# Kehitetään varsinainen algoritmi (ilman 1. vaiheen resämpläystä)

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

}
})
aika/3600
save.image(file="rao_blackwellized_PMC_28102014_T10_M500.Rdata") # Ajo kesti 4h!

load(file="rao_blackwellized_PMC_28102014_T10_M500.Rdata")

# Lasketaan asymptoottisesti harhaton estimaattori
# 1. momentti
first_mom<-c(0,0)
dim(res_theta_array)
dim(w)
first_mom[1]<-sum(res_theta_array[,,1]*w)/Tt
first_mom[2]<-sum(res_theta_array[,,2]*w)/Tt

# 2. momentti
second_mom<-c(0,0)
second_mom[1]<-sum((res_theta_array[,,1]^2)*w)/Tt
second_mom[2]<-sum((res_theta_array[,,2]^2)*w)/Tt
second_mom

variance<-second_mom-first_mom^2
std<-sqrt(variance)

# odotusarvo
first_mom
# keskihajonta
std
# luottoväli/luottamusväli
first_mom-1.96*std
first_mom+1.96*std
# oikeat luvut
c(mu,age_coef)

pick<-sample.int(Tt*M, prob=w,replace=T)
hist(as.vector(res_theta_array[,,1])[pick],100)
abline(v=mu,lwd=3)
hist(as.vector(res_theta_array[,,2])[pick],100)
abline(v=age_coef,lwd=3)
save.image(file="pmc_rao_blackwell_30min_run_05112014.Rdata")
