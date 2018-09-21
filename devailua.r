# testataan miten saadaan R:ssä kirjoittu funktio kutsutta C++:sta
# Juho Kopra / 20.9.2018
library(Rcpp)
library(inline)
library(mvtnorm)

theta_dens<-function(theta,d,vcov) {
  require(mvtnorm)
  dmvnorm(theta,mean=d,sigma=vcov,log=T)
}

src <- '
  RNGScope scp;
  Rcpp::Function td("theta_dens");
  return td(theta,d,vcov);
'
fun <- cxxfunction(signature(theta_dens="function",theta="double",d="double",vcov="matrix"), src, plugin="Rcpp")
set.seed(42)
fun(theta_dens,c(-0.6,1.07),c(0.12,0.85),vcov=matrix(c(1,0.12,0.12,1),2,2))
theta_dens(c(-0.6,1.07),c(0.12,0.85),vcov=matrix(c(1,0.12,0.12,1),2,2))
