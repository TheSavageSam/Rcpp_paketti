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
txt2 <- '#include <Rcpp.h>
#include "cppSample.h"
#include "rinitials.h"

RcppExport SEXP popmc(SEXP df_, SEXP parlist_) {

try {          // or use BEGIN_RCPP macro
//take parameters in
Rcpp::DataFrame df = Rcpp::DataFrame(df_);
Rcpp::List parlist = Rcpp::List(parlist_);

//parameters
int nsim = Rcpp::as<int>(parlist["nsim"]);
double sigmax_init = Rcpp::as<double>(parlist["sigmax_init"]);
double initial_value = Rcpp::as<double>(parlist["init_sim"]);
double home_eff_init = Rcpp::as<double>(parlist["home_eff_init"]);
double scaler_par = 0.95; // Parametrien arvoja pienennetään tällä kertoimella jokaista päivää kohti
// Käytä: double R_pow(double x, double y); => ::R_pow(x,y);

//data
std::vector<int> goals = df["goals"];
std::vector<double> day = df["day"];
std::vector<double> home = df["home"];

// lasketaan montako päivää kahden ottelun välillä on:
std::vector<double> daysdifference(day.size());
std::vector<double> scale(day.size());
std::vector<double> rate(day.size());
daysdifference[0] = 1;
scale[0] = scaler_par;
rate[0] = sigmax_init * scale[0];
for(int i = 1; i < day.size(); i++) {
daysdifference[i] = day[i]-day[i-1];
scale[i] = ::R_pow(scaler_par,daysdifference[i]);
rate[i] = sigmax_init * scale[i];
}

int ncol = goals.size()+1;
int nrow = nsim;
double tarkkuus = 1/sigmax_init;
Rcpp::NumericMatrix zline_mat(nrow,ncol);
Rcpp::NumericMatrix w_tilde_log(nrow,ncol);
Rcpp::NumericMatrix w_tilde(nrow,ncol);
Rcpp::NumericMatrix w_line_mat(nrow,ncol);
Rcpp::NumericMatrix z_mat(nrow,ncol);

std::vector<double> w_sum(nrow);
//Asetetaan nollia täyteen: Tarviiko tehdä erikseen vai onko valmiiksi?

for(int i = 0; i < nrow; i++) {
double initvalues = rinitials(initial_value*rate[0], 1/rate[0]);
zline_mat(i,0) = initvalues;
z_mat(i,0) = initvalues; //Tämä on sämplättyjen matriisi, mutta käytetään silti tätä näin syntaktisista syistä.

double tmp_w = ::Rf_dgamma(initvalues,3.0*rate[0], 1/rate[0], true);
w_tilde_log(i,0) = tmp_w;
w_tilde(i,0) = 1;
w_sum[i] = 0;
}


for(int j = 1; j < ncol; j++) {
for(int i = 0; i < nrow; i++) {
double tmp = ::Rf_rgamma(0.5*(z_mat(i,j-1)+goals[j-1])*rate[j-1], 1/rate[j-1]);
double proposal = tmp; //Rcpp::as<double>(tmp);
zline_mat(i,j) = proposal;

// probability weights of proposal distribution
double loglik_proposal = ::Rf_dgamma(tmp,0.5*(z_mat(i,j-1)+goals[j-1])*rate[j-1], 1/rate[j-1], true);

// TODO: probability weights of parameters change in time -distribution
double loglik_model = ::Rf_dgamma(proposal, z_mat(i,j-1)*rate[j-1], 1/rate[j-1], 1);

double heff = 1.0;
if (home[j-1] == 1) {
heff = 1.11;
} else {
heff = 2 - 1.11;
}
// probability weights of p(observations | parameters) -distribution (poisson)
double loglik_obsmodel = ::Rf_dpois(goals[j-1], heff * proposal, 1); // lisää data tähän goals[j-1]

w_tilde_log(i,j) = loglik_model + loglik_obsmodel - loglik_proposal;
w_sum[i] = w_tilde_log(i,j); // TODO: Viimeksi tässä ei ollut summaa?
w_tilde(i,j) = ::exp(w_tilde_log(i,j));

}

std::vector<double> w_tilde_exp(nrow);
std::transform(w_sum.begin(),w_sum.end(),w_tilde_exp.begin(),exp); // lasketaan log-painoista tavan painot
std::vector<double> w_line(nrow);
// Suoritetaan normeeraus käyttäen STL-funktiota. Normeerattu todennäköisyys on nyt vektori z.
double summa = std::accumulate(w_tilde_exp.begin(),w_tilde_exp.end(),0.0);
//if (summa > 0) {
std::transform(w_tilde_exp.begin(), w_tilde_exp.end(), w_line.begin(), std::bind1st(std::multiplies<double>(),1/summa));
//} else {
//  z[0] <- 1; // TODO: Toimiikohan?
//}
//TODO: tarkasta, että summa ei ole 0

//sämpläys
std::vector<double> sample_from(nrow);  // = Rcpp::as<std::vector<double> >(zline_mat(Rcpp::_,j));
for(int i = 0; i < nrow; i++) {
sample_from[i] = zline_mat(i,j);
w_line_mat(i,j) = w_line[i];
}
std::vector<double> sampled_values = cppSample(nrow, sample_from, w_line);
for(int i = 0; i < nrow; i++) {
z_mat(i,j) = sampled_values[i];
}
if (j%5==0) {
  Rprintf("\\v%d ajanhetkeä simuloitu...\\n",j);
  R_FlushConsole();
  R_ProcessEvents();
}

}

return(Rcpp::List::create(Rcpp::Named("Zline.mat") = zline_mat,
Rcpp::Named("log(Wtilde.mat)") = w_tilde_log,
Rcpp::Named("Wtilde.mat") = w_tilde,
Rcpp::Named("Wline.mat") = w_line_mat,
Rcpp::Named("Z.mat") = z_mat));

} catch( std::exception &ex ) {    // or use END_RCPP macro
forward_exception_to_r( ex );
} catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}'

rfunc_txt <- '
popmc <- function(df_,par_) {
y <- .Call("popmc",df_=df_,parlist_=par_,PACKAGE="PopMC")
for(i in 1:length(y)) {
y[[i]] <- y[[i]][,-1]
}
y
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
cat(rinitials_txt,file=paste("./",package_name,"/src/rinitials.cpp",sep=""))
cat(rinitialsh_txt,file=paste("./",package_name,"/src/rinitials.h",sep=""))
cat(txt2,file=paste("./",package_name,"/src/popmc.cpp",sep=""))
cat(rfunc_txt,file=paste("./",package_name,"/R/popmc.R",sep=""))
cat(mavfunc_txt,file=paste("./",package_name,"/R/plot_pf_results.R",sep=""))
cat(plotfunc_txt,file=paste("./",package_name,"/R/mean_and_var.R",sep=""))

compileAttributes(pkgdir = "./PopMC/",verbose=TRUE)

## compile in between
# cd /home/juho/Asiakirjat/UEFhommat/Tutkimusprojekti/Rcpp_paketti/
# sudo R CMD build PopMC
# sudo R CMD INSTALL PopMC

library(PopMC)

df <- data.frame(goals=c(3,3,5,3,2,0,1,0,3,1),day=round(cumsum(rnorm(10,3.5,sd=0.8))), home=rep(c(0,1),5))
df$htid <- rep(1:3,3,1)
?expand.grid
tmp.grid<-expand.grid(1:4,1:4)
pick.id <- c(1:12)[tmp.grid[,1] != tmp.grid[,2]]

#generate two games:
genRound <- function(teamNames,startDate) {
  nteams <- length(teamNames)
  teamids <- 1:nteams
  samp <- sample(x=teamids,size=nteams,replace=FALSE)
  names(samp) <- paste(c("home","away"),sort(rep(1:(nteams/2),2)),sep="") #TODO: Tee yleisempi versio käyttäen pastea
  pick.home <- seq(1,nteams,2)
  pick.away <- seq(2,nteams,2)
  data.frame("htid" = samp[pick.home],"atid" = samp[pick.away],hometeam = teamNames[samp[pick.home]],awayteam=teamNames[samp[pick.away]],
             homegoals = round(rnorm(nteams/2,2.87,1.2)), awaygoals = round(rnorm(nteams/2,2.35,1.2)), day = as.Date(startDate,"%d-%m-%Y")+round(rnorm(nteams/2,2.4,sd=0.8))
  )
}

set.seed(42)
genSeason <- function(teamNames,seasonStartDate,nrounds) {
  df <-genRound(teamNames=teamNames,startDate=seasonStartDate) #init
  for(i in 2:nrounds) {
    df <- rbind(df,genRound(teamNames=teamNames,startDate=max(df$day[c(length(df$day)-1,length(df$day))])))
  }
  rownames(df) <- 1:dim(df)[1]
  df
}

df <- genSeason(teamNames=c("JYP","Ässät","KalPa","HIFK"),seasonStartDate="02-09-2013",nrounds=5)
# Season generated

#rownames(df) <- paste(format(as.Date("16.9.2012","%d.%m.%Y")+df$day,"%d.%m."),c("Ässät-JYP","JYP-Ässät","Ässät-JYP","JYP-Ässät","Ässät-JYP","JYP-Ässät","Ässät-JYP","JYP-Ässät","Ässät-JYP","JYP-Ässät"),sep=" ")

#nopeustesti
library(PopMC)

nsim <- 60
df <- data.frame(goals=rpois(nsim,3.0),day=round(cumsum(rnorm(nsim,3.5,sd=0.8))), home=rep(c(0,1),nsim/2))
system.time({
  y <- popmc(df_=df,par=list(nsim = 1000, init_sim = 1, sigmax_init = 1, home_eff_init = 1.11))
})
##with covariates:
#df$cov <- rep(1,length(df[,1]))
#system.time({
#  y <- popmccov(df_=df,par=list(nsim = 1000, init_sim = 1, sigmax_init = 1, home_eff_init = 1.11))
#})

mav <- mean_and_var(y$Wline.mat,y$Zline.mat)
plot_pf_results(df=df,mav,c(0.025,0.5,0.975),main="JYPin maaliodotusarvo 2012")

#TODO: Tee predict-funtkio, jolle annetaan eri parametrissa otteluit joita halutaan ennustaa.
# predict-funktio ennustaa maaliodotusarvoa tuleviin otteluihin.
pf.predict <- function(mav,pred_data,par) {
  # a = moa[i-1] * sigma_init * 0.95^daysdiff[i], b.rate = sigma_init * 0.95^daysdiff[i]
  
  # 1. Laske parametrit skaalattuna s.e. niissä huomioidaan sekä koti/vierasetu että vastustajan maaliodotusarvo.
}
pf.predict(mav=mav,pred_data=df_pred,par=list(sigmax_init = 5, home_eff_init = 1.11))

