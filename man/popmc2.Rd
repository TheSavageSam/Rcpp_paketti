\name{PopMCarmatest}
\alias{PopMCarmatest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PopMC2
}
\description{
Calculates a Population Monte Carlo algorithm.
}
\usage{
popmc2(df_, par_, theta_dens_, theta_prop_, d_k_log_, d_g_log_)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_}{
df
}
  \item{par_}{
par
}
  \item{theta_dens_}{
theta_dens_ function
}
  \item{theta_prop_}{
theta_prop funtiocn
}
  \item{d_k_log_}{
fucntion
}
  \item{d_g_log_}{
function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (df_, par_, theta_dens_, theta_prop_, d_k_log_, d_g_log_) 
{
    .Call(`_popmc2`, df_, par_, theta_dens_, theta_prop_, 
        d_k_log_, d_g_log_)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
