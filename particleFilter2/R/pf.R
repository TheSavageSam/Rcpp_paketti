
pf <- function(df_,par_) {
y <- .Call("pf",df_=df_,parlist_=par_,PACKAGE="particleFilter2")
for(i in 1:length(y)) {
y[[i]] <- y[[i]][,-1]
}
y
}