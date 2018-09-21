
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
}