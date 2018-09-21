
popmc <- function(df_,par_,theta_dens_,theta_prop_,d_k_log_,d_g_log_) {
#	mis<-apply(df_,2,function(x) any(is.na(x)))
#	newname<-paste(names(mis),"_mis",sep="")[mis]
#	newvalues<-as.integer(is.na(df_[,mis]))
#	df_[,dim(df_)[2]+1]<-newvalues
#	names(df_)<-c(names(df_),newname)
	y <- .Call("popmc",df_=df_,ymis_index_=which(is.na(df_$smokes_hav)),parlist_=par_,theta_dens_=theta_dens_,theta_prop_=theta_prop_,d_k_log=d_k_log_,d_g_log=d_g_log_,PACKAGE="PopMC")
#	for(i in 1:(length(y)-1)) {
#		y[[i]] <- y[[i]][,-1]
#	}
	y
}