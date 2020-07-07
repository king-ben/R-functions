#input mean and standard deviation in real space, outputs mean and sd in log space
convert.ln <- function(mean, sd){
	log(mean)−0.5*log(((sd/mean)^2)+1) -> lmean
	sqrt(log(((sd/mean)^2)+1)) -> lsd
	return(c(lmean, lsd))
}

