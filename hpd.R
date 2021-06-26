hpd <- function(posterior, p=0.95){
	unname(sort(posterior)) -> sp
	round(length(sp)*p) -> n
	ints <- vector()
	for(i in 1:ceiling(length(sp)*(1-p))){
		append(ints, sp[i+(n-1)]-sp[i]) -> ints
		}
	c(sp[which(ints == min(ints))], sp[which(ints == min(ints))+(n-1)]) -> ci
	if(length(ci)>2){
		list() -> ci2
		ci[1:(length(ci)/2)] -> mins
		unique(mins) -> mins
		ci[(length(ci)/2+1):length(ci)] -> maxs
		unique(maxs) -> maxs
		for(i in 1:length(mins)){
			c(mins[i], maxs[i]) -> ci2[[i]]
		}
		return(ci2)
	}
	else{
		return(ci)
	}
}



#Here is the code used by BEAST/tracer

#public static double[] HPDInterval(double proportion, double[] x, int[] indices) {
#
#        double minRange = Double.MAX_VALUE;
#        int hpdIndex = 0;
#
#        final int diff = (int) Math.round(proportion * (double) x.length);
#        for (int i = 0; i <= (x.length - diff); i++) {
#           final double minValue = x[indices[i]];
#            final double maxValue = x[indices[i + diff - 1]];
#            final double range = Math.abs(maxValue - minValue);
#THIS IF STATEMENT MEANS THAT ONLY THE LOWEST POSSIBLE HPD IS TAKEN IF THERE IS MORE THAN ONE (ie becasue 
#            if (range < minRange) {
#                minRange = range;
#                hpdIndex = i;
#            }
#        }
#
#        return new double[]{x[indices[hpdIndex]], x[indices[hpdIndex + diff - 1]]};
#    }


