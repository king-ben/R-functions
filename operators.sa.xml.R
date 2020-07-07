operators.xml.sa <- function(taxondates, outputfile, treeid="tree", weight=10, offset=FALSE, offsetid="offset_tree"){
	
	######### INITIALISE ########
	
	#Read in the table
	read.table(taxondates, header=T) -> taxa
	
	#####GENERATE XML ITEMS###########
	
	operator <- newXMLNode("operator")
	xmlAttrs(operator) = c(spec = "SampledNodeDateRandomWalker", windowSize="10", tree=paste("@", treeid, sep=""), weight=weight)
		#if using offset, add treeWOffset attribute
	if(offset==TRUE){
		xmlAttrs(operator) = c(treeWOffset=paste("@", offsetid, sep=""))
	}
	taxonset <- newXMLNode("taxonset", parent=operator)
	xmlAttrs(taxonset) = c(spec="TaxonSet")
	for(i in 1:nrow(taxa)){
	
		#For each taxon in the site...
		taxon <- newXMLNode("taxon", parent=taxonset)
		xmlAttrs(taxon) = c(idref=as.character(taxa$taxon[i]))
		
		dates <- newXMLNode("samplingDates", parent=operator)
		xmlAttrs(dates) = c(id=paste("samplingDate", i, sep=""), spec="beast.evolution.tree.SamplingDate", taxon=paste("@", taxa$taxon[i], sep=""), upper=taxa$upper[i], lower=taxa$lower[i])
		
	}
	
	###### GENERATE OPUTPUT XML BLOCK
	
	sink(file=outputfile)
	print(operator)
	sink()
}	










