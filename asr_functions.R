#ALL FUNCTIONS ASSUME FIRST SITE IN EACH PARTITION IS AN ASCERTAINMENT COLUMN

# extracts partitions from a nexus file
# reads the assumptions block and converts it to a 3-cloumn dataframe
# columns are partition name, from index and to index
find_partitions <- function(file){
  nex <- readLines(file)
  from <- grep("begin assumptions", nex, ignore.case=T)+1
  from <- from[grep("charset", nex[from])]
  if(length(from)==0){stop("assumptions block with character sets not found")}
  to <- which(nex=="end;")
  to <- to[to>from][1]-1
  cb <- nex[from:to]
  cb <- gsub(".*charset ", "", cb)
  cb <- gsub(";", "", cb)
  cb <- gsub(" = ", "\t", cb)
  cb <- gsub("=", "\t", cb)
  cb <- gsub("-", "\t", cb)
  cbd <- sapply(cb, strsplit, split="\t")
  cbd <- t(as.data.frame(cbd))
  rownames(cbd) <- NULL
  colnames(cbd) <- c("concept", "from", "to")
  cbd <- transform(cbd, from=as.numeric(from), to=as.numeric(to))
  return(cbd)
}

#makes a vector of character state labels read from nexus characterstatelabels block
find_site_names <- function(file){
  nex <- readLines(file)
  if(length(grep("begin charstatelabels", nex, ignore.case = T))>0){
    from <- grep("begin charstatelabels", nex, ignore.case=T)+1
    if(length(from)==0){stop("characterstatelabels block not found")}
    to <- which(nex=="end;")
    to <- to[to>from][1]-1
    cb <- nex[from:to]
    cb <- gsub("    ", "", cb)
    cb <- gsub(".*\\s", "", cb)
    return(cb)
  }else{
    from <- grep("charstatelabels", nex, ignore.case=T)+1
    if(length(from)==0){stop("characterstatelabels block not found")}
    to <- which(nex==";")
    to <- to[to>from][1]-1
    cb <- nex[from:to]
    cb <- gsub("    ", "", cb)
    cb <- gsub(".*\\s", "", cb)
    cb <- gsub(",", "", cb)
    return(cb)
  }
}


# Makes a table of rate parameter links
# binning is by default 1-10, 11-20 (excl ascertainment)
rate_links <- function(wp, cutoff=10){
  n <- wp[,3]-wp[,2]
  wp <- wp[order(n),]
  n <- n[order(n)]
  names(n) <- wp[,1]
  b <- wp[,1:2]
  for(i in 1:ceiling(max(n)/cutoff)){
    b[which(n>(cutoff*(i-1)) & n<(cutoff*i+1)),2] <- names(n)[which(n>(cutoff*(i-1)) & n<(cutoff*i+1))][1]
  }
  names(b)[2] <- "rate"
  return(b)
}

order_partitions <- function(wp){
  n <- wp[,3]-wp[,2]
  wp <- wp[order(n),]
  return(wp)
}

#makes logger of likelihood of every cognate set
siteprobs_blocks <- function(parts, names, logevery=1000){
  logger <- newXMLNode("logger")
  xmlAttrs(logger) <- c(id="Sitelikelihoods", fileName="site_likelihoods.txt", logEvery=logevery, mode="tree")
  for(i in 1:nrow(parts)){
    log <- newXMLNode("log", parent=logger)
    pnam <- parts[i,1]
    first <- parts[i, 2]
    last <- parts[i, 3]
    cnam <- names[first:last, 2]
    
    xmlAttrs(log) <- c(id=paste("sitelik", pnam, sep="."), spec="babel.util.SiteLikelihoodLogger", likelihood=paste("@treeLikelihood", pnam, sep="."), value=paste(cnam, collapse=" "))
  }
  return(logger)
}

# makes ancestral state reconstruction logger
asr_blocks <- function(parts, names=NULL, links, id="AncestralSequenceLogger", logevery=1000, taxonset, logOrigin=FALSE, fileName="asr_logger.txt"){
  logger <- newXMLNode("logger")
  xmlAttrs(logger) <- c(id=id, fileName=fileName, logEvery=logevery, mode="tree")
  for(i in 1:nrow(parts)){
    log <- newXMLNode("log", parent=logger)
    pnam <- parts[i,1]
    link <- links[which(links[,1]==pnam),2]
    first <- parts[i, 2]
    last <- parts[i, 3]
    xmlAttrs(log) <- c(id=paste(id, pnam, sep="."), spec="beast.evolution.likelihood.AncestralStateLogger", data=paste("@orgdata", pnam, sep="."), siteModel=paste("@SiteModel.s:", link, sep=""), branchRateModel="@RelaxedClock.c:clock", tree="@Tree.t:tree", taxonset=paste("@", taxonset, sep=""))
    if(is.null(names)){
      xmlAttrs(log) <- c(value=paste(paste0(pnam, ".ascertainment"), paste0(pnam, ".", 1:(last-first), collapse=" ")))
    } else{
      cnam <- names[first:last]
      xmlAttrs(log) <- c(id=paste(id, pnam, sep="."), spec="beast.evolution.likelihood.AncestralStateLogger", data=paste("@orgdata", pnam, sep="."), siteModel=paste("@SiteModel.s:", link, sep=""), branchRateModel="@RelaxedClock.c:clock", tree="@Tree.t:tree", value=paste(cnam, collapse=" "), taxonset=paste("@", taxonset, sep=""))
    }
    if(logOrigin==TRUE){
      xmlAttrs(log) <- c(logParent="true", logMRCA="false")
    }
  }
  return(logger)
}

#calculates synonymy at ancestral state reconstruction node
# input is the logger output by the beast analysis after it has had collapse_covarion.py run on it
# asr logger is read in with fread
asr_synonymy <- function(asr, wp, burnin=0.1, thinfactor=10){
  b <- round(burnin*nrow(asr))
  n <- nrow(asr)
  nc <- ncol(asr)
  asr <- asr[seq(b, n, by=thinfactor)]
  n <- nrow(asr)
  syn <- matrix(NA, nrow = n, ncol = nrow(wp))
  for(i in 1:nrow(wp)){
    for(j in 1:n){
      t <- asr[j, wp[i,2]:wp[i,3]]
      syn[j,i] <- sum(t)
    }
  }
  return(syn)
}

#find changes on a branch. Takes asr from the nodes at either end of the branch
get_innovations <- function(asr_mrca, asr_parent, burnin=0.1, siglevel=0.5, mode="innovations"){
  b <- round(burnin*nrow(asr_mrca))
  n <- nrow(asr_mrca)
  asr_mrca <- asr_mrca[b:n,]
  asr_parent <- asr_parent[b:n,]
  
  if(mode=="innovations"){
    asr_comp <-ifelse(asr_parent==0 & asr_mrca==1, 1, 0)
  }
  if(mode=="losses"){
    asr_comp <-ifelse(asr_parent==1 & asr_mrca==0, 1, 0)
  }
  
  asr_mean <- apply(asr_comp, 2, mean)
  sig <- asr_mean[asr_mean>siglevel]
  sig <- sig[order(sig, decreasing=T)]
  return(sig)
}

# find innovations on terminal branches i.e. the data at the tips is known and fixed
# takes a vector of tip states
# MUST be in same order as asr file suggested code:
# wn <- find_site_names("../make_asr_xml/ie.nex")
# cognateorder <- match(names(asr_T), wn)
# taxon_data <- nex[[taxon]]
# taxon_data <- taxon_data[cognateorder]

get_innovations_terminal <- function(taxon_data, asr_parent, burnin=0.1, siglevel=0.5, mode="innovations"){
  b <- round(burnin*nrow(asr_parent))
  n <- nrow(asr_parent)
  asr_parent <- asr_parent[b:n,]
  asr_parent <- as.data.frame(asr_parent)
  
  if(mode=="innovations"){
    asr_parent_test <- asr_parent[,which(taxon_data==1)]
    asr_comp <- 1-asr_parent_test
  }
  if(mode=="losses"){
    asr_comp <- asr_parent[,which(taxon_data==0)]
  }
  
  asr_mean <- apply(asr_comp, 2, mean)
  sig <- asr_mean[asr_mean>siglevel]
  sig <- sig[order(sig, decreasing=T)]
  return(sig)
}

