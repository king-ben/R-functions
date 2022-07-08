#Returns the height of a taxon
#requires phytools library
tip.height <- function(tree, tip){return(max(nodeHeights(tree))-nodeheight(tree, which(tree$tip.label == tip)))}

tipHeights <- function(tree){
  max <- max(nodeHeights(tree))
  th <- vector()
  for(i in 1:length(tree$tip.label)){
    th[i] <- nodeheight(tree, i)
  }
  return(max-th)
}

node.height <- function(node, tree){
  max(nodeHeights(tree)) - nodeheight(tree, node)
}

nodeheights <- function(tree){
  edg <- match((tree$Nnode+2):(2*tree$Nnode+1), tree$edge[,1])
  pt <- nodeHeights(tree)
  ht <- pt[edg,1]
  ht <- max(pt)-ht
}

