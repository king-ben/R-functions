library("caper")
monophyl.multi <- function(trees, clades, exclude.rogue=F, rogue.taxa, burnin=0.1){
  0 -> count
  trees[((round(length(trees)*burnin))+1):length(trees)] -> treesb
  if(exclude.rogue == TRUE){
    list() -> treesc
    for(i in 1:length(treesb)){
      drop.tip(treesb[[i]], rogue.taxa) -> treesc[[i]]
    }
  }
  if(exclude.rogue == FALSE)
    treesb -> treesc
  vector() -> fails
  vector() -> passes
  for(i in 1:length(treesc)){
    1 -> state
    
    for(j in 1:length(clades)){
      getMRCA(treesc[[i]], clades[[j]]) -> node
      if(length(clade.members(node, treesc[[i]])) > length(clades[[j]]))
        0 -> state
    }
    if(state ==1){
      count + 1 -> count
      append(passes, as.numeric(i)) -> passes
    }
    if(state ==0)
      append(fails, as.numeric(i)) -> fails
  }
  list() -> results
  count/length(treesc)-> results$posterior
  treesb[passes] -> results$true.trees
  treesb[fails] -> results$false.trees
  passes + (length(trees) - length(treesb)) -> results$truegens
  fails + (length(trees) - length(treesb)) -> results$falsegens
  return(results)
}

