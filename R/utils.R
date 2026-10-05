is_legitimate_set <- function(set, graph){
  # the set is a sub-dag of the original DAG
  # the set is legitimate if all ancestors of each member of the set are also members
  ## of the set
  # set: set of pvalues to check, has to be named vertices
  # graph: original graph
  check <- sapply(set%>%names, function(x){
    p <- subcomponent(graph, x, mode="in") #set of all ancestors, including the node itself
    if(all(names(p) %in% names(set)))
      return(TRUE) else return(FALSE)
  }
  )
  
  if(all(check))
    return(list("is.legitimate"=TRUE, "check" = check)) else 
      return(list("is.legitimate" = FALSE, "check" = check))
}

is_self_consistent_set <- function(pvals, alpha,N){
  max_thresh <- length(pvals)*alpha/N
  if(max(pvals)<=max_thresh)
    return(TRUE) else return(FALSE)
}

BH.proc <- function(pvals, alpha,n=length(pvals)){
  # wrapper around p.adjust to apply the BH procedure
  # returns number of doscoveries and names of rejected p-values
  adjp <- p.adjust(pvals,method="BH",n)
  rej.pvals <- which(adjp<=alpha)
  return(list("N_Discoveries"=max(0, length(rej.pvals)), "Rejected_pvalues" = names(rej.pvals)))
}

BY.proc <- function(pvals, alpha, n=length(pvals)){
  # wrapper around p.adjust to apply the BY procedure
  # returns number of doscoveries and names of rejected p-values
  adjp <- p.adjust(pvals,method="BY",n)
  rej.pvals <- which(adjp<=alpha)
  return(list("N_Discoveries"=max(0, length(rej.pvals)), "Rejected_pvalues" = names(rej.pvals)))
}



# G is a list with sublists
## Pa - Parents
#convert a graph to focusedBH type
graph_as_list <- function(graph){
  m <- length(graph)
  Pa <- lapply(1:m, function(x) {y=find.parents(graph,x) %>% names
  return(which(V(graph)$name %in% y))})
  C <- lapply(1:m, function(x) {y=find.children(graph,x) %>% names
  return(which(V(graph)$name %in% y))})
  G <- list("m"=m, "Pa"=Pa, "C"=C)
  return(G)
}

Yekutieli.proc <- function(decisions, theta){
  #decisions: a vector of N T/F according as the corresponding p-value is rejected/not rejected by Yekutieli's method
  #theta:  a vector of true states (0-true null, 1 - false null) of the N hypotheses
  #returns fdp and power
  fdp <- sum(theta[decisions==TRUE]==0)/max(sum(decisions),1)
  pow <- sum(theta[decisions==TRUE]==1)/sum(theta)
  No_rejection <- sum(decisions)
  Rejections <- which(decisions)
  return(list("fdp"=fdp, "pow" = pow, "No_rejection"=No_rejection, "Rejections"=Rejections))
}
