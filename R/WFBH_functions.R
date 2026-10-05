
calculate.weights <- function(graph, pvalues, fbh = FALSE, lambda=0.5){
  #if fbh=TRUE, calculate weights for Focused Storey BH
  depths <- find.depth(graph)
  D <- depths$Max_depth
  nd <- n.d(graph)
  hd <- sapply(depths$Depths, length)
  weights <- vector("list", D)
  if(fbh==TRUE){
    pvals <- unlist(pvalues)
    wt <- (sum(pvals>= lambda) +1)/((1-lambda)*length(pvals))
    weights <- lapply(depths$Depths, function(x) rep((1/wt), length(x))%>% setNames(nm=names(x)) )
  } else {
    if(hd[1]>1)
      pi_hat1 <- (sum(pvalues[names(depths$Depths[[1]])] >lambda)+1)/((1-lambda)*hd[1]) else pi_hat1 <- 1
    weights[[1]] <-  rep((1/pi_hat1), times=length(depths$Depths[[1]]))%>% setNames(nm=names(depths$Depths[[1]]))
    
    if(D>1){
      for(d in 2:D){
      w_di_hat <- sapply(depths$Depths[[d]]%>%names, function(x){
      parents <- find.parents(graph, x)
      children_at_d <- lapply(parents, function(y){
        children <- find.children(graph, y)%>%names
        children[children %in% names(depths$Depths[[d]]) ]%>% return()
      })
      
      pi_hat <- lapply(children_at_d, function(xy) { 
        numerator <- sum((pvalues[xy] >= lambda)) + 1
        denominator <- (1-lambda)*length(xy)
        return(numerator/denominator)})
      
      w_i_hat <- mapply(function(x,y) {
          bestcase <- (1/(1-lambda))*nd[d]/hd[d]
          if(hd[d]==1)
            return(1) else 
            if(hd[d]>1 & bestcase>1){
              pi_hat0 <- (sum(pvalues[names(depths$Depths[[d]])] >lambda)+1)/((1-lambda)*hd[d])
            return(1/pi_hat0)} else{
              w <- hd[d]/(x*length(y)*nd[d])
              return(w)} 
      },
      pi_hat, children_at_d, SIMPLIFY=TRUE) #check if need to keep SIMPLIFY=F as before
      return(max(w_i_hat))
    })
    
    weights[[d]] <- w_di_hat
      }
    }
  }
  
  return(weights)
}

new.struct.bh.proc <- function(graph, weighted=FALSE, type, alpha=0.05, 
                               lambda=0.5, fbh=FALSE){
  # Implements the FBH and WFBH procedures with DAG structure filter
  # Inputs: graph; weighted, indicating whether weights are used; type, specifying the
  # smoothing method; FDR threshold alpha; lambda, used in the WFBH weights (default = 0.5);
  # and fbh, indicating whether FBH is implemented.
  n.disc=0; rej.pvals = NULL; stopcond=FALSE
  theta <- V(graph)$theta 
  pvals <- V(graph)$pvalues
  if (type == "simple")
    pvals <- V(graph)$pvalues else if (type == "fisher")
      pvals <- V(graph)$pvalues.fisher else if (type == "stouffer")
        pvals <- V(graph)$pvalues.stouffer else if (type == "simes")
          pvals <- V(graph)$pvalues.simes else if (type == "tippett")
            pvals <- V(graph)$pvalues.tippett
  
  names(pvals) <- names(theta) <- names(V(graph))
  
  if(fbh & weighted)
    keep_indices <- pvals <= lambda 
  
  if (weighted){
    w <- unlist(calculate.weights(graph, pvalues=pvals, fbh, lambda=lambda))
    pvals <- pvals/w[names(pvals)]
  }
  names(pvals) <- names(V(graph))
  
  if(fbh & weighted)
    pvals <- pvals[keep_indices]
  BH.results <- BH.proc(pvals, alpha=alpha, n=length(graph))
  rej.pvals <- V(graph)[BH.results$Rejected_pvalues]
  n.disc=length(rej.pvals)
    
  
  if(n.disc>0){
    
    legitimate.set <- is_legitimate_set(rej.pvals, graph)
    
    while(!(legitimate.set$is.legitimate)){
      
      R.minus <-  rej.pvals[which(legitimate.set$check==TRUE)] 
      if(length(R.minus)>0){
        self.consistent <- is_self_consistent_set(pvals[names(R.minus)], alpha, N = length(graph))
        if(self.consistent){
          n.disc=length(R.minus)
          rej.pvals=R.minus
          stopcond=TRUE
        }else{
          BH.results = BH.proc(pvals=pvals[names(R.minus)], 
                               alpha = (alpha * length(R.minus))/length(graph)) #length(graph) was length(pvals)
          rej.pvals <- R.minus[BH.results$Rejected_pvalues]
          n.disc=length(rej.pvals)
          if(n.disc==0){
            rej.pvals=NULL
            break} else
              legitimate.set <- is_legitimate_set(rej.pvals, graph)
        }
      } else {
        n.disc=0
        rej.pvals=NULL
        break
      }
      if(stopcond) break
    }
  } 
  false_rej <- sum(theta[rej.pvals] == 0)
  true_rej <- sum(theta[rej.pvals] == 1)
  fdp <- false_rej/max(n.disc,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = n.disc, 
              "Rejections" = rej.pvals%>% na.omit))
}

new.struct.by.proc <- function(graph, weighted=FALSE, type, alpha=0.05, lambda=0.5, fbh=FALSE){
  # Implements the reshaped FBH and WFBH procedures with DAG structure filter
  # Inputs: graph; weighted, indicating whether weights are used; type, specifying the
  # smoothing method; FDR threshold alpha; lambda, used in the WFBH weights (default = 0.5);
  # and fbh, indicating whether FBH is implemented.
  n.disc=0; rej.pvals = NULL; stopcond=FALSE
  theta <- V(graph)$theta 
  pvals <- V(graph)$pvalues
  if (type == "simple")
    pvals <- V(graph)$pvalues else if (type == "fisher")
      pvals <- V(graph)$pvalues.fisher else if (type == "stouffer")
        pvals <- V(graph)$pvalues.stouffer else if (type == "simes")
          pvals <- V(graph)$pvalues.simes else if (type == "tippett")
            pvals <- V(graph)$pvalues.tippett
  
  names(pvals) <- names(theta) <- names(V(graph))
  
  if(fbh & weighted)
    keep_indices <- pvals <= lambda  
  
  if (weighted){
    w <- unlist(calculate.weights(graph, pvalues=pvals,  fbh, lambda=lambda))
    pvals <- w[names(pvals)]*pvals
  }
  names(pvals) <- names(V(graph))
  
  if(fbh & weighted)
    pvals <- pvals[keep_indices]
  
  BY.results <- BY.proc(pvals, alpha=alpha, n=length(graph))
  rej.pvals <- V(graph)[BY.results$Rejected_pvalues]
  n.disc=length(rej.pvals)
  
 
  if(n.disc>0){
    
    legitimate.set <- is_legitimate_set(rej.pvals, graph)
    
    while(!(legitimate.set$is.legitimate)){
      
      R.minus <-  rej.pvals[which(legitimate.set$check==TRUE)] 
      if(length(R.minus)>0){
        self.consistent <- is_self_consistent_set(pvals[names(R.minus)], alpha, N = length(graph))
        if(self.consistent){
          n.disc=length(R.minus)
          rej.pvals=R.minus
          stopcond=TRUE
        }else{
          BY.results = BY.proc(pvals=pvals[names(R.minus)], 
                               alpha = (alpha * length(R.minus))/length(graph))
          rej.pvals <- R.minus[BY.results$Rejected_pvalues]
          n.disc=BY.results$N_Discoveries
          if(n.disc==0){
            rej.pvals=NULL
            break} else
              legitimate.set <- is_legitimate_set(rej.pvals, graph)
        }
      } else {
        n.disc=0
        rej.pvals=NULL
        break
      }
      if(stopcond) break
    }
  } 
  false_rej <- sum(theta[rej.pvals] == 0)
  true_rej <- sum(theta[rej.pvals] == 1)
  fdp <- false_rej/max(n.disc,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = n.disc, 
              "Rejections" = rej.pvals%>% na.omit))
}


is_legitimate_set_outernode <- function(set, graph){
  # the set is a sub-dag of the original DAG
  # the set is legitimate if all ancestors of each member of the set are also members
  ## of the set
  # set: set of pvalues to check, has to be named vertices
  # graph: original graph
  check <- sapply(set%>%names, function(x){
    p <- subcomponent(graph, x, mode="out") #set of all descendants, including the node itself
    if((names(p)[1] %in% names(set)) & (length(names(p[-1]) %in% names(set))==0))
      return(TRUE) else return(FALSE)
  }
  )
  
  if(all(check))
    return(list("is.legitimate"=TRUE, "check" = check)) else 
      return(list("is.legitimate" = FALSE, "check" = check))
}

outernode.struct.bh.proc <- function(graph, weighted=FALSE, type, alpha=0.05, 
                                     lambda=0.5, fbh=FALSE){
  # Implements the FBH and WFBH procedures with outer node filter
  # Inputs: graph; weighted, indicating whether weights are used; type, specifying the
  # smoothing method; FDR threshold alpha; lambda, used in the WFBH weights (default = 0.5);
  # and fbh, indicating whether FBH is implemented.
  n.disc=0; rej.pvals = NULL; stopcond=FALSE
  theta <- V(graph)$theta 
  pvals <- V(graph)$pvalues
  if (type == "simple")
    pvals <- V(graph)$pvalues else if (type == "fisher")
      pvals <- V(graph)$pvalues.fisher else if (type == "stouffer")
        pvals <- V(graph)$pvalues.stouffer else if (type == "simes")
          pvals <- V(graph)$pvalues.simes else if (type == "tippett")
            pvals <- V(graph)$pvalues.tippett
  
  names(pvals) <- names(theta) <- names(V(graph))
  
  if (weighted){
    w <- unlist(calculate.weights(graph, pvalues=pvals, fbh, lambda=lambda))
    pvals <- w[names(pvals)]*pvals
  }
  names(pvals) <- names(V(graph))
  
  BH.results <- BH.proc(pvals, alpha=alpha)
  rej.pvals <- V(graph)[BH.results$Rejected_pvalues]
  n.disc=BH.results$N_Discoveries
  
  if(n.disc>0){
    
    legitimate.set <- is_legitimate_set_outernode(rej.pvals, graph)
    
    while(!(legitimate.set$is.legitimate)){
      
      R.minus <-  rej.pvals[which(legitimate.set$check==TRUE)] 
      if(length(R.minus)>0){
        self.consistent <- is_self_consistent_set(pvals[R.minus], alpha, N = length(pvals))
        if(self.consistent){
          n.disc=length(R.minus)
          rej.pvals=R.minus
          stopcond=TRUE
        }else{
          BH.results = BH.proc(pvals=pvals[R.minus], 
                               alpha = (alpha * length(R.minus))/length(pvals))
          rej.pvals <- R.minus[BH.results$Rejected_pvalues]
          n.disc=BH.results$N_Discoveries
          if(n.disc==0){
            rej.pvals=NULL
            break} else
              legitimate.set <- is_legitimate_set_outernode(rej.pvals, graph)
        }
      } else {
        n.disc=0
        rej.pvals=NULL
        break
      }
      if(stopcond) break
    }
  } 
  false_rej <- sum(theta[rej.pvals] == 0)
  true_rej <- sum(theta[rej.pvals] == 1)
  fdp <- false_rej/max(n.disc,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = n.disc, 
              "Rejections" = rej.pvals%>% na.omit))
}


outernode.struct.by.proc <- function(graph, weighted=FALSE, type, alpha=0.05, 
                                     lambda=0.5, fbh=FALSE){
  # Implements the reshaped FBH and WFBH procedures with outer node filter
  # Inputs: graph; weighted, indicating whether weights are used; type, specifying the
  # smoothing method; FDR threshold alpha; lambda, used in the WFBH weights (default = 0.5);
  # and fbh, indicating whether FBH is implemented.
  n.disc=0; rej.pvals = NULL; stopcond=FALSE
  theta <- V(graph)$theta 
  pvals <- V(graph)$pvalues
  if (type == "simple")
    pvals <- V(graph)$pvalues else if (type == "fisher")
      pvals <- V(graph)$pvalues.fisher else if (type == "stouffer")
        pvals <- V(graph)$pvalues.stouffer else if (type == "simes")
          pvals <- V(graph)$pvalues.simes else if (type == "tippett")
            pvals <- V(graph)$pvalues.tippett
  
  names(pvals) <- names(theta) <- names(V(graph))
  
  if (weighted){
    w <- unlist(calculate.weights(graph, pvalues=pvals, fbh, lambda=lambda))
    pvals <- w[names(pvals)]*pvals
  }
  names(pvals) <- names(V(graph))
  
  BY.results <- BY.proc(pvals, alpha=alpha)
  rej.pvals <- V(graph)[BY.results$Rejected_pvalues]
  n.disc=BY.results$N_Discoveries
  
  if(n.disc>0){
    
    legitimate.set <- is_legitimate_set_outernode(rej.pvals, graph)
    
    while(!(legitimate.set$is.legitimate)){
      
      R.minus <-  rej.pvals[which(legitimate.set$check==TRUE)] 
      if(length(R.minus)>0){
        self.consistent <- is_self_consistent_set(pvals[R.minus], alpha, N = length(pvals))
        if(self.consistent){
          n.disc=length(R.minus)
          rej.pvals=R.minus
          stopcond=TRUE
        }else{
          BY.results = BY.proc(pvals=pvals[R.minus], 
                               alpha = (alpha * length(R.minus))/length(pvals))
          rej.pvals <- R.minus[BY.results$Rejected_pvalues]
          n.disc=BY.results$N_Discoveries
          if(n.disc==0){
            rej.pvals=NULL
            break} else
              legitimate.set <- is_legitimate_set_outernode(rej.pvals, graph)
        }
      } else {
        n.disc=0
        rej.pvals=NULL
        break
      }
      if(stopcond) break
    }
  } 
  false_rej <- sum(theta[rej.pvals] == 0)
  true_rej <- sum(theta[rej.pvals] == 1)
  fdp <- false_rej/max(n.disc,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = n.disc, 
              "Rejections" = rej.pvals%>% na.omit))
}

