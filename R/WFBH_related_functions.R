
library(igraph)


find.depth <- function(graph){
  #Calculates levels of a graph and the members in each level
  # Returns a list containing the maximum depth, node names grouped by depth,
  # and a data frame with each node's depth (d), position within that depth (i),
  # and node value (Member.Hypotheses).
  tot.nodes <- length(graph)
  if(tot.nodes>0)
    levels   = list() else stop("Empty graph!")
  i=1
  df <- data.frame("d" = 0, "i" = 0, "Member Hypotheses"=0)
  while(length(graph) > 0){
    levels[[i]] <- which(igraph::degree(graph, mode="in")==0)
    graph <- delete_vertices(graph, levels[[i]])
    df_level <- data.frame("d" = i, "i" = 1:length(levels[[i]]), "Member Hypotheses"  = levels[[i]])
    df <- rbind(df, df_level)
    i = i+1
  }
  df <- df[-1,]

  df <- df[order(df$d),]
  return(list("Max_depth" =length(levels), "Depths" = levels, "df"=df))
}  

find.parents <- function(graph, vertex){
  #Finds parents of a given node of a graph
  root <- which(V(graph) == V(graph)[vertex])
  return(neighbors(graph,vertex, mode="in"))
  
}
find.children <- function(graph, vertex){
  #Finds children of a given node of a graph
  root <- which(V(graph) == V(graph)[vertex])
  return(neighbors(graph,vertex, mode="out"))
  
}

m.d <- function(graph){ #calculates effective number of nodes at depth d
  depths <- find.depth(graph)
  n_parents <- lapply(depths$Depths, function(x) 
    sapply(names(x), function(y) {
      l <- find.parents(graph,y)%>% length
      if(l==0) # Just making sure the first layer has 1 parent each, so no NaN in results
        l=1
      return(l)
      }))
  eff.nodes <- sapply(n_parents, function(x) sum(1/x))
  return(eff.nodes)
}
n.d <- function(graph){
  #finds the number of groups at each depth of the graph
  #This is defined as the number of nodes at depths 1:d-1 that have at least one child in depth d
  #returns a vector
  depths <- find.depth(graph)
  D <- depths$Max_depth
  n.d <- c()
  n.d[1] <- 1 #at depth 1, there is one group of nodes
  if(D>1){
    for(i in 2:D){
      parents <- unlist(depths$Depths[1:(i-1)])
      children_at_d <- sapply(names(parents), function(x) {
        children <- find.children(graph, x)
        sum(names(children) %in% names(depths$Depths[[i]]))%>%return()
      })
      n.d[i] <- sum(children_at_d>0)
    }}
  return(n.d)
}

find.leaves <- function(graph, vertex=NA, all=FALSE){
  #vertex = Name of the vertex
  # Set all = TRUE if all leaf nodes of the graph are required
  roots <- V(graph)$name[degree(graph = graph, v = V(graph), mode = "in") == 0]
  
  #Determine the name(s) of the leaf vertex/vertices.
  leafs <- V(graph)$name[degree(graph = graph, v = V(graph), mode = "out") == 0]
  
  #Calculate the graph depth. That is, the largest path length of all shortest paths from root 
  #to leaf vertices.
  sh_paths <- sapply(roots, function(x) 
    shortest_paths(graph = graph, from = x, to = leafs)$vpath )
  max_sh_path_length <- max(sapply(X = sh_paths, FUN = length))
  
  #If 'vertex' is a leaf node itself, then return the name of that vertex.
  if(!is.na(vertex)){
    if(vertex %in% leafs){
      return(vertex)
      } else {
        #If 'vertex' is not a leaf node, then determine the subset of all vertices that are in the 
        #neighborhood of 'vertex', excluding 'vertex' itself ('mindist' = 1). The maximum order or 
        #'depth' is specified by ('max_sh_path_length')
        vertices_subset <- neighborhood(graph = graph, 
                                        order = max_sh_path_length - 1, 
                                        nodes = V(graph)[V(graph)$name == vertex], 
                                        mode = "out",
                                        mindist = 1)
        #Extract the names of the vertices in 'vertices_subset'
        vertices_subset_names <- names(unlist(vertices_subset))
        
        #The overlap/intersection of vertex names between 'vertices_subset_names' and 'leafs' gives 
        #all leaf nodes that are part of the subtree with root vertex 'vertex'.
        result <- intersect(x = vertices_subset_names, leafs)
        
        return(result)
        }
    } else if (is.na(vertex) & isTRUE(all))
      return(leafs)
}

samples.fun <- function(graph, mu.vector, alt.indices, d=1.5,
                        case_name=c("global", "incremental", "decremental"), 
                        rho, mu=2) {
  #According to the case, assigns mu to the alternate indices and then finds the samples
  deps <- find.depth(graph)
  D <- deps$Max_depth
  if(case_name =="incremental"){
    mu.vector[alt.indices] <- sapply(alt.indices, function(x) 
      mu + d*(D-sapply(deps$Depths, function(y) x %in% names(y)) %>% which()) )  
  } else if(case_name =="decremental"){
    mu.vector[alt.indices] <- sapply(alt.indices, function(x) 
      mu + d*(sapply(deps$Depths, function(y) x %in% names(y)) %>% which()-1) )  
  }else if(case_name=="global"){
    mu.vector[alt.indices] = mu}
    x <- rnorm(1)
  samples <- (mu.vector + ((1-rho)*rnorm(length(mu.vector))) +(rho*x)) %>% return()
}

altindices <- function(adj_mat, leaf_alt_density){
  # Identifies the indices of non-null nodes in the graph.
  # Inputs: the graph adjacency matrix and the proportion of non-null leaves.
  # Returns: a list containing the graph and the indices of all non-null nodes.
  graph <- graph_from_adjacency_matrix( adj_mat )
  deps <- find.depth(graph)
  D <- deps$Max_depth
  leafs <- find.leaves(graph, vertex=NA, all=TRUE)
  L <- length(leafs)
  alt.indices <- leafs[sample(1:L, round(leaf_alt_density*L))]
  for(d in (D-1):1){
    a <- sapply(deps$Depths[[d]]%>%names, function(x){
      children <- find.children(graph, x)%>%names
      if(any(children %in% alt.indices)) 
        alt.indices <- c(x, alt.indices)
      return(alt.indices)})
    alt.indices <- Reduce(union, a)
  }
  return(list("graph"=graph, "alt.indices"=alt.indices))
  
}

create_graph <- function( graph, mu,
                          rho=0, alt.indices, 
                          case = c("global", "incremental", "decremental"),
                          d=1.5){
  # Inputs: a graph, the non-null mean, the common correlation coefficient,
  # the indices of non-null nodes, and the setup (global, incremental, or decremental).
  # d specifies the change in the non-null mean between consecutive depths
  # in the incremental and decremental setups.
  # Computes p-values according to models specified in Section 4 and Web Appendix G of the manuscript
  # Returns the graph with smoothed p-values, effective m and l as in Ramdas et al(2019) 
  # and the original p-values ina list.
  N <- length(graph)
  theta <- mu.vector <- rep(0,N)
  pvalues <- rep(NA,N)
  V(graph)$name <- names(theta) <- names(mu.vector) <- 
    names(pvalues) <- names(V(graph))
  theta[alt.indices] <- 1
  
  pvalues[1:N] <- 1-pnorm(samples.fun(graph, mu.vector=mu.vector, alt.indices=alt.indices, d=d, case_name = case,rho=rho, mu=mu))
  pvalues.fisher <- find.fisher.pvalues(graph, pvalues)
  pvalues.simes <- find.simes.pvalues(graph, pvalues)
    
  graph <- set_vertex_attr(graph, name="pvalues", value = pvalues)
  graph <- set_vertex_attr(graph, "pvalues.fisher", value = pvalues.fisher)
  graph <- set_vertex_attr(graph, "pvalues.simes", value = pvalues.simes)
  graph <- set_vertex_attr(graph, "theta", value = theta)
  graph <- set_vertex_attr(graph, "effective_m", value = sapply(effective.m(graph), unlist)%>% unlist)
  graph <- set_vertex_attr(graph, "effective_l", value = sapply(effective.leaves(graph), unlist)%>% unlist)

  return(list("graph"=graph, "pvalues"=pvalues))
}


find.fisher.pvalues <- function(graph, pvalues){
  # Computes Fisher's smoothed p-values from the graph and original p=values
  pvalues.subtree <- lapply(V(graph), function(x) subcomponent(graph, x, mode= "out")) %>%
    lapply(function(x) pvalues[names(x)]) 
  pvalues.subtree <- pvalues.subtree[sapply(pvalues.subtree, length)>0] # just removing the empty lists
  fishers.pvalues <- sapply(pvalues.subtree, function(x) {
    if(length(x)==1) # make sure to do no smoothing on the leaf pvalues
      return(x) else
        return((-2*sum(log(x))) %>% pchisq(df = (2*length(x)), lower.tail=F))
    })
  names(fishers.pvalues) <- names(pvalues.subtree)
  return(fishers.pvalues)
}

find.simes.pvalues <- function(graph, pvalues){
  # Computes Simes's smoothed p-values from the graph and original p=values
  pvalues.subtree <- lapply(V(graph), function(x) subcomponent(graph, x, mode= "out")) %>%
    lapply(function(x) pvalues[names(x)]) 
  pvalues.subtree <- pvalues.subtree[sapply(pvalues.subtree, length)>0] # just removing the empty lists
  simes.pvalues <- sapply(pvalues.subtree, function(x) {
    if(length(x)==1) # make sure to do no smoothing on the leaf pvalues
      return(x) else{
        sp <- min(length(x)*sort(x)/(1:length(x)))
        return(sp)}
    })
  names(simes.pvalues) <- names(pvalues.subtree)
  return(simes.pvalues)
}

find.bonferroni.pvalues <- function(graph, pvalues){
  # Computes Bonferroni's smoothed p-values from the graph and original p=values
  pvalues.subtree <- lapply(V(graph), function(x) subcomponent(graph, x, mode= "out")) %>%
    lapply(function(x) pvalues[names(x)]) 
  pvalues.subtree <- pvalues.subtree[sapply(pvalues.subtree, length)>0] # just removing the empty lists
  tip.pvalues <- sapply(pvalues.subtree, function(x) min(x)*length(x))
  names(tip.pvalues) <- names(pvalues.subtree)
  return(tip.pvalues)
}

effective.m <- function(graph){
  # Returns effective number of nodes from bottom to top in a list as used in the DAGGER procedure in Ramdas et al(2019) 
  depths <- find.depth(graph)
  D <- depths$Max_depth
  m <- vector("list", D)
  
  m[[D]] <- rep(1, length(depths$Depths[[D]])) #assign m=1 to each leaf-node
  names(m[[D]]) <- depths$Depths[[D]]%>%names
  
  if(D>1){
    for(d in (D-1):1){
      nodes_at_d <- depths$Depths[[d]]%>%names
      names(nodes_at_d) <- nodes_at_d
      children <- lapply(nodes_at_d, function(x) find.children(graph,x) )
      #names(children) <- nodes_at_d
      m_children <- lapply(children, function(x) {
        munl <- unlist(m[d+1:D]) #vector of all m values found before dth step
        mch <- munl[names(munl) %in% names(x)] 
        if (length(mch)==0) {mch=0 #if a leaf is found, set m value=0, later 1 will be added
        names(mch) ="leaf"}
        return(mch)
      })%>% setNames(nodes_at_d)
      
      n_parents <- lapply(children, function(x) {
        p <- sapply(names(x), function(y) find.parents(graph, y) %>% length )
        if(length(p)==0) {p=1; names(p)="leaf"}
        return(p)
        #return(setNames(p, x) )
      }) %>%  setNames(nodes_at_d)
      m[[d]] <- 1+ mapply(function(x,y) 
        sum(as.numeric(x)/as.numeric(y[names(x)])), m_children, n_parents)
    }}
  return(m)}

effective.leaves <- function(graph){
  # Returns effective number of leaves from bottom to top in a list as used in the DAGGER procedure in Ramdas et al(2019) 
  tot.nodes <- length(graph)
  if(tot.nodes>0)
    leaves  = c() else stop("Empty graph!")
  depths <- find.depth(graph)
  D <- depths$Max_depth
  l <- vector("list", D)
  l[[D]] <- rep(1, length(depths$Depths[[D]])) #assign l=1 to each leaf-node
  names(l[[D]]) <- depths$Depths[[D]]%>%names
  
  if(D>1){
    for(d in (D-1):1){
      nodes_at_d <- depths$Depths[[d]]%>%names
      children <- lapply(nodes_at_d, function(x) find.children(graph,x) )
      names(children) <- nodes_at_d
      n_leaves <- lapply(children, function(x) { 
        lunl <- unlist(l[d+1:D]) #vector of all leaves before dth step
        lch <- lunl[names(lunl) %in% names(x)]
        if (length(lch)==0) {lch=1 #if a leaf is found, set l value=1
        names(lch) = "leaf"}
        return(lch)}  )
      
      n_parents <- lapply(children, function(x) {
        p <- sapply(names(x), function(y) find.parents(graph, y) %>% length )
        if(length(p)==0) {p=1; names(p)="leaf"}
        return(p)
      }) %>%  setNames(nodes_at_d)
      l[[d]] <-  mapply(function(x,y) 
        sum(as.numeric(x)/as.numeric(y[names(x)])), n_leaves, n_parents)
    }
  }
  
  return(l)
} 


stepup.proc <- function (pvals, m, l, L, R_1_d, alpha) {
  # Stepup procedure as used in the DAGGER procedure in Ramdas et al(2019) 
  n_d <- length(pvals)
  levs <- alpha* (l/L)* (m + R_1_d + (1:n_d) - 1)/m
  R <- max(c(which(ifelse(sort(pvals) < levs, 1, 0) == 1), 
             0))
  rej.pvals <- order(pvals)[1:R]
  if(R==0)
    Discoveries = NA else Discoveries = sort(rej.pvals)
  list(No_Discoveries = R, Discoveries=Discoveries)
}  
stepup.proc.by <- function (pvals, m, d, hd, l, L, R_1_d, alpha) {
  
  # Stepup procedure with reshaping as used in the DAGGER procedure in Ramdas et al(2019) 
  n_d <- length(pvals)
  s <- sapply(m, function(x) (x+d-1): (x+hd-1))
  if(!is.null(dim(s)))
  k <- apply(s, 2, function(x) sum(1/x)) else k <- sum(1/s)
  levs <- alpha* (l/L)* (m + R_1_d + (1:n_d) - 1)/(m*k)
  R <- max(c(which(ifelse(sort(pvals) < levs, 1, 0) == 1), 
             0))
  rej.pvals <- order(pvals)[1:R]
  if(R==0)
    Discoveries = NA else Discoveries = sort(rej.pvals)
  list(No_Discoveries = R, Discoveries=Discoveries)
}  
dagger.proc <- function(graph, type, alpha){
  # Implements the DAGGER procedure in Ramdas et al(2019) 
  depths <- find.depth(graph)
  theta <- V(graph)$theta
  effective_m <- V(graph)$effective_m
  effective_l <- V(graph)$effective_l
  if (type == "simple")
    pvalues <- V(graph)$pvalues else if (type == "fisher")
      pvalues <- V(graph)$pvalues.fisher else if (type == "simes")
          pvalues <- V(graph)$pvalues.simes else if (type == "bonferroni")
            pvalues <- V(graph)$pvalues.bonferroni
  names(pvalues) <- names(theta) <- names(effective_m) <- names(effective_l) <- names(V(graph))
  
  R.list <- c()
  L <- find.leaves(graph, vertex=NA, all=TRUE) %>% length
  
  nodes_at_d <- depths$Depths[[1]] %>% names
  R_1_d <- 0
  apply_stepup <- stepup.proc(pvals = pvalues[nodes_at_d], 
                              m = effective_m[nodes_at_d],
                              l = effective_l[nodes_at_d],
                              L,
                              R_1_d, alpha =alpha)
  detected_nodes <- nodes_at_d[apply_stepup$Discoveries]
  R.list <- union(R.list, detected_nodes%>% na.omit())
  R_1_d <- length(R.list)
  
  
  for(d in 2:depths$Max_depth){
    if(apply_stepup$No_Discoveries ==0){
      #print(paste("No more discoveries at depth", d-1))
      break} else {
        nodes_at_d <- lapply(detected_nodes, function(x){
          children <- find.children(graph, x)%>%names
          children[children %in% names(depths$Depths[[d]]) ]%>% return()
          }) 
        
        check_rej_parents <- lapply(nodes_at_d, function(x){
          sapply(x, function(y){
            parents <- find.parents(graph, y)
            sapply(parents%>%names, function(z) z %in% R.list) %>% prod})
        })
        nodes_at_d <- unlist(nodes_at_d)[unlist(check_rej_parents)>0]%>%unique()
        apply_stepup <- stepup.proc(pvals = pvalues[nodes_at_d], 
                                    m = effective_m[nodes_at_d],
                                    l = effective_l[nodes_at_d],
                                    L=L,
                                    R_1_d=R_1_d, alpha = alpha)
        
        detected_nodes <- nodes_at_d[apply_stepup$Discoveries]
        R.list <- union(R.list, detected_nodes%>% na.omit())
        R_1_d <- length(R.list)
      }
  }
  false_rej <- sum(theta[R.list] == 0)
  true_rej <- sum(theta[R.list] == 1)
  fdp <- false_rej/max(R_1_d,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = R_1_d, "Rejections" = R.list%>% na.omit))
}


dagger.proc.by <- function(graph, type, alpha){
  # Implements the reshaped DAGGER procedure in Ramdas et al(2019) 
  depths <- find.depth(graph)
  theta <- V(graph)$theta
  effective_m <- V(graph)$effective_m
  effective_l <- V(graph)$effective_l
  if (type == "simple")
    pvalues <- V(graph)$pvalues else if (type == "fisher")
      pvalues <- V(graph)$pvalues.fisher else if (type == "simes")
          pvalues <- V(graph)$pvalues.simes else if (type == "bonferroni")
            pvalues <- V(graph)$pvalues.bonferroni
  names(pvalues) <- names(theta) <- names(effective_m) <- names(effective_l) <- names(V(graph))
  
  R.list <- c()
  L <- find.leaves(graph, vertex=NA, all=TRUE) %>% length
  
  nodes_at_d <- depths$Depths[[1]] %>% names
  R_1_d <- 0
  d=1
  apply_stepup <- stepup.proc.by(pvals = pvalues[nodes_at_d], 
                              m = effective_m[nodes_at_d],
                              d=1,
                              hd = depths$Depths[[d]]%>%length(),
                              l = effective_l[nodes_at_d],
                              L,
                              R_1_d = R_1_d, alpha =alpha)
  detected_nodes <- nodes_at_d[apply_stepup$Discoveries]
  R.list <- union(R.list, detected_nodes%>% na.omit())
  R_1_d <- length(R.list)
  
  
  for(d in 2:depths$Max_depth){
    if(apply_stepup$No_Discoveries ==0){
      #print(paste("No more discoveries at depth", d-1))
      break} else {
        nodes_at_d <- lapply(detected_nodes, function(x){
          children <- find.children(graph, x)%>%names
          children[children %in% names(depths$Depths[[d]]) ]%>% return()
        }) 
        
        check_rej_parents <- lapply(nodes_at_d, function(x){
          sapply(x, function(y){
            parents <- find.parents(graph, y)
            sapply(parents%>%names, function(z) z %in% R.list) %>% prod})
        })
        nodes_at_d <- unlist(nodes_at_d)[unlist(check_rej_parents)>0]%>%unique()
        apply_stepup <- stepup.proc.by(pvals = pvalues[nodes_at_d], 
                                    m = effective_m[nodes_at_d],
                                    d=d,
                                    hd = depths$Depths[[d]]%>%length(),
                                    l = effective_l[nodes_at_d],
                                    L=L,
                                    R_1_d=R_1_d, alpha = alpha)
        
        detected_nodes <- nodes_at_d[apply_stepup$Discoveries]
        R.list <- union(R.list, detected_nodes%>% na.omit())
        R_1_d <- length(R.list)
        }
      }
  
  false_rej <- sum(theta[R.list] == 0)
  true_rej <- sum(theta[R.list] == 1)
  fdp <- false_rej/max(R_1_d,1)
  pow <- true_rej/sum(theta)
  return(list("fdp" = fdp, "pow"=pow, "No_rejection" = R_1_d, "Rejections" = R.list%>% na.omit))
}
#_______________________________________________________________________________#

