# Generates simulated data corresponding to FDR and Power
# for wide tree considered in Section 4.
# Vary non-null proportion in [0.01, 1] while independent samples are drawn with rho=0
set.seed(250)

library(igraph)
library(dplyr)
library(gtools)
source("utils.R")
source("WFBH_related_functions.R")
load("widetree_adj_mat.RData")

n.rep <- 200
type.list <- c("simple", "fisher")

case.list <- c("global", "incremental", "decremental")
density.list <- seq(0.1, 1, length.out=10) #proportion of non-null leaves
mu=2 
rho=0
alpha=0.05 #FDR threshold
lambda=0.5 # lambda for WFBH


fdp.array.dagger <- pow.array.dagger <- array(NA, dim=c(length(type.list), 
                                                        length(case.list), 
                                                        length(density.list), 
                                                        n.rep))

fdp.array.fbh <- pow.array.fbh <- array(NA, dim=c(length(type.list), 
                                                  length(case.list),
                                                  length(density.list),
                                                  n.rep))
fdp.array.wfbh <- pow.array.wfbh <- array(NA, dim=c(length(type.list), 
                                                    length(case.list),
                                                    length(density.list),
                                                    n.rep))
fdp.array.star <- pow.array.star <- array(NA, dim=c(length(type.list), 
                                                    length(case.list),
                                                    length(density.list),
                                                    n.rep))

fdp.array.y <- pow.array.y <- array(NA, dim=c(length(type.list), 
                                              length(case.list),
                                              length(density.list),
                                              n.rep))
for(i in 1:length(case.list)){
  case = case.list[i]
  for(j in 1:length(density.list)){
    density=density.list[j]
    alt_indices <- suppressWarnings(altindices(adj_mat=adj_mat, leaf_alt_density=density))
    f = function() {
      lapply(1:n.rep, function(x){
        library(igraph)
        library(mvtnorm)
        library(multtest)
        library(structSSI)
        library(MASS)
        source("utils.R")
        source("Yekutieli.R")
        source("WFBH_related_functions.R")
        source("WFBH_functions.R")
        source("STAR_DAG.R")
        source("summarize_methods.R")
        source("expr_func.R")
        g <- create_graph(graph=alt_indices$graph, mu=mu, rho=rho, 
                          alt.indices=alt_indices$alt.indices, case=case)
        d.simple <- dagger.proc(g$graph, type = "simple", alpha)
        d.fisher <- dagger.proc(g$graph, type = "fisher", alpha)
        
        wfbh.simple <- new.struct.bh.proc(graph=g$graph, type = "simple", weighted=TRUE, fbh=FALSE, alpha=alpha, lambda=lambda)
        wfbh.fisher <- new.struct.bh.proc(graph=g$graph, type = "fisher", weighted=TRUE, fbh=FALSE, alpha=alpha, lambda=lambda)
        
        fbh.simple <- new.struct.bh.proc(graph=g$graph, type = "simple", weighted=FALSE, alpha=alpha, fbh=TRUE, lambda=lambda)
        fbh.fisher <- new.struct.bh.proc(graph=g$graph, type = "fisher", weighted=FALSE, alpha=alpha, fbh=TRUE, lambda=lambda)

        #Yekutieli Procedure
        unadj.p.values <- V(g$graph)$pvalues 
        names(unadj.p.values) <- V(g$graph)$name
        G <- graph_as_list(g$graph)
        G1 <- convert_to_Yekutieli_format(G)
        y.simple <- Yekutieli(unadj.p.values, G1, (alpha/(2*1.44))) %>% Yekutieli.proc(theta=V(g$graph)$theta)
        
        unadj.p.values <- V(g$graph)$pvalues.fisher 
        names(unadj.p.values) <- V(g$graph)$name
        G <- graph_as_list(g$graph)
        G1 <- convert_to_Yekutieli_format(G)
        y.fisher <- Yekutieli(unadj.p.values, G1, (alpha/(2*1.44))) %>% Yekutieli.proc(theta=V(g$graph)$theta)
        
        #STAR
        graph2=g$graph
        V(graph2)$name <- 1:length(graph2)
        fun <- create.fun("SeqStep", 0.5)
        STARDAG.simple <- STAR.DAG(pvals=unname(g$pvalues), DAG=graph2, print.quiet = TRUE,
                                   alpha.list = alpha, prop.high = (1/length(graph2)), prop.low=(1/length(graph2)),
                                   fun.list = fun)
        STARDAG.simple$pow <- pow.STAR(theta = V(graph2)$theta, mask= STARDAG.simple$mask)
        
        STARDAG.fisher <- STAR.DAG(pvals=unname(V(graph2)$pvalues.fisher), DAG=graph2, print.quiet = TRUE,
                                   alpha.list = alpha, prop.high = (1/length(graph2)), prop.low=(1/length(graph2)),
                                   fun.list = fun)
        STARDAG.fisher$pow <- pow.STAR(theta = V(graph2)$theta, mask= STARDAG.fisher$mask)
        
        return(list("d.simple" = d.simple, 
                    "fbh.simple" = fbh.simple,
                    "wfbh.simple" = wfbh.simple,
                    "star.simple" = STARDAG.simple,
                    "yekutieli.simple" = y.simple,
                    "d.fisher" = d.fisher, 
                    "fbh.fisher" = fbh.fisher,
                    "wfbh.fisher" = wfbh.fisher,
                    "star.fisher" = STARDAG.fisher,
                    "yekutieli.fisher" = y.fisher
        ))})
    }
    
    result <- f() %>% suppressMessages()
    fdp.array.dagger[1,i,j,] <- sapply(result, function(x) x$d.simple$fdp)
    pow.array.dagger[1,i,j,] <- sapply(result, function(x) x$d.simple$pow)
    
    fdp.array.dagger[2,i,j,] <- sapply(result, function(x) x$d.fisher$fdp)
    pow.array.dagger[2,i,j,] <- sapply(result, function(x) x$d.fisher$pow)
    
    fdp.array.fbh[1,i,j,] <- sapply(result, function(x) x$fbh.simple$fdp)
    pow.array.fbh[1,i,j,] <- sapply(result, function(x) x$fbh.simple$pow)
    
    fdp.array.wfbh[1,i,j,] <- sapply(result, function(x) x$wfbh.simple$fdp)
    pow.array.wfbh[1,i,j,] <- sapply(result, function(x) x$wfbh.simple$pow)
    
    fdp.array.fbh[2,i,j,] <- sapply(result, function(x) x$fbh.fisher$fdp)
    pow.array.fbh[2,i,j,] <- sapply(result, function(x) x$fbh.fisher$pow)
    
    fdp.array.wfbh[2,i,j,] <- sapply(result, function(x) x$wfbh.fisher$fdp)
    pow.array.wfbh[2,i,j,] <- sapply(result, function(x) x$wfbh.fisher$pow)
    
    fdp.array.star[1,i,j,] <- sapply(result, function(x) x$star.simple$fdp)
    pow.array.star[1,i,j,] <- sapply(result, function(x) x$star.simple$pow)
    
    fdp.array.star[2,i,j,] <- sapply(result, function(x) x$star.fisher$fdp)
    pow.array.star[2,i,j,] <- sapply(result, function(x) x$star.fisher$pow)

    fdp.array.y[1,i,j,] <- sapply(result, function(x) x$yekutieli.simple$fdp)
    pow.array.y[1,i,j,] <- sapply(result, function(x) x$yekutieli.simple$pow)
    
    fdp.array.y[2,i,j,] <- sapply(result, function(x) x$yekutieli.fisher$fdp)
    pow.array.y[2,i,j,] <- sapply(result, function(x) x$yekutieli.fisher$pow)
    print(j)
  }
  print(i)
}

FDR.fbh <- apply(fdp.array.fbh, c(1,2,3), mean)
FDR.wfbh <- apply(fdp.array.wfbh, c(1,2,3), mean) 
FDR.dagger <- apply(fdp.array.dagger, c(1,2,3), mean)
FDR.star <- apply(fdp.array.star, c(1,2,3), mean)
FDR.y <- apply(fdp.array.y, c(1,2,3), mean)
Power.fbh <- apply(pow.array.fbh, c(1,2,3), mean)
Power.wfbh <- apply(pow.array.wfbh, c(1,2,3), mean)
Power.dagger <- apply(pow.array.dagger, c(1,2,3), mean)
Power.star <- apply(pow.array.star, c(1,2,3), mean)
Power.y <- apply(pow.array.y, c(1,2,3), mean)
dimnames(FDR.fbh) <- dimnames(FDR.wfbh) <- dimnames(FDR.dagger) <- dimnames(FDR.y) <-dimnames(FDR.star) <-  
  dimnames(Power.fbh) <- dimnames(Power.wfbh) <- dimnames(Power.dagger) <- dimnames(Power.y) <-dimnames(Power.star) <-  
  list(paste(type.list), paste(case.list), paste(density.list))
#_____________________________________________________________________________________#


pow.simple.data <- pow.fisher.data <- c()
fdr.simple.data <- fdr.fisher.data <- c()
for(i in 1:length(case.list)){
  melt.fdrsimple.data <- melt.powsimple.data <-c()
  melt.fdrfisher.data <- melt.powfisher.data <-c()
  
  nam.fdrsimple <- paste0("simplefdr.", case.list[i])
  assign(nam.fdrsimple, rbind("FBH" = FDR.fbh[1,i,], 
                              "WFBH" = FDR.wfbh[1,i,],
                              "DAGGER" = FDR.dagger[1,i,],
                              "STAR" = FDR.star[1,i,],
                              "Yk" = FDR.y[1,i,])%>% round(digits=3) %>% t %>%
           'rownames<-' (round(density.list, digits=3)) %>% as.data.frame %>% tibble::rownames_to_column( "p"))
  
  nam.fdrfisher <- paste0("fisherfdr.", case.list[i])
  assign(nam.fdrfisher, rbind("FBH" = FDR.fbh[2,i,], 
                              "WFBH" = FDR.wfbh[2,i,], 
                              "DAGGER" = FDR.dagger[2,i,],
                              "STAR" = FDR.star[2,i,],
                              "Yk" = FDR.y[2,i,])%>% round(digits =3) %>% t  %>%
           'rownames<-' (round(density.list, digits=3))%>% as.data.frame %>% tibble::rownames_to_column( "p"))
  
  
  #_____________________________________________________________________________________#
  nam.powsimple <- paste0("simplepow.", case.list[i])
  assign(nam.powsimple, rbind("FBH" = Power.fbh[1,i,],
                              "WFBH" = Power.wfbh[1,i,], 
                              "DAGGER" = Power.dagger[1,i,],
                              "STAR" = Power.star[1,i,],
                              "Yk" = Power.y[1,i,])%>% round(digits =3) %>% t %>%
           'rownames<-' (round(density.list, digits=3))%>% as.data.frame %>% tibble::rownames_to_column( "p"))
  
  nam.powfisher <- paste0("fisherpow.", case.list[i])
  assign(nam.powfisher, rbind("FBH" = Power.fbh[2,i,],
                              "WFBH" = Power.wfbh[2,i,],
                              "DAGGER" = Power.dagger[2,i,],
                              "STAR" = Power.star[2,i,],
                              "Yk" = Power.y[2,i,])%>% round(digits =3) %>% t %>%
           'rownames<-' (round(density.list, digits=3))%>% as.data.frame %>% tibble::rownames_to_column( "p"))
  
  #_____________________________________________________________________________________#
  # 5 methods. change 5 to the number of methods being used
  nam.fdrsimple.melt <- paste0("fdrmelt.simple.", case.list[i])
  assign(nam.fdrsimple.melt, melt(get(nam.fdrsimple)) %>% 
           `colnames<-`(c("p", "Method", "FDR")) %>% 
           cbind("Smooth"=rep("Original", (5*dim(get(nam.fdrsimple))[1]))))
  
  nam.fdrfisher.melt <- paste0("fdrmelt.fisher.", case.list[i])
  assign(nam.fdrfisher.melt, melt(get(nam.fdrfisher)) %>% 
           `colnames<-`(c("p", "Method", "FDR")) %>% 
           cbind("Smooth"=rep("Fisher", (5*dim(get(nam.fdrfisher))[1]))))
  
  nam.powsimple.melt <- paste0("powmelt.simple.", case.list[i])
  assign(nam.powsimple.melt, melt(get(nam.powsimple)) %>% 
           `colnames<-`(c("p", "Method", "Power")) %>% 
           cbind("Smooth"=rep("Original", (5*dim(get(nam.powsimple))[1]))))
  
  nam.powfisher.melt <- paste0("powmelt.fisher.", case.list[i])
  assign(nam.powfisher.melt, melt(get(nam.powfisher)) %>% 
           `colnames<-`(c("p", "Method", "Power")) %>% 
           cbind("Smooth"=rep("Fisher", (5*dim(get(nam.powfisher))[1]))))
  
  #_______________________________________________________________________________________________________#
  
  
  melt.fdrsimple.data <- rbind(melt.fdrsimple.data, get(nam.fdrsimple.melt))
  fdr.simple.data <- rbind(fdr.simple.data, melt.fdrsimple.data)
  melt.powsimple.data <- rbind(melt.powsimple.data, get(nam.powsimple.melt))
  pow.simple.data <- rbind(pow.simple.data, melt.powsimple.data)
  
  melt.fdrfisher.data <- rbind(melt.fdrfisher.data, get(nam.fdrfisher.melt))
  fdr.fisher.data <- rbind(fdr.fisher.data, melt.fdrfisher.data)
  melt.powfisher.data <- rbind(melt.powfisher.data, get(nam.powfisher.melt))
  pow.fisher.data <- rbind(pow.fisher.data, melt.powfisher.data)
  
}
save(fdr.simple.data, fdr.fisher.data, pow.simple.data, pow.fisher.data, file="widetree_plotdata.RData")