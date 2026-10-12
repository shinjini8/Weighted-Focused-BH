# Code to generate figures 1 and 5 and Table 1 from Section 5 of the paper.
library(phyloseq)
library(igraph)
library(structSSI)
library(multtest)
library(MASS)

source("utils.R")
source("Yekutieli.R")
source("WFBH_functions.R")
source("WFBH_related_functions.R")

# STAR implementation:
# Lei et al. (2019), "A General Interactive Framework for False Discovery
# Rate Control under Structural Constraints"
#
# Original implementation available at:
# https://github.com/lihualei71/STAR
#
# The required STAR functions are sourced from the original implementation.

source("STAR_DAG.R")
source("summarize_methods.R")
source("expr_func.R")

#____________________________________________________________________________#

data(GlobalPatterns)
chlamydiae_data <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae") 

#___________________________________________________________________________________#
abundances <- otu_table(chlamydiae_data)
environments <- sample_data(chlamydiae_data)$SampleType
ch.tree <- as_edgelist(as.igraph(phy_tree(chlamydiae_data)))
ch.pval <- treePValues(ch.tree, abundances, environments)
pvalues <- pvalues.y <- ch.pval
theta <- rep(1, length(ch.pval))
graph <- as.igraph(phy_tree(chlamydiae_data))
depths <- find.depth(graph)
vs <- depths$Depths %>% unlist %>% names()
pvalues <- pvalues[vs]; names(theta) <- vs
pvalues.y <- pvalues.y[vs]

#Example code to find Fisher's smoothed p-values
pvalues.fisher <- find.fisher.pvalues(graph, pvalues[names(V(graph))])

graph <- set_vertex_attr(graph, name="pvalues", value = pvalues[names(V(graph))])

graph <- set_vertex_attr(graph, name="pvalues.y", value = pvalues.y[names(V(graph))])

graph <- set_vertex_attr(graph, "theta", value = theta)
graph <- set_vertex_attr(graph, "effective_m", value = unlist(effective.m(graph))[names(V(graph))])
graph <- set_vertex_attr(graph, "effective_l", value = unlist(effective.leaves(graph))[names(V(graph))])

alpha=0.2; lambda=alpha
d.simple <- dagger.proc(graph, type = "simple", alpha)

wsbh.simple <- new.struct.bh.proc(graph=graph, type = "simple", weighted=TRUE, alpha=alpha, fbh=FALSE, lambda=lambda)

sbh.simple <- new.struct.bh.proc(graph=graph, type = "simple", weighted=FALSE, alpha=alpha, fbh=TRUE, lambda=lambda)

#Yekutieli's Procedure
unadj.p.values <- V(graph)$pvalues.y %>%setNames(nm=names(V(graph))) 
G <- graph_as_list(graph)
G1 <- convert_to_Yekutieli_format(G)
y.simple <- Yekutieli(unadj.p.values, G1, (alpha/(2))) %>% Yekutieli.proc(theta=V(graph)$theta)


#BH
adjp <- p.adjust(pvalues, method="BH")
bh.simple <- sum(adjp<=alpha)
bh.simple.rej <- which(adjp<=alpha)
rm(adjp)

#STAR
graph2=graph
V(graph2)$name <- 1:length(graph2)
fun <- create.fun("SeqStep", pstar=0.5)
STARDAG.simple <- STAR.DAG(pvals=V(graph2)$pvalues, DAG=graph2, 
                           print.quiet = TRUE,
                           type="naive",
                           alpha.list = alpha,prop.high = 0.01, prop.low=0.01,
                           fun.list = fun)

nrej <- c(d.simple$No_rejection, wsbh.simple$No_rejection, sbh.simple$No_rejection, y.simple$No_rejection, 
          bh.simple, STARDAG.simple$num_rej)
names(nrej) <- c("Dagger", "WSBH", "SBH", "Yekutieli", "BH", "STAR")

nrej

## EXTRA GRAPHS ##
BHd <- bh.simple.rej%>%names
SBHd <- sbh.simple$Rejections%>%names
WSBHd <- wsbh.simple$Rejections%>%names
dd <- d.simple$Rejections
ykd <- V(graph)$name[y.simple$Rejections]
star <- V(graph)$name[STARDAG.simple$mask==TRUE]

alldisc <- list(BHd, SBHd, WSBHd, dd, ykd, star)
common_all <- base::Reduce(base::intersect, alldisc)
all_disc_union <- base::Reduce(base::union, alldisc)
membership <- sapply(alldisc, function(s) all_disc_union %in% s)
rownames(membership) <- all_disc_union; colnames(membership) <- c("BH", "FBH", "WFBH", "DAGGER", "Yk", "STAR")
membership

mems <- apply(membership, 1, function(x) {
  paste(names(x)[x], collapse = ", ")
})

groups <- split(rownames(membership), mems)

colours <- rep(NA, 41)
names(colours) <- names(pvalues)

colours[groups[[1]]] <- "deeppink3"
colours[groups[[2]]] <- "goldenrod1"
colours[groups[[3]]] <- "blue4"
colours[groups[[6]]] <- "green2"
colours[groups[[4]]] <- "slateblue1"
colours[groups[[5]]] <- "yellow"
colours[groups[[7]]] <- "darkgreen"

sum(is.na(colours))
graph <- set_vertex_attr(graph, "color", value = colours[names(V(graph))])
par(mar=c(0.0,0,1,0)+.0001)
plots <- plot(graph, layout = layout_as_tree, vertex.size = 5,vertex.label=V(graph)$pvalues%>%round(digits=3),
              vertex.label.dist = 1.2, edge.color="black",
              edge.arrow.size=0.22, edge.arrow.width=0.8, margin=0, vertex.color=V(graph)$color)

legend("bottomleft",                       # Position of the legend
       legend = c("No rejection", "STAR only",
                  "STAR, BH", "STAR, BH, WFBH", 
                  "STAR, BH, WFBH, FBH, Yk", "STAR, BH, WFBH, FBH", 
                  "STAR, BH, WFBH, DAGGER", 
                  "STAR, BH, WFBH, FBH,  DAGGER, Yk"),  # Labels for the legend
       pch=rep(21,7), pt.cex=1.5, col=rep("black",7),
       pt.bg = c("white", "darkgreen",
                 "slateblue1", "green2", "blue4", "goldenrod1", "yellow", 
                 "deeppink3"),  # Colors matching the node colors
       title = "",cex = 0.8,                                 # Make the text in the legend smaller
       title.cex = 0.9,     )  # Title of the legend


#-----------------------------------------------------------------------------------------------------#

# Filtering graph
colours2 <- rep(NA, 41); shapes <- rep("circle", 41)
names(colours2) <- names(shapes) <- names(pvalues)

colours2[WSBHd] <- "blue4"

bhd2 <- BHd[!(BHd %in% WSBHd)]
colours2[bhd2] <- "slateblue1"

leafs <- find.leaves(graph, all=TRUE)
outernodes <- leafs[leafs %in% BHd]
shapes[outernodes] <- "square"

sum(is.na(colours2))

graph <- set_vertex_attr(graph, "colorfilter", value = colours2[names(V(graph))])
graph <- set_vertex_attr(graph, "shape", value = shapes[names(V(graph))])
graph <- set_vertex_attr(graph, "labels", value = labels[names(V(graph))])

par(mar=c(0.0,0,1,0)+.0001)
plots <- plot(graph, layout = layout_as_tree, vertex.size = 5,vertex.label=NA, #V(graph)$labels,
              vertex.label.dist = 1, edge.color="black",
              edge.arrow.size=0.22, edge.arrow.width=0.8, margin=0, 
              vertex.color=V(graph)$colorfilter,
              vertex.shapes=V(graph)$shape)
