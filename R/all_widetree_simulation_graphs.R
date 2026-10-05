# Generates FDR and power plots for the wide-tree simulations in Section 4.
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)

#_______GRAPHS FOR WIDE TREE INDEPENDENCE ________#
load("widetree_plotdata.RData")
colours <- c("green2", "darkgreen", "royalblue1", "darkorchid", "goldenrod2")
#...... SIMPLE ORIGINAL P-VALUES ........#
widetree.fdr.simple.glo.ind <- fdr.simple.data[1:50,]
widetree.fdr.simple.inc.ind <- fdr.simple.data[51:100,]
widetree.fdr.simple.dec.ind <- fdr.simple.data[101:150,]

widetree.pow.simple.glo.ind <- pow.simple.data[1:50,]
widetree.pow.simple.inc.ind <- pow.simple.data[51:100,]
widetree.pow.simple.dec.ind <- pow.simple.data[101:150,]


widetree.fdr.fisher.glo.ind <- fdr.fisher.data[1:50,]
widetree.fdr.fisher.inc.ind <- fdr.fisher.data[51:100,]
widetree.fdr.fisher.dec.ind <- fdr.fisher.data[101:150,]

widetree.pow.fisher.glo.ind <- pow.fisher.data[1:50,]
widetree.pow.fisher.inc.ind <- pow.fisher.data[51:100,]
widetree.pow.fisher.dec.ind <- pow.fisher.data[101:150,]

glo.fdr <- rbind(widetree.fdr.simple.glo.ind, widetree.fdr.fisher.glo.ind)
dec.fdr <- rbind(widetree.fdr.simple.dec.ind, widetree.fdr.fisher.dec.ind)
inc.fdr <- rbind(widetree.fdr.simple.inc.ind, widetree.fdr.fisher.inc.ind)

glo.pow <- rbind(widetree.pow.simple.glo.ind, widetree.pow.fisher.glo.ind)
dec.pow <- rbind(widetree.pow.simple.dec.ind, widetree.pow.fisher.dec.ind)
inc.pow <- rbind(widetree.pow.simple.inc.ind, widetree.pow.fisher.inc.ind)

# FDR plots
plot.widetree.fdr.glo.ind <- ggplot(glo.fdr, aes(x=p, y=FDR, color = Method, 
                                                 group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "FDR")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          
    legend.title = element_text(size = 8),         
    axis.title = element_text(size = 8),           
    axis.text = element_text(size = 5),            
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="Wide Tree", subtitle="Global")

plot.widetree.fdr.dec.ind <- ggplot(dec.fdr, aes(x=p, y=FDR, color = Method, group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "FDR")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          
    legend.title = element_text(size = 8),         
    axis.title = element_text(size = 8),           
    axis.text = element_text(size = 5),            
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="", subtitle="Decremental")

plot.widetree.fdr.inc.ind <- ggplot(inc.fdr, aes(x=p, y=FDR, color = Method, group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "FDR")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          
    legend.title = element_text(size = 8),         
    axis.title = element_text(size = 8),           
    axis.text = element_text(size = 5),            
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="", subtitle="Incremental")



# Power plots
plot.widetree.pow.glo.ind <- ggplot(glo.pow, aes(x=p, y=Power, color = Method, 
                                            group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "Power")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          
    legend.title = element_text(size = 8),         
    axis.title = element_text(size = 8),           
    axis.text = element_text(size = 5),            
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="Wide Tree", subtitle="Global")

plot.widetree.pow.dec.ind <- ggplot(dec.pow, aes(x=p, y=Power, color = Method, group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "Power")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          
    legend.title = element_text(size = 8),         
    axis.title = element_text(size = 8),           
    axis.text = element_text(size = 5),            
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="", subtitle="Decremental")

plot.widetree.pow.inc.ind <- ggplot(inc.pow, aes(x=p, y=Power, color = Method, group=interaction(Smooth, Method))) +
  geom_line(aes(linetype = Smooth))+labs( x="p", y = "Power")+ 
  ylim(0,1) + scale_colour_manual(values = colours)+
  scale_linetype_manual(name = "Smoothing", values=c(2,1))+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 8),          # Increase legend text size
    legend.title = element_text(size = 8),         # Increase legend title size
    axis.title = element_text(size = 8),           # Increase axis title size
    axis.text = element_text(size = 5),            # Increase axis tick labels size
    plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),            # Increase plot title size
    plot.title.position = "panel", plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5),
    panel.grid.minor = element_blank())+
  labs(title="", subtitle="Incremental")


