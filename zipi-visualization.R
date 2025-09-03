rm(list=ls())

setwd("D:/phD/16s/16s_files/zipi/gc/")



zi_pi<-read.csv(file = "ZiPi.csv")

head(zi_pi)


library(ggplot2)
library(ggrepel)

zi_pi$module <- factor(zi_pi$module, levels = c( "3", "2", "4"))  #factor


p2 <- ggplot(zi_pi, aes(p, z)) +
  geom_point(aes(color = module), alpha = 1, size = 4) +
  scale_color_manual(values = c("#F9A220"  ,  "#00AC11"  ,  "#CA272C"))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), 
        legend.key = element_blank(),
        legend.text = element_text(color="black",size=16,family = "serif",face = "bold"),
        legend.title = element_text(color="black",size=18,family = "serif",face = "bold"),
        axis.text = element_text(color="black",size=16,family = "serif",face = "bold"),
        axis.title = element_text(color="black",size=18,family = "serif",face = "bold")) +
  labs(x = 'Among-module connectivity (Pi)', y = 'Within-module connectivity (Zi)', color = 'Type') +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+ #Define the range and scale of the x-axis  
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1))+ #Define the range and scale of the y-axis
  geom_vline(xintercept = 0.62, linetype="dotted",size = 0.6) +
  geom_hline(yintercept = 2.5, linetype="dotted",size = 0.6) +
  geom_label_repel(data = zi_pi, aes(label = label),
                   size = 3,
                   fontface="bold",
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "black",
                   show.legend = FALSE,
                   max.overlaps = 10000)
p2
  
  
  
  
  
        
        