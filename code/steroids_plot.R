setwd("/Users/tili/Desktop/STAT_ovary/data/")
df <- read.csv("Steroids_CE groups.csv", header = TRUE, sep = ";", dec = ",")

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggbreak)

theme <- theme(text=element_text(size=15),
               axis.text.x=element_text(color="black", angle = 90, vjust = 0.3, hjust = 0.3),
               axis.text.y=element_text(color="black"),
               axis.ticks.x.bottom=element_blank(),
               #axis.ticks.length.x=unit(-.25, "cm"),
               strip.background = element_rect(colour=NA, fill=NA),
               plot.background = element_blank(),
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(),
               panel.spacing.x=unit(.5, "lines"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               axis.line.x = element_line(size=0.4),
               axis.line.y = element_line(size=0.4),
               legend.key.size = unit(2, 'lines'),
               panel.background = element_blank()
)

df <- df[df$Number.of.follicles != 0,]
normalized <- data.frame(df$Biopsy.code) 

for (i in 6:15) {
  temp <- df[,c(3, i)]
  temp <- transform(temp, nor = temp[,2] / temp[,1]) 
  normalized <- cbind(normalized, temp$nor)
}

colnames(normalized) <- colnames(df)[c(2, 6:15)]
normalized$IVA <- df$IVA 
normalized$IVAs <- ifelse(str_detect(normalized$IVA, "0"), "No", "Yes")

plot_list <- list()
for (i in 2:11) {
  temp <- normalized[,c(13,i)]
  temp[is.infinite(temp[,2]),2] <- 0
  temp$IVAs <- factor(temp$IVAs, levels = c("No", "Yes"))
  plot_list[[colnames(temp)[2]]] <- temp %>% ggplot(aes(x = IVAs, y = log10(.[,2]+1), fill=IVAs)) + geom_boxplot() + theme  + 
    ggtitle(colnames(temp)[2]) + ylab("Log10(Concentration+1)") + xlab("") + 
    guides(fill = guide_legend(title = "Groups"))  +
    scale_fill_manual(values = c("#fb8072", "#8dd3c7"))
}
wrap_plots(plot_list[c(10,7,8,5,6,4)], nrow = 2, ncol = 3) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

for (i in 2:11) {
  temp <- normalized[,c(13,i)]
  temp[is.infinite(temp[,2]),2] <- 0
  temp$IVAs <- factor(temp$IVAs, levels = c("No", "Yes"))
  print(colnames(temp)[2])
  print(t.test(log10(temp[,2] + 1) ~ temp$IVAs))
}

