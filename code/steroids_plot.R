setwd("/Users/tili/Desktop/STAT_ovary/")
df <- read.csv("steroid.csv", header = TRUE, sep = ";", dec = ",")

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

df <- na.omit(df)

plot_list <- list()
for (i in 5:14) {
  temp <- df[,c(1:4,i)]
  temp$IVA <- factor(temp$IVA, levels = c("No", "Yes"))
  plot_list[[colnames(temp)[5]]] <- temp %>% ggplot(aes(x = IVA, y = log10(.[,5]+1), fill=IVA)) + geom_boxplot() +
    theme  + 
    ggtitle(colnames(temp)[5]) + ylab("Log10(Concentration+1)") + xlab("") + 
    guides(fill = guide_legend(title = "Groups"))  +
    scale_fill_manual(values = c("#8dd3c7","#fb8072"))
}
wrap_plots(plot_list, nrow = 2, ncol = 5) + plot_layout(guides = "collect")

for (i in 5:14) {
  temp <- df[,c(1:4,i)]
  temp$IVA <- factor(temp$IVA, levels = c("No", "Yes"))
  print(colnames(temp)[5])
  model <- lm(log10(temp[,5] + 1) ~ IVA + Age, data = temp)
  print(summary(model))
}

