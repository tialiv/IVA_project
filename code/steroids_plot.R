
setwd("/Users/tianyi.li.2/Desktop/STAT_ovary/data/")
#df <- read.csv("Steroids_CE groups.csv", header = TRUE, sep = ";", dec = ",") # Old data
df <- read.csv("/Users/tianyi.li.2/Desktop/STAT_ovary/data/E18_steroid_results.csv", 
               stringsAsFactors = F,
               header = TRUE, sep = ";", dec = ",")

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


##################### not used in the ms ########################
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
df$IVA <- ifelse(str_detect(df$IVA, "0"), "Fragmentation", "IVA")

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


########################### Used in the ms ######################
# not normalized to follicle number
plot_list <- list()
for (i in 6:15) {
  temp <- df[,c(5,i)]
  temp[is.infinite(temp[,2]),2] <- 0
  temp$IVA <- factor(temp$IVA, levels = c("Fragmentation", "IVA"))
  plot_list[[colnames(temp)[2]]] <- temp %>% ggplot(aes(x = IVA, y = log10(.[,2]+1), fill=IVA)) + geom_boxplot() + theme + 
    ggtitle(colnames(temp)[2]) + ylab("Log10(Concentration+1)") + xlab("") + 
    guides(fill = guide_legend(title = "Groups"))  + ylim(c(0, 2.5)) +
    scale_fill_manual(values = c("#8dd3c7", "#fb8072"))
}
wrap_plots(plot_list[c(10,7,8,5,6,4)], nrow = 2, ncol = 3) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

df <- df %>% mutate(gr = paste0(group, "_", Time))
df$gr <- factor(df$gr, levels = c("Fragmentation_Day_1", "IVA_Day_1",
                                  "Fragmentation_Day_3", "IVA_Day_3",
                                  "Fragmentation_Day_5", "IVA_Day_5",
                                  "Fragmentation_Day_7", "IVA_Day_7"))
p1 <- ggplot(df, aes(x = gr, y = log10(as.numeric(pregnenolone)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("Pregnenolone") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p2 <- ggplot(df, aes(x = gr, y = log10(as.numeric(progesterone)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("Progesterone") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p3 <- ggplot(df, aes(x = gr, y = log10(as.numeric(DHEA_normalized)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("DHEA") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p4 <- ggplot(df, aes(x = gr, y = log10(as.numeric(androstenedione)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("Androstenedione") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p5 <- ggplot(df, aes(x = gr, y = log10(as.numeric(T_normalized)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("Testosterone") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p6 <- ggplot(df, aes(x = gr, y = log10(as.numeric(E2_normalized)+1), fill=group)) + geom_boxplot() + theme  + 
  ggtitle("Estradiol") + ylab("Log10(Concentration+1)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

p1 + p2 + p3 + p4 + p5 + p6 + plot_annotation(tag_levels = "A") + plot_layout(ncol = 3, guides = "collect")


ggplot(df, aes(x = gr, y = (log10(as.numeric(E2_normalized)+1)/log10(as.numeric(T_normalized)+1)), fill=group)) + 
  geom_boxplot() + theme  + 
  ggtitle("Estradiol/Testosterone") + ylab("Percent (%)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))

ggplot(df, aes(x = gr, y = (log10(as.numeric(DHEA_normalized)+1)/log10(as.numeric(progesterone)+1)), fill=group)) + 
  geom_boxplot() + theme  + 
  ggtitle("DHEA/progesterone") + ylab("Percent (%)") + xlab("") + 
  guides(fill = guide_legend(title = "Groups")) + 
  scale_fill_manual(values = c("#8dd3c7", "#fb8072"))


################ Statistics ###############

# Average level and standard deviation of selected steroids
d5 <- df %>% select(gr, E2_normalized)
d5[,2] <- as.numeric(d5[, 2])
d5 %>% group_by(gr) %>% summarise(mean(E2_normalized))
d5 %>% group_by(gr) %>% summarise(sd(E2_normalized))

# One-way ANOVA
d5 <- df %>% select(gr, pregnenolone)
d5[,2] <- as.numeric(d5[, 2])
an <- aov(log10(d5$pregnenolone + 1) ~ d5$gr)  
TukeyHSD(an)  

# Two-way ANOVA with interaction
d5 <- df %>% select(group, Time, E2_normalized)
d5[,3] <- as.numeric(d5[, 3])
an <- aov(log10(d5$E2_normalized + 1) ~ d5$group + d5$Time + d5$group:d5$Time)  
TukeyHSD(an)  

###################################### Cytokine ##############################

d3 <- read.csv("/Users/tianyi.li.2/Desktop/STAT_ovary/data/cytokine.csv", sep = ";", stringsAsFactors = F)
d3$Value <- gsub(",", ".", d3$Value)
options(digits=6)
d3[,4] <- as.numeric(d3[, 4])

# Two-way ANOVA with interaction
d7 <- d3 %>% filter(Target == "IL2")
an <- aov(d7$Value ~ d7$Group + d7$Time + d7$Group:d7$Time)  
TukeyHSD(an) 

# Calculate mean values of the two technical replicates
d7 <- d7 %>% mutate(gr = paste0(Sample, "_", Group, "_", Time))
d8 <- d7 %>% group_by(gr) %>% summarise(mean(Value))
d8[, 3:5] <- str_split(d8$gr, "_", simplify = T)
colnames(d8)[c(2,4:5)] <- c("Value", "Group", "Time")

an <- aov(d7$Value ~ d7$Group + d7$Time + d7$Group:d7$Time)  
TukeyHSD(an)  

