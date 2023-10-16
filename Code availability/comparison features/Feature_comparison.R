CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)
rm(list=ls())

library(plyr)
library(ggplot2)
library(gridExtra)
library(FME)
library(data.table) # fread
library(multcomp)

## data
df <- fread("Features_delomi.csv")
df_2 <- subset(df, variant==2)
df_3 <- subset(df, variant==3)

###################################################################
### Distribution of features
###################################################################
qcolors <- c("black","#00a0e9","#f65294","#009944")
vlinesize <- 2
vlinealpha <- 0.5

Features<-df
df2<-df_2
df3<-df_3

xlabel<-"Duration"
ylabel<-"Density"
plt_dur <- ggplot()+
  geom_density(data=Features, aes(x=duration, alpha=0.5, fill=as.factor(variant))) +
  scale_fill_manual(values = qcolors[3:4]) +
  geom_vline(xintercept = mean(df2$duration), col=qcolors[3], size=vlinesize, alpha=vlinealpha) +
  geom_vline(xintercept = mean(df3$duration), col=qcolors[4], size=vlinesize, alpha=vlinealpha) +
  xlab(xlabel) + ylab(ylabel)+
  scale_x_continuous(expand = c(0, 0), breaks=seq(5,30,by=5), limits=c(4.5,30.5))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,0.2,by=0.1), limits=c(0,0.21))+
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none')

plt_dur

xlabel <- "Peak viral load"
ylabel <- "Density"
plt_peak <- ggplot()+
  geom_density(data=Features, aes(x=peak, alpha=0.5, fill=as.factor(variant))) +
  scale_fill_manual(values = qcolors[3:4]) +
  geom_vline(xintercept = mean(df2$peak), col=qcolors[3], size=vlinesize, alpha=vlinealpha) +
  geom_vline(xintercept = mean(df3$peak), col=qcolors[4], size=vlinesize, alpha=vlinealpha) +
  xlab(xlabel) + ylab(ylabel)+
  scale_x_continuous(expand = c(0, 0), breaks=seq(4,10,by=1), limits=c(3.8,10.1))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,0.6,by=0.3), limits=c(0,0.65))+
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none') 
plt_peak

xlabel <- "Peak day(day)"
ylabel <- "Density"
plt_pday <- ggplot()+
  geom_density(data=Features, aes(x=peak_day, alpha=0.5, fill=as.factor(variant))) +
  scale_fill_manual(values = qcolors[3:4]) +
  geom_vline(xintercept = mean(df2$peak_day), col=qcolors[3], size=vlinesize, alpha=vlinealpha) +
  geom_vline(xintercept = mean(df3$peak_day), col=qcolors[4], size=vlinesize, alpha=vlinealpha) +
  xlab(xlabel) + ylab(ylabel)+
  scale_x_continuous(expand = c(0, 0), breaks=seq(2,6,by=1), limits=c(1.8,6.2))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,1.00,by=0.5), limits=c(0,1.2))+
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none')

plt_pday
