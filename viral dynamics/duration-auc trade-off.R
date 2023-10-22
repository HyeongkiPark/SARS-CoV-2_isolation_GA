CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)
rm(list=ls())

## Data
aucdata<-read.table("Features_3VOCs.csv", sep = ",", comment.char = "", header = T)
aucdata$auc <- round(aucdata$auc, 3)
fit_mixture <- c(198.58, 1.169)

## Function
ODEs<- function(pars) {
  D50 <- as.numeric(pars[1])
  Dk <- as.numeric(pars[2])
  Dmax <- 80
  aucscale <- seq(0,150,0.001)
 
  out <- Dmax*aucscale^Dk/(aucscale^Dk+D50^Dk)
  out2<-cbind(auc=aucscale,duration=out)
  as.data.frame(out2)
}
fitted <- ODEs(pars=fit_mixture)

dd_Ta<-data.frame(time=aucdata$auc, value=aucdata$duration, variant=aucdata$variant)
ds_Ta <- data.frame(time=fitted$auc, value=fitted$duration)

## Plot
xlabel <- "Cumulative log-trasformed viral load"
ylabel <- "Duration"
plt <- ggplot() + 
  geom_point(data=dd_Ta,aes(x=time,y=value,colour=as.factor(variant)),size=4,shape=19)+
  geom_path(data=ds_Ta,aes(x=time,y=value),color="gray",lwd=1) +
  scale_colour_manual(values=c('black',"#00a0e9","#f65294"))+
  xlab(xlabel) + ylab(ylabel) + 
  scale_y_continuous(breaks=seq(0,40,by=5), limits=c(0,40))+
  scale_x_continuous(breaks=seq(10,150,by=20),limits=c(10,150))+
  theme(axis.text = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position='none',axis.title.y = element_text(family="Helvetica"),axis.title.x = element_text(family="Helvetica"))

plt
