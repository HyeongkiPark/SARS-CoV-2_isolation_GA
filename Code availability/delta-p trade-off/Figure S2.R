CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)
rm(list=ls())

## data
fitall <- read.table("Features_3VOCs.csv", sep = ",", comment.char = "", header = T)
fit_mixture <-c(D50=100, Dmax=1.57)

prealpha<-subset(fitall, fitall$variant == 0)
alpha<-subset(fitall, fitall$variant == 1)
delta<-subset(fitall, fitall$variant == 2)

fitall<-rbind(prealpha,alpha,delta)

## function
ODEs<- function(pars) {
  D50 <- as.numeric(pars[1])
  Dmax <- as.numeric(pars[2])
  aucscale <- seq(0,650000,1)
  
  out <- Dmax*aucscale/(aucscale+D50)
  out2<-cbind(auc=aucscale, duration=out)
  as.data.frame(out2)
}

fitted <- ODEs(pars=fit_mixture)

## Plot
dd_Ta <- data.frame(time=log10(fitall$p), value=fitall$delta, variant=fitall$variant)
ds_Ta <- data.frame(time=log10(fitted$auc), value=fitted$duration)

xlabel <- expression(p ~ '(log10)')
ylabel <- expression(delta)
plt <- ggplot() + 
  geom_path(data=ds_Ta,aes(x=time,y=value),color="black",lwd=1) +
  geom_point(data=dd_Ta,aes(x=time,y=value,colour=as.factor(variant)),size=3,shape=19)+
  scale_colour_manual(values=c('black',"#00a0e9","#f65294"))+
  xlab(xlabel) + ylab(ylabel) + 
  scale_x_continuous(breaks = scales::breaks_extended(n = 10))+
  scale_y_continuous(breaks = scales::breaks_extended(n = 10)) +
  theme(axis.text = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position='none',axis.title.y = element_text(family="Helvetica"),axis.title.x = element_text(family="Helvetica"))

plt

