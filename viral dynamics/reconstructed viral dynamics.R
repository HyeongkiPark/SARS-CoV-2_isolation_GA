CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)
rm(list=ls())

library(plyr)
library(ggplot2)
library(gridExtra)
library(FME)
library(data.table) # fread

##################################################################
### Functions
##################################################################
Tmin <- 0.0
Tmax <- 40
step_size <- 0.01
stime <- seq(Tmin,Tmax,step_size)

Covfun <- function(pars){
  
  beta <- as.numeric(pars["beta1"]) * (10^(-5))
  r <- as.numeric(pars["r"])
  delta <- as.numeric(pars["delta"])
  
  derivs<-function(time, y, pars){
    with(as.list(c(pars, y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  y<-c(Ta=1, V=0.01)
  
  times<-c(seq(0,40,0.01))
  out<-lsoda(y=y, parms=pars, times=times, func=derivs, rtol=0.00004, atol=0.00000000000001)
  out2<-cbind(time=out[,1],aV=((log10(out[,3]))))
  as.data.frame(out2)
}

sensrfun2 <- function(pars) {
  return( Covfun(pars) )
}

##################################################################
### Data: 3VOCs
##################################################################
PopPars_3VOC <- fread("monolix_estimated_parameter/3VOC/populationParameters.txt")
IndPars_3VOC <- fread("monolix_estimated_parameter/3VOC/estimatedIndividualParameters.txt")

pars_3VOC_0 <- c(beta1 = PopPars_3VOC$value[7], r = PopPars_3VOC$value[1], delta=PopPars_3VOC$value[4])
pars_3VOC_1 <- c(beta1 = PopPars_3VOC$value[7]*exp(PopPars_3VOC$value[8]),
            r = PopPars_3VOC$value[1]*exp(PopPars_3VOC$value[2]),
            delta=PopPars_3VOC$value[4]*exp(PopPars_3VOC$value[5]))
pars_3VOC_2 <- c(beta1 = PopPars_3VOC$value[7]*exp(PopPars_3VOC$value[9]),
            r = PopPars_3VOC$value[1]*exp(PopPars_3VOC$value[3]),
            delta=PopPars_3VOC$value[4]*exp(PopPars_3VOC$value[6]))
tau_3VOC = PopPars_3VOC$value[10]
tau1_3VOC = PopPars_3VOC$value[10]*exp(PopPars_3VOC$value[11])
tau2_3VOC = PopPars_3VOC$value[10]*exp(PopPars_3VOC$value[12])

IndPars_3VOC <- IndPars_3VOC[,c(12,10,11,18)]
names(IndPars_3VOC) <- c("beta1", "r", "delta", "variant")

IndPars_3VOC_0 <- subset(IndPars_3VOC, variant==0)[,c(1:3)]
IndPars_3VOC_1 <- subset(IndPars_3VOC, variant==1)[,c(1:3)]
IndPars_3VOC_2 <- subset(IndPars_3VOC, variant==2)[,c(1:3)]


##################################################################
### Data: 3VOCs with immunity
##################################################################
PopPars_3VOC_im <- fread("monolix_estimated_parameter/3VOC_with_prior_immunity/populationParameters.txt")
IndPars_3VOC_im <- fread("monolix_estimated_parameter/3VOC_with_prior_immunity/estimatedIndividualParameters.txt")

pars_3VOC_im_0 <- c(beta1 = PopPars_3VOC_im$value[7], r = PopPars_3VOC_im$value[1], delta=PopPars_3VOC_im$value[4])
pars_3VOC_im_1 <- c(beta1 = PopPars_3VOC_im$value[7]*exp(PopPars_3VOC_im$value[8]),
                 r = PopPars_3VOC_im$value[1]*exp(PopPars_3VOC_im$value[2]),
                 delta=PopPars_3VOC_im$value[4]*exp(PopPars_3VOC_im$value[5]))
pars_3VOC_im_2 <- c(beta1 = PopPars_3VOC_im$value[7]*exp(PopPars_3VOC_im$value[9]),
                 r = PopPars_3VOC_im$value[1]*exp(PopPars_3VOC_im$value[3]),
                 delta=PopPars_3VOC_im$value[4]*exp(PopPars_3VOC_im$value[6]))
tau_3VOC_im = PopPars_3VOC_im$value[10]
tau1_3VOC_im = PopPars_3VOC_im$value[10]*exp(PopPars_3VOC_im$value[11])
tau2_3VOC_im = PopPars_3VOC_im$value[10]*exp(PopPars_3VOC_im$value[12])

IndPars_3VOC_im <- IndPars_3VOC_im[,c(12,10,11,18)]
names(IndPars_3VOC_im) <- c("beta1", "r", "delta", "variant")

IndPars_3VOC_im_0 <- subset(IndPars_3VOC_im, variant==0)[,c(1:3)]
IndPars_3VOC_im_1 <- subset(IndPars_3VOC_im, variant==1)[,c(1:3)]
IndPars_3VOC_im_2 <- subset(IndPars_3VOC_im, variant==2)[,c(1:3)]


######################
### Data: omi, delta
######################
PopPars <- fread("monolix_estimated_parameter/del_omi/populationParameters.txt")
IndPars <- fread("monolix_estimated_parameter/del_omi/estimatedIndividualParameters.txt")

pars_del <- c(beta1 = PopPars$value[5], r = PopPars$value[1], delta=PopPars$value[3])
pars_omi <- c(beta1 = PopPars$value[5]*exp(PopPars$value[6]),
              r = PopPars$value[1]*exp(PopPars$value[2]),
              delta=PopPars$value[3]*exp(PopPars$value[4]))

tau <- PopPars$value[7]
tau1 <- PopPars$value[7]*exp(PopPars$value[8])


IndPars <- IndPars[,c(12,10,11,18)]
names(IndPars) <- c("beta1", "r", "delta", "variant")

IndPars_del <- subset(IndPars, variant==2)[,c(1:3)]
IndPars_omi <- subset(IndPars, variant==3)[,c(1:3)]


##################################################################
##### Run plot code: 3VOCs with/without immunity (Figure1, S7)
##################################################################
qcolors <- c("black","#00a0e9","#f65294","#009944")

## For Figure 1
pop_pars0 <- pars_3VOC_0
pop_pars1 <- pars_3VOC_1
pop_pars2 <- pars_3VOC_2
ind_pars0 <- IndPars_3VOC_0
ind_pars1 <- IndPars_3VOC_1
ind_pars2 <- IndPars_3VOC_2

## For Figure S7
# pop_pars0 <- pars_3VOC_0_im
# pop_pars1 <- pars_3VOC_1_im
# pop_pars2 <- pars_3VOC_2_im
# ind_pars0 <- IndPars_3VOC_0_im
# ind_pars1 <- IndPars_3VOC_1_im
# ind_pars2 <- IndPars_3VOC_2_im

DL <- 0
sR0 <- sensRange(func=sensrfun2, parms=pop_pars0, parInput=ind_pars0)
sR0 <- sR0[complete.cases(sR0), ]
name_lst <- colnames(sR0)
gids <- grep("aV",name_lst)
mat <- sR0[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<DL] <- DL
ylow[ylow< DL] <- DL
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted <- Covfun(pop_pars0)
d0 <- data.frame(x=times,y=(fitted$aV))
dc0 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

sR1 <- sensRange(func=sensrfun2, parms=pop_pars1, parInput=ind_pars1)
sR1 <- sR1[complete.cases(sR1), ]
name_lst <- colnames(sR1)
gids <- grep("aV",name_lst)
mat <- sR1[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<DL] <- DL
ylow[ylow< DL] <- DL
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted <- Covfun(pop_pars1)
d1 <- data.frame(x=times,y=(fitted$aV))
dc1 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)


sR2 <- sensRange(func=sensrfun2, parms=pop_pars2, parInput=ind_pars2)
sR2 <- sR2[complete.cases(sR2), ]
name_lst <- colnames(sR2)
gids <- grep("aV",name_lst)
mat <- sR2[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<DL] <- DL
ylow[ylow< DL] <- DL
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted <- Covfun(pop_pars2)
d2 <- data.frame(x=times,y=(fitted$aV))
dc2 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

xlabel <- "Days after infection"
ylabel <- "Viral RNA load (copies/ml)"
plt <- ggplot() + 
  geom_ribbon(data=dc0,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[1],alpha=0.2) + # alpha=0.2: NOT supported in EPS
  geom_ribbon(data=dc1,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[2],alpha=0.2) + 
  geom_ribbon(data=dc2,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[3],alpha=0.2) + 
  geom_path(data=d0,aes(x=x,y=y),color=qcolors[1],lwd=2.5) +
  geom_path(data=d1,aes(x=x,y=y),color=qcolors[2],lwd=2.5) +
  geom_path(data=d2,aes(x=x,y=y),color=qcolors[3],lwd=2.5) +
  xlab(xlabel) + ylab(ylabel) +
  scale_x_continuous(breaks=c(seq(0,40,by=4)), labels = as.character(c(seq(0,40,by=4))),limits=c(0,40)) +
  scale_y_continuous(breaks=seq(0,10,by=1),labels = expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10), limit=c(0,10)) +
  theme(axis.text = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position='none',
        axis.title.y = element_text(family="Helvetica"),axis.title.x = element_text(family="Helvetica"))

plt


##################################################################
##### Run plot code: delta, omicron (Figure 4)
##################################################################
qcolors2 <- c("#f65294","#009944")

pop_pars0 <- pars_del
pop_pars1 <- pars_omi

ind_pars0 <- IndPars_del
ind_pars1 <- IndPars_omi

DL <- 0
sR0 <- sensRange(func=sensrfun2, parms=pop_pars0, parInput=ind_pars0)
sR0 <- sR0[complete.cases(sR0), ]
name_lst <- colnames(sR0)
gids <- grep("aV",name_lst)
mat <- sR0[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<DL] <- DL
ylow[ylow< DL] <- DL
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted <- Covfun(pop_pars0)
d0 <- data.frame(x=times,y=(fitted$aV))
dc0 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

sR1 <- sensRange(func=sensrfun2, parms=pop_pars1, parInput=ind_pars1)
sR1 <- sR1[complete.cases(sR1), ]
name_lst <- colnames(sR1)
gids <- grep("aV",name_lst)
mat <- sR1[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<DL] <- DL
ylow[ylow< DL] <- DL
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted <- Covfun(pop_pars1)
d1 <- data.frame(x=times,y=(fitted$aV))
dc1 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

xlabel <- "Days after infection"
ylabel <- "Viral RNA load (copies/ml)"
plt <- ggplot() + 
  geom_ribbon(data=dc0,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors2[1],alpha=0.2) + # alpha=0.2: NOT supported in EPS
  geom_ribbon(data=dc1,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors2[2],alpha=0.2) + 
  geom_path(data=d0,aes(x=x,y=y),color=qcolors2[1],lwd=2.5) +
  geom_path(data=d1,aes(x=x,y=y),color=qcolors2[2],lwd=2.5) +
  xlab(xlabel) + ylab(ylabel) +
  scale_x_continuous(breaks=c(seq(0,40,by=4)), labels = as.character(c(seq(0,40,by=4))),limits=c(0,40)) +
  scale_y_continuous(breaks=seq(0,10,by=1),labels = expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10), limit=c(0,10)) +
  theme(axis.text = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position='none',
        axis.title.y = element_text(family="Helvetica"),axis.title.x = element_text(family="Helvetica"))

plt

##################################################################
##### Run plot code: 3VOCs individual plot (Figure S1)
##################################################################
qcolors <- c("black","#00a0e9","#f65294","#009944")
IndiPars <- fread("monolix_estimated_parameter/3VOC/estimatedIndividualParameters.txt")

plot_list <- vector(mode="list",length=length(IndiPars$id))
xlabel <- "Days since SARS-CoV-2 peak"
ylabel <- "Viral RNA load\n(copies/mL)"
for(i in 1:length(IndiPars$id)){
  infectiontime <- IndiPars$tau_mode[i]
  newtimes <- c(seq(-infectiontime, 40-infectiontime, 0.01))
  sR_pars <- c(beta = IndiPars$beta1_mode[i]*10^(-5), r = IndiPars$r_mode[i], delta = IndiPars$delta_mode[i])
  fitted <- Covfun(sR_pars)
  fitted1 <- data.frame(time = newtimes, aV = fitted$aV)
  id_number <- IndiPars$id[i]
  curvecolor <- qcolors[as.numeric(IndiPars$variant[i]) + 1]
  
  plt <- ggplot() +
    geom_line(data = fitted1, aes(x=time, y=aV), lwd=1, color = curvecolor) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(breaks=seq(-20,40,by=10),labels = expression(-20,-10,0,10,20,30,40),limits=c(-20,40)) +
    scale_y_continuous(breaks=seq(1,9,by=2),labels = expression(10^1,10^3,10^5,10^7,10^9),limits=c(0,10)) +
    theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position='none', axis.title.y = element_text(size=7,family="Helvetica"), axis.title.x = element_text(size=9,family="Helvetica"))
  
  plot_list[[i]] <- plt
}

plot_list[[3]]

##################################################################
##### Run plot code: delta and omicron individual plot (Figure S6)
##################################################################
qcolors <- c("#f65294","#009944")
IndiPars <- fread("monolix_estimated_parameter/del_omi/estimatedIndividualParameters.txt")

plot_list <- vector(mode="list",length=length(IndiPars$id))
xlabel <- "Days since SARS-CoV-2 peak"
ylabel <- "Viral RNA load\n(copies/mL)"
for(i in 1:length(IndiPars$id)){
  infectiontime <- IndiPars$tau_mode[i]
  newtimes <- c(seq(-infectiontime, 40-infectiontime, 0.01))
  sR_pars <- c(beta = IndiPars$beta1_mode[i]*10^(-5), r = IndiPars$r_mode[i], delta = IndiPars$delta_mode[i])
  fitted <- Covfun(sR_pars)
  fitted1 <- data.frame(time = newtimes, aV = fitted$aV)
  id_number <- IndiPars$id[i]
  curvecolor <- qcolors[as.numeric(IndiPars$variant[i]) - 1]
  
  plt <- ggplot() +
    geom_line(data = fitted1, aes(x=time, y=aV), lwd=1, color = curvecolor) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(breaks=seq(-20,40,by=10),labels = expression(-20,-10,0,10,20,30,40),limits=c(-20,40)) +
    scale_y_continuous(breaks=seq(1,9,by=2),labels = expression(10^1,10^3,10^5,10^7,10^9),limits=c(0,10)) +
    theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position='none', axis.title.y = element_text(size=7,family="Helvetica"), axis.title.x = element_text(size=9,family="Helvetica"))
  
  plot_list[[i]] <- plt
}

plot_list[[3]]

