CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)
rm(list=ls())

## Data
cdfpar <- fread("monolix_estimated_parameter/3VOC/populationParameters.txt")
fitall <- fread("monolix_estimated_parameter/3VOC/estimatedIndividualParameters.txt")

prealpha<-subset(fitall, fitall$variant == 0)
alpha<-subset(fitall, fitall$variant == 1)
delta<-subset(fitall, fitall$variant == 2)
fitall<-rbind(prealpha,alpha,delta)

Tmin <- 0.0
Tmax <- 40
step_size <- 0.01
stime <- seq(Tmin,Tmax,step_size)

## Function
Covfun<-function(pars){
  beta <- as.numeric(pars[1])
  r <- as.numeric(pars[2])
  delta <- as.numeric(pars[3])
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  y<-c(Ta=1,V=0.01)
  
  times<-c(seq(0,40,0.01))
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0.00000000000001)
  pro<-vector(mode="list",length=length(times))
  # for (i in 1:length(times)){
  #   out[i,3]=max(out[i,3],100)
  # }
  # 
  for (i in 1:length(times)){
    viral<-log10(out[i,3])
    #viral<-(out[i,3])
    if (viral <5) {
      pro[i]<-0
    } else if  (viral>=5 & viral <7) {
      pro[i]<-0.12
    } else if (viral>=7 & viral <10) {
      pro[i]<-0.15
    } else if (viral>=10) {
      pro[i]<-0.24
    }
  }
  out2<-cbind(time=out[,1],aV=pro,reV=((log10(out[,3]))))
  as.data.frame(out2)
}

## Calculate transmission probability and make data file
prealphapara<- c(beta=cdfpar$value[7]*(10^(-5)),r=cdfpar$value[1],delta=cdfpar$value[4])
alphapara<- c(beta=cdfpar$value[7]*exp(cdfpar$value[8])*(10^(-5)),r=cdfpar$value[1]*exp(cdfpar$value[2]),delta=cdfpar$value[4]*exp(cdfpar$value[5]))
deltapara<- c(beta=cdfpar$value[7]*exp(cdfpar$value[9])*(10^(-5)),r=cdfpar$value[1]*exp(cdfpar$value[3]),delta=cdfpar$value[4]*exp(cdfpar$value[6]))

fittedprealpha<-Covfun(prealphapara)
fittedalpha<-Covfun(alphapara)
fitteddelta<-Covfun(deltapara)

aaa1<-fitall[1:86,]
aaa2<-fitall[87:145,]
aaa3<-fitall[146:225,]

ds_prealpha<-data.frame(time=as.numeric(fittedprealpha$time),value=as.numeric(fittedprealpha$aV))
ds_alpha<-data.frame(time=as.numeric(fittedalpha$time),value=as.numeric(fittedalpha$aV))
ds_delta<-data.frame(time=as.numeric(fitteddelta$time),value=as.numeric(fitteddelta$aV))


p1<-subset(ds_prealpha,ds_prealpha$value==max(ds_prealpha$value))$time
p2<-subset(ds_alpha,ds_alpha$value==max(ds_alpha$value))$time
p3<-subset(ds_delta,ds_delta$value==max(ds_delta$value))$time


assign(paste0('parset',1),data.frame(cbind(beta=aaa1$beta_mode[1:length(aaa1$id)],r=aaa1$r_mode[1:length(aaa1$id)],delta=aaa1$delta_mode[1:length(aaa1$id)])))
assign(paste0('parset',2),data.frame(cbind(beta=aaa2$beta_mode[1:length(aaa2$id)],r=aaa2$r_mode[1:length(aaa2$id)],delta=aaa2$delta_mode[1:length(aaa2$id)])))
assign(paste0('parset',3),data.frame(cbind(beta=aaa3$beta_mode[1:length(aaa3$id)],r=aaa3$r_mode[1:length(aaa3$id)],delta=aaa3$delta_mode[1:length(aaa3$id)])))

sensrfun2 <- function(pars) {
  return( Covfun(pars) )
}

sR_pars1 <- prealphapara
sR_pars2 <- alphapara
sR_pars3 <- deltapara
allpar1<-parset1
allpar2<-parset2
allpar3<-parset3
sR1 <- sensRange(func=sensrfun2,parms=sR_pars1,parInput=allpar1)
sR2 <- sensRange(func=sensrfun2,parms=sR_pars2,parInput=allpar2)
sR3 <- sensRange(func=sensrfun2,parms=sR_pars3,parInput=allpar3)

sR1<-sR1[complete.cases(sR1), ]
name_lst <- colnames(sR1)
gids <- grep("aV",name_lst)
mat <- sR1[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,1)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<0] <- 0
ylow[ylow< 0] <- 0
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted<-Covfun(sR_pars1)
d1<-data.frame(x=times,y=(fitted$aV))
dc1 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

sR2<-sR2[complete.cases(sR2), ]
name_lst <- colnames(sR2)
gids <- grep("aV",name_lst)
mat <- sR2[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<0] <- 0
ylow[ylow< 0] <- 0
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted<-Covfun(sR_pars2)
d2<-data.frame(x=times,y=(fitted$aV))
dc2 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

sR3<-sR3[complete.cases(sR3), ]
name_lst <- colnames(sR3)
gids <- grep("aV",name_lst)
mat <- sR3[,gids]
times <- seq(Tmin,Tmax,length=length(gids))
ymean <- ( apply(mat,2,mean) )
yCIlow <- ( apply(mat,2,function(x){quantile(x,0.025)}) )
yCIhigh <- ( apply(mat,2,function(x){quantile(x,0.975)}) )
ylow <- as.numeric(yCIlow)
yup  <- as.numeric(yCIhigh)
yup[yup<0] <- 0
ylow[ylow< 0] <- 0
xrange <- c(Tmin,Tmax)
yrange <- c(min(ylow),max(yup))
xn <- length(times)
yn <- length(ymean)

labels <- gl(3,xn,label=c("mean","5%","95%"))
labels2 <- gl(1,xn,label="mean")
x <- rep(times,3)
fitted<-Covfun(sR_pars3)
d3<-data.frame(x=times,y=(fitted$aV))
dc3 <- data.frame(x=times,mean=ymean,ymin=ylow,ymax=yup)

## Plot
qcolors <- c("black","#00a0e9","#f65294")
xlabel <- "Days after infection"
ylabel <- "Transmission probability"
plt <- ggplot() + 
  geom_ribbon(data=dc1,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[1],alpha=0.2) + # alpha=0.2: NOT supported in EPS
  geom_ribbon(data=dc2,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[2],alpha=0.2) + 
  geom_ribbon(data=dc3,aes(x=x,ymin=ymin,ymax=ymax),fill=qcolors[3],alpha=0.2) + 
  geom_path(data=ds_alpha,aes(x=time,y=value),color=qcolors[2],lwd=3) +
  geom_path(data=ds_delta,aes(x=time,y=value),color=qcolors[3],lwd=3) +
  geom_path(data=ds_prealpha,aes(x=time,y=value),color=qcolors[1],lwd=3,linetype="dashed") +
  xlab(xlabel) + ylab(ylabel) +
  scale_x_continuous(breaks=c(seq(0,20,by=4)),labels = as.character(c(seq(0,20,by=4))),limits=c(0,20)) +
  scale_y_continuous(breaks=seq(0,0.2,by=0.05),labels = expression(0,0.05,0.1,0.15,0.2),limits=c(0,0.2)) +
  theme(axis.text = element_text(colour = "black"),axis.ticks = element_line(colour = "black"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position='none',axis.title.y = element_text(family="Helvetica"),axis.title.x = element_text(family="Helvetica"))

plt



