# This code calculates optimal designs for discriminating competitive 
# and noncompetitive inhibition models using the Ds optimal method, sequentially.
# here the competitive model is the data generator


par(mfrow=c(1,1))
rm(list=ls())
set.seed(123456789)
library(BB)
library(mgcv)
library(stats)
library(plotfunctions)
library(numDeriv)
library(MASS) #generalized inverse

xstart = matrix(c(0.00, 0.00,
                  0.00, 30.00,
                  0.00, 60.00,
                  15.00, 0.00,
                  15.00, 30.00,
                  15.00, 60.00,
                  30.00, 0.00,
                  30.00, 30.00,
                  30.00, 60.00), nrow=9,ncol=2,byrow=TRUE)


xstart1=xstart[,1]
xstart2=xstart[,2]

n <- nrow(xstart)
wstart=c(rep((1/n),n))


#prior parameter values resulting from parameter est. in 120 real data

#prior values
th0start <- c(7.298, 4.386, 2.582)
s0start <- c(0.114, 0.233, 0.145)
th1start <- c(8.696, 8.066, 12.057)
s1start<- c(0.222, 0.488, 0.671)

thEstart <- c(7.4253, 4.6808, 3.0581,0.9636)
sEstart <- c(0.1298, 0.2724, 0.2815,0.0191)


TH0START<-c(7.298, 4.386, 2.582)
TH1START<- c(8.696, 8.066, 12.057)
SIG0<-0.1553
SIG1<-0.2272

THESTART <- c(7.4253, 4.6808, 3.0581,0.9636)
SIGE<-0.1526


x1g=seq(0,30,length=31)
x2g=seq(0,60,length=61)
grid=expand.grid(x1g,x2g)


iter=0
max.iter=500


colfunc <- colorRampPalette(c("lightblue","cornflowerblue", "darkblue"))
colorgrad <- colfunc(max.iter)

THETA.HAT=matrix(0,max.iter,4)
THETA.HATBBo=matrix(0,max.iter,4)
S1SEq= S0SEq = th1l=th0l=th0u=th1u=matrix(0,max.iter,3)
SESEq=matrix(0,max.iter,4)
num.max=numeric(max.iter)
norm.sum.crit.ds=sum.crit.ds=max.sen.crit.ds=sigmaEhat=sigma1hat=numeric(max.iter)

####################################################

# Model 0: competitive model as a function of parameters and design points
f0 <- function(S, I, th){
  V <- th[1]
  Km <- th[2]
  Kic <- th[3]
  V * S / (Km * (1 + I/Kic) + S)
}

# Model 1: non competitive model as a function of parameters and design points
f1 <- function(S, I, th){
  V <- th[1]
  Km <- th[2]
  Kin <- th[3]
  V * S / ((Km + S) * (1 + I/Kin))
}

# Model 2  (encompassing)
fE <- function(S, I, th){
  V <- th[1]
  Km <- th[2]
  Klm <- th[3]
  Lmbda <- th[4]
  V * S / ( Km *(1+I/Klm) + S * (1 + (1-Lmbda) * I/Klm))
}

#-------------------------------------------------------
SE1 <- 0
S11 <- 0
########################################################

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(xstart1, xstart2, type = "n",xlab="x1",ylab="x2",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5); grid()
points(xstart1, xstart2,col=colorgrad[1], pch=16,cex=0.8)

gradientLegend(valRange=c(1,max.iter),color = c("lightblue","cornflowerblue", "darkblue"),
               nCol = 3,inside = FALSE, pos=.825,dec = 0,n.seg=3,cex.main=2, cex.lab=1.5, cex.axis=1.5)
#text(x=31.5, y = 62, lables=expression("Tttitle"))


# data generator model (competitive model)
y <- f0(xstart1, xstart2, TH0START)+rnorm(nrow(xstart), sd = SIG0)
for(i in 1:length(y)){
  if(y[i]<0) y[i]<-0
}

basia <- data.frame(y,xstart)
#------------------------------------------------------------------
# iterates the whole sequential procedure to calculate optimal designs
while(iter< max.iter){
  iter=iter+1
  print(iter)
  #-----------------------------------------------------------------------------------------------------------------------
  #                                                        Theta.HAT
  #-----------------------------------------------------------------------------------------------------------------------  
  
  # parameter estimation: parameters are updated at each step
  mEpar.est<-nls(y~b1*basia[,2]/(b2*(1+basia[,3]/b3)+basia[,2]*(1+(1-b4)*basia[,3]/b3)),algorithm="port",lower=c(0.001,0.001,0.001,0.0001),upper=c(Inf,Inf,Inf,0.9636), start= list(b1=thEstart[1],b2=thEstart[2],b3=thEstart[3],b4=thEstart[4]), data=basia)
  THETA.HAT[iter,]=thEstart<-summary(mEpar.est)$parameters[,1]
  SESEq[iter,] = sEstart<-summary(mEpar.est)$parameters[,2]
  sigmaEhat[iter]<-sqrt(sum(residuals(mEpar.est)^2)/(nrow(basia)-4))
  
  
  #-----------------------------------------------------------------------------------------------------------------------
  #                                                        X.HAT
  #-----------------------------------------------------------------------------------------------------------------------  
  
  #-----------------------------------------------------
  
  # calculates the first derivative of the parameters of the encompassing model
  # it is a function of parameter estimates and design supports
  F.deriv <- function(thE,support){
    
    FE3.full <- matrix(0, nrow = nrow(support), ncol = 3)
    FE.full <- matrix(0, nrow = nrow(support), ncol = 4)
    
    x1 <- support[,1]
    x2 <- support[,2]
    
    
    a <- thE[1] * x1 
    b <- 1 + x2 / thE[3]
    cc <- 1 + (1-thE[4]) * x2 / thE[3]
    d <- (b * thE[2]) + (x1 * cc)
    FE.full[,1] <- x1 / d
    FE.full[,2] <- -a / d^2 * b
    FE.full[,3] <- a / d^2 * ((thE[2] + x1 * (1-thE[4])) * x2 / thE[3]^2)
    FE.full[,4] <- a / d^2 * (x1 * x2 / thE[3] )
    
    FE3.full <- FE.full[,1:3]
    
    return(list(FE3.full=FE3.full,FE.full=FE.full))
  }
  
  deriv.result <- F.deriv(thEstart,xstart)
  FE.full.mat <- deriv.result$FE.full
  FE3.full.mat <- deriv.result$FE3.full
  
  # full information matrix
  M.mat <- t(FE.full.mat) %*% FE.full.mat
  
  # information matrix of the nuisance parameters
  ME3.mat <- t(FE3.full.mat) %*% FE3.full.mat  
  #ME3.mat <- M.mat[1:3,1:3]
  
  # Ds variance function to calculate the next support point
  ds.funct <- function(xx1,xx2){
    
    f0.full <- matrix(0, nrow = 1, ncol = 3)
    f1.full <- matrix(0, nrow = 1, ncol = 3)
    fE3.full <- matrix(0, nrow = 1, ncol = 3)
    fE.full <- matrix(0, nrow = 1, ncol = 4)
    
    
    a <- thEstart[1] * xx1 
    b <- 1 + xx2 / thEstart[3]
    cc <- 1 + (1-thEstart[4]) * xx2 / thEstart[3]
    d <- (b * thEstart[2]) + (xx1 * cc)
    fE.full[1] <- xx1 / d
    fE.full[2] <- -a / d^2 * b
    fE.full[3] <- a / d^2 * ((thEstart[2] + xx1 * (1-thEstart[4])) * xx2 / thEstart[3]^2)
    fE.full[4] <- a / d^2 * (xx1 * xx2 / thEstart[3] )
    
    fE3.full <- matrix(fE.full[1,1:3],1,3)
    
    xx=matrix(c(xx1,xx2),1,2)
    
    d<-fE.full%*%ginv(M.mat)%*%t(fE.full)-fE3.full%*%ginv(ME3.mat)%*%t(fE3.full)
    d
  }
  d.value <- mapply(ds.funct,grid[,1],grid[,2])
  
  max.sen.crit.ds[iter]<-max(d.value)
  xnew<-c(grid[which.max(d.value),1],grid[which.max(d.value),2])
  
  #-----------------------------------------------------
  xstart<-rbind(xstart,xnew)
  
  deriv.result1 <- F.deriv(thEstart,xstart)
  FFE.full.mat <- deriv.result1$FE.full
  FFE3.full.mat <- deriv.result1$FE3.full
  
  
  MM.mat <- t(FFE.full.mat)  %*% FFE.full.mat
  MME3.mat <- t(FFE3.full.mat)  %*%  FFE3.full.mat
  #MME3.mat <- MM.mat[1:3,1:3]
  Cri.val <- det(MM.mat)/det(MME3.mat)
  sum.crit.ds[iter] <- Cri.val
  norm.sum.crit.ds[iter] <- Cri.val/nrow(xstart)
  
  
  ynew <- f0(xnew[1], xnew[2], TH0START)+rnorm(1, sd = SIG0)
  if(ynew<0) ynew<-0
  
  y <- c(y,ynew) 
  
  
  cexxx=(length(which(xstart[,1]==xnew[1]&xstart[,2]==xnew[2])))^(1/2)
  points(xnew[1], xnew[2], type = "p",col=colorgrad[iter], pch=16,cex=0.8*cexxx,cex.main=2, cex.lab=1.5, cex.axis=1.5 )
  
  basia<-data.frame(y,xstart)
  
  
} 
#----------------------------------------------------------------
design1.uniq=uniquecombs(xstart)
count1=attr(design1.uniq,"index")
count2=unique(sort(count1))
#----------------------------------------------------------------
counter=counter2=numeric(nrow(design1.uniq))
for(j in 1:length(count2)){
  m=which(count1==count2[j])
  #weight[j]=sum(design[m,3])
  counter[j]=length(m)
  counter2[j]=length(m)/max.iter
}
#----------------------------------------------------------------
design.opt=cbind(design1.uniq,counter,counter2)
design.optF <- as.data.frame(design.opt)
design.optF
#################################
cc1<-which(y==0)
sum(cc1) 

cc2<-numeric(length(THESTART))
for(k in 1:length(THESTART)){
  ccc<-which(THETA.HAT[,k]==0.001)
  cc2[k]<-sum(ccc) 
} 
cc2

##################################################################
quanDs <- c(0.25,0.5,0.75,0.9,0.95)
quanDs


earlyiterdsup <- numeric(length(quanDs))
for(i in 1:length(quanDs)){
  earlyiterdsup[i] <- which(norm.sum.crit.ds >= max(norm.sum.crit.ds)*quanDs[i])[1]
}
earlyiterdsup
#[1]  2  6 30 45 60
#----------------------------------------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/maxsen.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), max.sen.crit.ds , type = "p",xlab="iter",ylab="Maximum sensitivity",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/maxsen.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), max.sen.crit.ds , type = "p",xlab="iter",ylab="Maximum sensitivity",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#-------------------------------
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/criterion.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),sum.crit.ds , type = "p",xlab="iter",ylab="Criterion",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/criterion.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),sum.crit.ds , type = "p",xlab="iter",ylab="Criterion",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#-------------------------------
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/normcrit.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),norm.sum.crit.ds  , type = "p",ylim=c(0,4),xlab="iter",ylab="Norm.Criterion Ds",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=max(norm.sum.crit.ds ),lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/normcrit.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),norm.sum.crit.ds  , type = "p",ylim=c(0,4),xlab="iter",ylab="Norm.Criterion Ds",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=max(norm.sum.crit.ds ),lty=2,col="red")

dev.off()
#--------------------------------#-------------------------------

#------------------------------------------------

#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/sigmaE.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigmaEhat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[2])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG0,lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/sigmaE.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigmaEhat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[2])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG0,lty=2,col="red")

dev.off()
#--------------------------------#-------------------------------
#################################################################
#if we dont restrict lambda (in extimation part) we go far above the theoretical bound

SIG0
min(sigmaEhat)
max(sigmaEhat)
max(norm.sum.crit.ds)
max(norm.sum.crit.ds)/13

par(mfrow=c(1,1))


quantile(norm.sum.crit.ds)
per0.5 <- quantile(norm.sum.crit.ds)[3]-quantile(norm.sum.crit.ds)[1]
per0.5 / (max(norm.sum.crit.ds)-min(norm.sum.crit.ds))

which(norm.sum.crit.ds >= 3*norm.sum.crit.ds[1])[1]
length(which(norm.sum.crit.ds >= 3*norm.sum.crit.ds[1]))
boxplot(norm.sum.crit.ds)


#################################
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/sigmaparsE.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),SESEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,4] , type = "p",xlab="iter",pch=16,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#abline(h=SIG0,lty=2,col="red")

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep="")),
         expression(paste(hat(sigma[lambda]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen","black"),pch=16,bty="n",cex=1.5, pt.cex = 1)


dev.off()

which(SESEq[,1] <= SIGE)[1]
which(SESEq[,2] <= SIGE)[1]
which(SESEq[,3] <= SIGE)[1]
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/sigmaparsE.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),SESEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),SESEq[,4] , type = "p",xlab="iter",pch=16,col="black",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#abline(h=SIG0,lty=2,col="red")

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep="")),
         expression(paste(hat(sigma[lambda]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen","black"),pch=16,bty="n",cex=1.5, pt.cex = 1)


dev.off()
#--------------------------------#-------------------------------

#################################
TH0START
min(THETA.HAT[,1])
max(THETA.HAT[,1])

min(THETA.HAT[,2])
max(THETA.HAT[,2])

min(THETA.HAT[,3])
max(THETA.HAT[,3])

min(THETA.HAT[,4])
max(THETA.HAT[,4])

min(THETA.HAT)
max(THETA.HAT)
TH0START

#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th21.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH0START[1],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th21.eps")


par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH0START[1],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th22.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH0START[2],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th22.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH0START[2],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th23.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH0START[3],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th23.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH0START[3],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th24.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,4] ,ylim=c(0,1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(lambda)),pch=16,col="blue")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th24.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,4] ,ylim=c(0,1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(lambda)),pch=16,col="blue")


dev.off()
#--------------------------------
#################################

#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th21s.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=THESTART[1],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th21s.eps")


par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=THESTART[1],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th22s.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=THESTART[2],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th22s.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=THESTART[2],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th23s.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=THESTART[3],lty=2,col="red")


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Dsf0/think-Ds/th23s.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=THESTART[3],lty=2,col="red")


dev.off()
#--------------------------------



