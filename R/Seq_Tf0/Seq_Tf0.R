# This code calculates optimal designs for discriminating competitive 
# and noncompetitive inhibition models using the Hunter and Reiner (1965) optimal method, sequentially.
# here the competitive model is the data generator

rm(list=ls())
par(mfrow=c(1,1))
set.seed(123456789)
library(BB)
library(mgcv)
library(bvls)
library(stats)
library(plotfunctions)

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
th0start <- c(7.298, 4.386, 2.582)
s0start <- c(0.114, 0.233, 0.145)
th1start <- c(8.696, 8.066, 12.057)
s1start<- c(0.222, 0.488, 0.671)

TH0START<-c(7.298, 4.386, 2.582)
TH1START<- c(8.696, 8.066, 12.057)
SIG0<-0.1553
SIG1<-0.2272
lambda0<-1
lambda1<-0

x1g=seq(0,30,length=31)
x2g=seq(0,60,length=61)
grid=expand.grid(x1g,x2g)


iter=0
max.iter=500

colfunc <- colorRampPalette(c("lightblue","cornflowerblue", "darkblue"))
colorgrad <- colfunc(max.iter)

THETA.HAT=matrix(0,max.iter,6)
S1SEq= S0SEq =th1l=th0l=th0u=th1u=matrix(0,max.iter,3)
num.max=numeric(max.iter)
norm.sum.crit.t=sum.crit.t=max.sen.crit.t=sigma0hat=sigma1hat=cexxx=numeric(max.iter)

# Model 0  (comp)
f0 <- function(S, I, th){
  V <- th[1]
  Km <- th[2]
  Kic <- th[3]
  V * S / (Km * (1 + I/Kic) + S)
}

# Model 1  (non comp)
f1 <- function(S, I, th){
  V <- th[1]
  Km <- th[2]
  Kin <- th[3]
  V * S / ((Km + S) * (1 + I/Kin))
}


par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(xstart1, xstart2, type = "n",xlab="x1",ylab="x2",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5); grid()
points(xstart1, xstart2,col=colorgrad[1], pch=16,cex=0.8)

gradientLegend(valRange=c(1,max.iter),color = c("lightblue","cornflowerblue", "darkblue"),
               nCol = 3,inside = FALSE, pos=.825,dec = 0,n.seg=3,cex.main=2, cex.lab=1.5, cex.axis=1.5)
#text(x=31.5, y = 62, lables=expression("Tttitle"))


y<-f0(xstart1, xstart2, TH0START)+rnorm(nrow(xstart), sd = SIG0)
for(i in 1:length(y)){
  if(y[i]<0) y[i]<-0
}

basia<-data.frame(y,xstart)

while(iter< max.iter){
  iter=iter+1
  print(iter)
  #plot(iter,xlim=c(iter-1,iter+1))
  #text(x=iter+1,y=iter,iter,col="red")
  
  
  m0par.est<-nls(y~b1*basia[,2]/(b2*(1+basia[,3]/b3)+basia[,2]),algorithm="port",lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,Inf), start= list(b1=th0start[1],b2=th0start[2],b3=th0start[3]), data=basia)
  THETA.HAT[iter,1:3]=th0start<-summary(m0par.est)$parameters[,1]
  S0SEq[iter,] =s0start<-summary(m0par.est)$parameters[,2]
  sigma0hat[iter]<-sqrt(sum(residuals(m0par.est)^2)/(nrow(basia)-3))
  
  #-------------------------------------------------
  m1par.est<-nls(y~b1*basia[,2]/((b2+basia[,2])*(1+basia[,3]/b3)),algorithm="port",lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,Inf), start= list(b1=th1start[1],b2=th1start[2],b3=th1start[3]), data=basia)
  THETA.HAT[iter,4:6]=th1start<-summary(m1par.est)$parameters[,1]
  S1SEq[iter,] =s1start<-summary(m1par.est)$parameters[,2]
  sigma1hat[iter]<-sqrt(sum(residuals(m1par.est)^2)/(nrow(basia)-3))
  
  #-------------------------------------------------  
  
  psi=function(xest11,xest22){
    xest1=xest11
    xest2=xest22
    eta11=(th0start[1]*xest1)/ ((th0start[2]*(1+(xest2/th0start[3]))) + (xest1*(1+(((1-lambda0)*xest2)/th0start[3]))) )
    eta22=(th1start[1]*xest1)/ ((th1start[2]*(1+(xest2/th1start[3]))) + (xest1*(1+(((1-lambda1)*xest2)/th1start[3]))) )
    eta33=(eta11-eta22)^2
    eta33  
  }
  
  
  psi.vec=mapply(psi,grid[,1],grid[,2])
  max.sen.crit.t[iter]<-max(psi.vec)
  
  xnew<-c(grid[which.max(psi.vec),1],grid[which.max(psi.vec),2])
  xstart<-rbind(xstart,xnew)
  
  eta0 <- f0(xstart[,1],xstart[,2],th0start)
  eta1 <- f1(xstart[,1],xstart[,2],th1start)
  eta3 <- sum((eta0-eta1)^2)
  sum.crit.t[iter] <- eta3
  norm.sum.crit.t[iter] <- eta3/nrow(xstart)
  
  
  ynew<-f0(xnew[1], xnew[2], TH0START)+rnorm(1, sd = SIG0)
  if(ynew<0) ynew<-0
  
  
  y<-c(y,ynew) 
  
  
  cexxx=(length(which(xstart[,1]==xnew[1]&xstart[,2]==xnew[2])))^(1/2)
  points(xnew[1], xnew[2], type = "p",col=colorgrad[iter], pch=16,cex=0.8*cexxx,cex.main=2, cex.lab=1.5, cex.axis=1.5 )
  
  
  basia<-data.frame(y,xstart)
  
  
  
}

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
###############################
cc1<-which(y==0)
sum(cc1) 

cc2<-numeric(6)
for(k in 1:6){
  ccc<-which(THETA.HAT[,k]==0.001)
  cc2[k]<-sum(ccc) 
} 
cc2
#############################

SIG0
min(sigma0hat)
max(sigma0hat)

SIG1
min(sigma1hat)
max(sigma1hat)

par(mfrow=c(1,1))
#des.plot(design.opt)
##########################################################
max(norm.sum.crit.t)
#[1] 0.2820719


quantile(norm.sum.crit.t)

#0%        25%        50%        75%       100% 
#0.07003184 0.27293980 0.27912167 0.28039853 0.28207192

per0.5 <- quantile(norm.sum.crit.t )[3]-quantile(norm.sum.crit.t )[1]
per0.5 / (max(norm.sum.crit.t )-min(norm.sum.crit.t ))

#compute
which(norm.sum.crit.t >= max(norm.sum.crit.t)/2)[1]

which(norm.sum.crit.t  >= 3*norm.sum.crit.t [1])[1]
length(which(norm.sum.crit.t  >= 3*norm.sum.crit.t[1]))
boxplot(norm.sum.crit.t )

quanT <-  c(0.25,0.5,0.75,0.9,0.95)
quanT


earlyitert <- numeric(length(quanT))
for(i in 1:length(quanT)){
  earlyitert[i] <- which(norm.sum.crit.t >= max(norm.sum.crit.t)*quanT[i])[1]
}
earlyitert
#[1]   2   6  16  73 105
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/maxsen.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), max.sen.crit.t , type = "p",xlab="iter",ylab="Maximum sensitivity",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/maxsen.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), max.sen.crit.t , type = "p",xlab="iter",ylab="Maximum sensitivity",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#-------------------------------
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/criterion.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),sum.crit.t , type = "p",xlab="iter",ylab="Criterion",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/criterion.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),sum.crit.t , type = "p",xlab="iter",ylab="Criterion",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)

dev.off()
#--------------------------------#-------------------------------
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/normcrit.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),norm.sum.crit.t  , type = "p",ylim=c(0,0.3),xlab="iter",ylab="Norm.Criterion T",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=max(norm.sum.crit.t ),lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/normcrit.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),norm.sum.crit.t  , type = "p",ylim=c(0,0.3),xlab="iter",ylab="Norm.Criterion T",pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=max(norm.sum.crit.t ),lty=2,col="red")

dev.off()
#--------------------------------#-------------------------------

#------------------------------------------------

#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigma0.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigma0hat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[0])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG0,lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigma0.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigma0hat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[0])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG0,lty=2,col="red")

dev.off()
#--------------------------------#-------------------------------
#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigma1.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigma1hat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[1])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG1,lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigma1.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter), sigma1hat,ylim=c(0,0.6), type = "p",xlab="iter",ylab=expression(hat(sigma[1])),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(h=SIG1,lty=2,col="red")

dev.off()
#--------------------------------#-------------------------------

#--------------------------------#-------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigmapars0.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),S0SEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S0SEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S0SEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#abline(h=SIG0,lty=2,col="red")

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen"),pch=16,bty="n",cex=1.5, pt.cex = 1)


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigmapars0.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),S0SEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S0SEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S0SEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)
#abline(h=SIG0,lty=2,col="red")

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen"),pch=16,bty="n",cex=1.5, pt.cex = 1)

dev.off()
#--------------------------------#-------------------------------
max(S1SEq)
max(S0SEq)

which(S0SEq[,1] <= SIG0)[1]
which(S0SEq[,2] <= SIG0)[1]
which(S0SEq[,3] <= SIG0)[1]

#########################################################################
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigmapars1.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),S1SEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S1SEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S1SEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen"),pch=16,bty="n",cex=1.5, pt.cex = 1)


dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/sigmapars1.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)


plot(seq(1:max.iter),S1SEq[,1] , type = "p",xlab="iter",ylab="Parameter S.E. estimation",ylim=c(0,3),pch=16,col="blue",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S1SEq[,2] , type = "p",xlab="iter",pch=16,col="red",cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(seq(1:max.iter),S1SEq[,3] , type = "p",xlab="iter",pch=16,col="darkgreen",cex.main=2, cex.lab=1.5, cex.axis=1.5)

lgtxt0=c(expression(paste(hat(sigma[V]),sep="")), 
         expression(paste(hat(sigma[M]),sep="")),
         expression(paste(hat(sigma[K]),sep=""))
)


legend(440,3,legend=lgtxt0,col=c("blue","red","darkgreen"),pch=16,bty="n",cex=1.5, pt.cex = 1)


dev.off()

#########################################################################
TH0START
min(THETA.HAT[,1])
max(THETA.HAT[,1])

min(THETA.HAT[,2])
max(THETA.HAT[,2])

min(THETA.HAT[,3])
max(THETA.HAT[,3])

TH1START
min(THETA.HAT[,4])
max(THETA.HAT[,4])

min(THETA.HAT[,5])
max(THETA.HAT[,5])

min(THETA.HAT[,6])
max(THETA.HAT[,6])

par(mfrow=c(1,1))
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th01.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH0START[1],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th01.eps")


par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,1] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH0START[1],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th02.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH0START[2],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th02.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,2] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH0START[2],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th03.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH0START[3],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th03.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,3] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH0START[3],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th11.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,4] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH1START[1],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th11.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,4] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[V])),pch=16,col="blue")
abline(h=TH1START[1],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th12.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,5] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH1START[2],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th12.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,5] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[M])),pch=16,col="blue")
abline(h=TH1START[2],lty=2,col="red")


dev.off()
#--------------------------------
jpeg("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th13.jpeg")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,6] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH1START[3],lty=2,col="red")

dev.off()
#--------------------------------
setEPS()
postscript("G:/Other computers/My Computer/3files/2021/3-March/seq-witheff/NewSequential-Thes/Tf0/Tf0ZThesis/th13.eps")

par(mgp=c(2.5, 1, 0)) 
par(mar = c(5, 5, 4, 2) + 0.1)

plot(seq(1:max.iter),THETA.HAT[,6] ,ylim=c(min(THETA.HAT),max(THETA.HAT)+1), type = "p",xlab="iter",
     cex.main=2, cex.lab=1.5, cex.axis=1.5,ylab=expression(hat(theta[K])),pch=16,col="blue")
abline(h=TH1START[3],lty=2,col="red")


dev.off()
#--------------------------------



