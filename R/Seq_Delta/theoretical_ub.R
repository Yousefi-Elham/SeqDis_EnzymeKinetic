# calculation of the theoretical 



par(mfrow=c(1,1))
rm(list=ls())
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
                  30.00, 60.00), nrow=9,ncol=2,byrow=TRUE )


xstart1=xstart[,1]
xstart2=xstart[,2]

th0start <- c(7.298, 4.386, 2.582)
s0start <- c(0.114, 0.233, 0.145)
th1start <- c(8.696, 8.066, 12.057)
s1start<- c(0.222, 0.488, 0.671)

TH0START<-c(7.298, 4.386, 2.582)
TH1START<- c(8.696, 8.066, 12.057)
SIG0<-0.1553
SIG1<-0.2272

x1g=seq(0,30,length=31)
x2g=seq(0,60,length=61)
grid=expand.grid(x1=x1g,x2=x2g)
grid=cbind(grid[,1],grid[,2])



iter=1
max.iter=500

design.optFdelta <- matrix( c(0,  0,       1,    0.002,
                              0, 30,       1,    0.002,
                              0, 60,       1,    0.002,
                              15,  0,       1,    0.002,
                              15, 30,       1,    0.002,
                              15, 60,       1,    0.002,
                              30,  0,      29,    0.058,
                              30, 30,       1,    0.002,
                              30, 60,       1,    0.002,
                              8, 15,       1,    0.002,
                              5,  0,       1,    0.002,
                              6, 10,       3,    0.006,
                              30, 15,       1 ,   0.002,
                              6, 12,     155 ,   0.310,
                              4,  0,       1,    0.002,
                              30, 17,       1,    0.002,
                              6, 11,       6,    0.012,
                              3,  0,     123,    0.246,
                              5, 10,       4,    0.008,
                              30, 18,       1,    0.002,
                              5, 11 ,     42,    0.084,
                              30, 23,      98,    0.196,
                              30, 24,       1,    0.002,
                              30, 21,       3,    0.006,
                              30, 22,      18,    0.036,
                              30, 20,       4,    0.008,
                              5, 12,       9 ,   0.018) ,nrow=27,ncol=4,byrow=T)

th0l <-th0start-s0start
th0u <-th0start+s0start
th1l <-th1start-s1start
th1u <-th1start+s1start



delta.func=function(xdes1,xdes2){
  
  
  N <- length(xdes1)
  
  # Calculate F0, F1, a0, a1 for all design points
  F0.full <- matrix(0, nrow = N, ncol = 3)
  a0.full <- rep(0, N)
  F1.full <- matrix(0, nrow = N, ncol = 3)
  a1.full <- rep(0, N)
  
  
  
  
  for(ii in 1:N){
    
    x1 <- xdes1[ii]
    x2 <- xdes2[ii]
    
    a <- th0start[1] * x1 
    b <- 1 + x2 / th0start[3]
    d <- b * th0start[2] + x1
    F0.full[ii,1] <- x1 / d
    F0.full[ii,2] <- -a / d^2 * b
    F0.full[ii,3] <- a / d^2 * (th0start[2] * x2 / th0start[3]^2)
    a0.full[ii] <- a / d - sum(F0.full[ii,] * th0start)
    
    a <- th1start[1] * x1 
    b <- 1 + x2 / th1start[3]
    d <- b * (th1start[2] + x1)
    F1.full[ii,1] <- x1 / d
    F1.full[ii,2] <- -a / d^2 * b
    F1.full[ii,3] <- a / d^2 * ((th1start[2] + x1) * x2 / th1start[3]^2)
    a1.full[ii] <- a / d - sum(F1.full[ii,] * th1start)
    
  }
  
  
  crit_delta_sq<-bvls(cbind(F0.full, -F1.full), a1.full - a0.full, c(th0l,th1l), c(th0u,th1u))$deviance
  
  crit_delta_sq
  
  
}

delta.func(design.optFdelta[,1],design.optFdelta [,2])/nrow(design.optFdelta)
