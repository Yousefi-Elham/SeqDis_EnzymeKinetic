xstart = matrix(c(30.00, 0.00, 0.086,
                  3.884, 0.00, 0.208,
                  30.000, 16.952,0.206,
                  3.884, 16.952, 0.500), nrow=4,ncol=3,byrow=TRUE)

xstart

thEstart <- c(7.4253, 4.6808, 3.0581,0.9636)
sEstart <- c(0.1298, 0.2724, 0.2815,0.0191)
TH0START<-c(7.298, 4.386, 2.582)
TH1START<- c(8.696, 8.066, 12.057)
TTTT <- c(TH1START,0)

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

deriv.result <- F.deriv(TTTT,xstart[,c(1,2)])
FE.full.mat <- deriv.result$FE.full
FE3.full.mat <- deriv.result$FE3.full


wstart <- xstart[,3]
W <- diag(wstart)
M.mat <- t(FE.full.mat) %*% W %*% FE.full.mat
ME3.mat <- t(FE3.full.mat) %*% W %*% FE3.full.mat
Cri.val <- det(M.mat)/det(ME3.mat)
Cri.val
#[1]  0.02527844