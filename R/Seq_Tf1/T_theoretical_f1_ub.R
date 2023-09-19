#Tf1
xstart = matrix(c(30.00, 0.00, 0.063,
                  3.214, 0.00, 0.063,
                  30.00, 21.413,0.31,
                  5.625, 11.152,0.564), nrow=4,ncol=3,byrow=TRUE)

xstart

th0start <- c(7.298, 4.386, 2.582)
s0start <- c(0.114, 0.233, 0.145)
th1start <- c(8.696, 8.066, 12.057)
s1start<- c(0.222, 0.488, 0.671)


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

eta0 <- f0(xstart[,1],xstart[,2],th0start)
eta1 <- f1(xstart[,1],xstart[,2],th1start)
Crit.val <- sum(xstart[,3]*((eta0-eta1)^2))
Crit.val 
#[1] 0.2677787