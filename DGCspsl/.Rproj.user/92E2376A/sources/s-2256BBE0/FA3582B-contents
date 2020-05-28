############################
# Generate data of x and z #
############################
get_data <- function(n, nb, tau0, tau1){
if("mvtnorm" %in% rownames(installed.packages()) == FALSE) {install.packages("mvtnorm")}
if("matrixcalc" %in% rownames(installed.packages()) == FALSE) {install.packages("matrixcalc")}
if("Matrix" %in% rownames(installed.packages()) == FALSE) {install.packages("Matrix")}
library(mvtnorm)
library(matrixcalc)
library(Matrix)
fisher<-function(x) return(1/2*log((1+x)/(1-x)))
inv.fisher<-function(x) return((exp(2*x)-1)/(exp(2*x)+1))
nb <- nb # number of gene
q <- choose(nb,2) # size of combination
n <- n
# setting dynamic correlation
rho<-inv.fisher(t(t(outer(z,tau1))+tau0))
# Generating z
z<-runif(n,0,1)
# Generating x
time1<-Sys.time()
x<-matrix(rep(0,nb*n),nrow=n)
for (i in 1:n){
  cor.mat <- diag(rep(0.5,nb))
  cor.mat[lower.tri(cor.mat)] <- rho[i,]
  cor.mat<-cor.mat+t(cor.mat)
  while(is.positive.definite(cor.mat)==F){
    z[i]<-runif(1)
    rho[i,]<-inv.fisher(t(t(outer(z[i],tau1))+tau0))
    cor.mat <- diag(rep(0.5,nb))
    cor.mat[lower.tri(cor.mat)] <- rho[i,]
    cor.mat<-cor.mat+t(cor.mat)
  }
  x[i,]<-rnorm(nb)%*%chol(cor.mat) # Cholesky decomposition
}
time2<-Sys.time()
cat("generating x needs", as.numeric(time2-time1),"seconds\n")
return(list(x, z))
}


