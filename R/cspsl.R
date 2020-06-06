

###########################
#  Correlated spsl model  #
###########################
if("rjags" %in% rownames(installed.packages()) == FALSE) {install.packages("rjags")}
library(rjags)
cspsl <- function(data,
                  z,
                  n.iter = 5000,
                  n.chains = 1) {
  fisher <- function(x)
    return(1 / 2 * log((1 + x) / (1 - x)))
  inv.fisher <- function(x)
    return((exp(2 * x) - 1) / (exp(2 * x) + 1))
  nb = ncol(data) # number of gene
  n = nrow(data)
  q = choose(nb, 2) # size of combination
  getY <- function(x, nb) {
    y = c()
    col = nb
    for (j in 1:(col - 1)) {
      for (k in (j + 1):col) {
        y = cbind(y, x[, j] * x[, k])
      }
    }
    return(y)
  }
  y <- getY(scale(x), nb = nb)
  omega <- solve(cor(y))
  # correlated spike and slab
  model1_stringM1 <- "
  model{
    for(i in 1:n){
      for(j in 1:q){
        y[i,j] ~ dnorm(1-2/(exp(2*mu[i,j])+1), prec[j])
        mu[i,j]<-tau0[j]+z[i]*tau1[j]
      }
    }
    for (j in 1:q) {
      tau1[j]~dnorm(0,s[j])
      s[j]<-1/(vars[j])
      vars[j]<-tau[j]*r[j]
      tau[j]<-1/invtau[j]
      invtau[j]~dgamma(5,50)
      r[j]<-w[j]+0.005*w[j]
      w[j]<-pnorm(p[j],0,1)
      prec[j]~dgamma(0.01,0.01)
      tau0[j]~dnorm(0,0.001)
    }
    p~dmnorm(rep(0,10),omega)
  }
  "
  jags_data = list(
    n = n,
    y = y,
    z = z,
    q = q,
    omega = omega
  )
  n.iter = n.iter
  time1 <- Sys.time()
  jags_model = try(jags.model(
    textConnection(model1_stringM1),
    data = jags_data,
    n.adapt = 0,
    n.chains = n.chains
  ))
  class(jags_model)
  update(jags_model, n.iter = n.iter)
  coda_sample = coda.samples(
    jags_model,
    c("tau0", "tau1", "w", "vars", "r"),
    n.iter = n.iter,
    n.burnin = floor(n.iter / 2)
  )
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(coda_sample, as.numeric(time, units = "secs")))
}
