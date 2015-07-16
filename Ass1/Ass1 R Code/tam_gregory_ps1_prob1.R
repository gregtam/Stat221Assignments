temp = read.csv("data1985_area2.csv",sep="\t")
ldata = read.table("dataLogisticNorm3D.txt", header = TRUE)

#1.3
dlogisticnorm = function(u,mu,alpha,beta)
{
  d = length(u)
  Sigma = matrix(-beta, nrow=d, ncol=d)
  for(i in 1:d)
    Sigma[i,i] = alpha
  x = log(u)
  ans = (2*pi)^(-d/2)*abs(det(Sigma))^(-1/2)*exp(-1/2*t(log(u/sum(u))-mu) %*% solve(Sigma) %*% (log(u/sum(u))-mu))
  as.numeric(ans)
}

dlogisticnorm(c(0.6,0.1,.03),c(2,3,4),20,17)


initlog = c(1,1,1,1,1,1)
optim(list(mu=c(1,1), alpha=2,beta=1), dlogisticnorm, control=list(fnscale=-1), method="L-BFGS-B")

optim(init, ll.thetaH, control=list(fnscale=-1), method = "L-BFGS-B", theta = optTheta,X=X,G=optG,index=index)$par




logisticnorm.mle = function(U)
{
  wrapper = function(param)
  {
    tempMu = param[3:length(param)]
    d=length(tempMu)
    
    tempAlpha = abs(param[1])
    tempBeta = abs(param[2])
    
    loglikeli=0
    for(row in 1:nrow(U))
      loglikeli = loglikeli + log(dlogisticnorm(as.numeric(U[row,]),tempMu, tempAlpha, tempBeta))
    loglikeli
  }
  tempparam = c(1,1,rep(1,ncol(ldata)))
  tempoptim = optim(tempparam, wrapper, control=list(fnscale=-1))$par
  
  mu.hat = tempoptim[3:length(tempoptim)]
  alpha.hat = tempoptim[1]
  beta.hat = tempoptim[2]
  list("mu.hat" = mu.hat, "alpha.hat" = alpha.hat, "beta.hat" = beta.hat)
}

logisticnorm.mle(ldata[1,])
logisticnorm.mle(ldata)




