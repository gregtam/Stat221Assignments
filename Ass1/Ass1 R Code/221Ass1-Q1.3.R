library(mvtnorm)
library(expm)

dlogisticnorm = function(u,mu,alpha,beta)
{
  #u_{d+1}
  u_d.1 = 1-sum(u)
  X = matrix(-beta, nrow=length(u), ncol=length(u))
  for(i in 1:length(u))
    X[i,i]=alpha
  E = eigen(X)
  V = E$values 
  Q = E$vectors 
  Y = Q%*%diag(1/sqrt(V))%*%t(Q) 
}


logisticnorm.mle = function(U)
{
  
}






X = matrix(c(1,0.4,0.4,1),nrow=2)
E = eigen(X) 
V = E$values 
Q = E$vectors 
Y = Q%*%diag(1/sqrt(V))%*%t(Q) 






