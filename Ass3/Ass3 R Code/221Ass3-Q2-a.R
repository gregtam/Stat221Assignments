library(mvtnorm)

random.orthogonal = function(p)
{
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}


#SGD
sample.data.sgd = function(dim.t, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = c(1,1,1,rep(0.02,97))
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.t, mean=rep(0, dim.p), sigma=diag(1,dim.p))
  theta = matrix(1, nrow=dim.t, ncol=dim.p)
  # Data generation
  
  gamma_t = 1/(1+0.02 * 1:dim.t)
  
  Sys.time()
  for(i in 2:dim.t)
    theta[i,] = theta[i-1,] - gamma_t[i-1]*2*A %*% (theta[i-1,] - X[i-1,])
  Sys.time()
  
  sgd_risk = sapply(1:dim.t, function(i) theta[i,] %*% A %*% theta[i,])
  return(list(Yplot=log10(sgd_risk)))
}
#ASGD Good
sample.data.asgd.good = function(dim.t, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = c(1,1,1,rep(0.02,97))
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.t, mean=rep(0, dim.p), sigma=diag(1,dim.p))
  theta = matrix(1, nrow=dim.t, ncol=dim.p)
  
  # Data generation
  gamma_t = (1+0.02 * 1:dim.t)^(-2/3)
  
  theta_mean_mat = matrix(1,nrow=dim.t,ncol=dim.p)
  Sys.time()
  for(i in 2:dim.t)
  {
    theta[i,] = theta[i-1,] - gamma_t[i-1]*2*A %*% (theta[i-1,] - X[i-1,])
    #calculates the next average based on previous average, size, and the next theta vector
    theta_mean_mat[i,] = (theta_mean_mat[i-1,]*(i-1)+theta[i,])/i
  }
  Sys.time()
  
  asgd_good_risk = sapply(1:dim.t, function(i) theta_mean_mat[i,] %*% A %*% theta_mean_mat[i,])
  return(list(Yplot=log10(asgd_good_risk)))
}
#ASGD Bad
sample.data.asgd.bad = function(dim.t, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = c(1,1,1,rep(0.02,97))
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.t, mean=rep(0, dim.p), sigma=diag(1,dim.p))
  theta = matrix(1, nrow=dim.t, ncol=dim.p)
  # Data generation
  
  gamma_t = (1+1:dim.t)^(-1/2)
  
  theta_mean_mat = matrix(1,nrow=dim.t,ncol=dim.p)
  Sys.time()
  for(i in 2:dim.t)
  {
    theta[i,] = theta[i-1,] - gamma_t[i-1]*2*A %*% (theta[i-1,] - X[i-1,])
    #calculates the next average based on previous average, size, and the next theta vector
    theta_mean_mat[i,] = (theta_mean_mat[i-1,]*(i-1)+theta[i,])/i
  }
  Sys.time()
  
  asgd_bad_risk = sapply(1:dim.t, function(i) theta_mean_mat[i,] %*% A %*% theta_mean_mat[i,])
  return(list(Yplot=log10(asgd_bad_risk)))
}
#Implicit
sample.data.implicit = function(dim.t, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = c(1,1,1,rep(0.02,97))
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.t, mean=rep(0, dim.p), sigma=diag(1,dim.p))
  theta = matrix(1, nrow=dim.t, ncol=dim.p)
  # Data generation
  
  gamma_t = (1+0.02 *1:dim.t)^(-1/2)
  
  Sys.time()
  for(i in 2:dim.t)
    theta[i,] = solve(diag(dim.p) + gamma_t[i-1]*2*A) %*% (theta[i-1,] + gamma_t[i-1]*2*A %*% X[i-1,])
  Sys.time()
  
  implicit_risk = sapply(1:dim.t, function(i) theta[i,] %*% A %*% theta[i,])
  return(list(Yplot=log10(implicit_risk)))
}
#Batch
sample.data.batch = function(dim.t, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = c(1,1,1,rep(0.02,97))
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.t, mean=rep(0, dim.p), sigma=diag(1,dim.p))
  theta = matrix(1, nrow=dim.t, ncol=dim.p)
  # Data generation
  
  gamma_t = (1+0.02 *1:dim.t)^(-1)
  
  Sys.time()
  theta[1,] = X[1,]
  for(i in 2:dim.t)
    theta[i,] = (theta[i-1,]*(i-1) + X[i,])/i
  Sys.time()
  
  batch_risk = sapply(1:dim.t, function(i) theta[i,] %*% A %*% theta[i,])
  return(list(Yplot=log10(batch_risk)))
}


start_time = Sys.time()
data.sgd = sample.data.sgd(dim.t = 1000000, dim.p = 100)
data.asgd.good = sample.data.asgd.good(dim.t = 1000000, dim.p = 100)
data.asgd.bad = sample.data.asgd.bad(dim.t = 1000000, dim.p = 100)
data.implicit = sample.data.implicit(dim.t = 1000000, dim.p = 100)
data.batch = sample.data.batch(dim.t = 1000000, dim.p = 100)
end_time = Sys.time()
total_time =  as.numeric(end_time-start_time,units="mins")


save(list=c("data.sgd","data.asgd.good","data.asgd.bad","data.implicit","data.batch","start_time","end_time","total_time"), file="odyssey/pset3/2a.rda")











