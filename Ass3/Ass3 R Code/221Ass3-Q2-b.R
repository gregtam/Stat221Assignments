library(mvtnorm)

random.orthogonal <- function(p)
{
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}
generate.A <- function(p)
{
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p)
  lambdas = seq(0.01, 1, length.out=p)
  A = Q %*% diag(lambdas) %*% t(Q)
  return(A)
}

sample.data <- function(dim.n, dim.p, model="gaussian")
{
  # Samples the covariates as normal with the specific correlation
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = seq(0.01, 1, length.out=dim.p)
  A = Q %*% diag(lambdas) %*% t(Q)
  
  X = rmvnorm(dim.n, mean=rep(0, dim.p), sigma=A)
  theta = matrix(1, ncol=1, nrow=dim.p)
  epsilon = rnorm(dim.n, mean=0, sd=1)
  # Data generation
  y = X %*% theta  + epsilon
  
  return(list(Y=y, X=X, A=A, theta=theta))
}

check.data <- function(data)
{
  # Do this to check the data object.
  # 
  nx = nrow(data$X)
  ny = length(data$Y)
  p = ncol(data$X)
  stopifnot(nx==ny, p==length(data$theta))
  lambdas = eigen(cov(data$X))$values
  print(lambdas)
  print(mean(data$Y))
  print(var(data$Y))
  print(1 + sum(cov(data$X)))
}

plot.risk <- function(data, est)
{
  # est = p x niters 
  est.bias = apply(est, 2, function(colum) 
    log10(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))
  plot(log10(1:length(est.bias)),est.bias, type="l", lty=3,xlim=range(2,5))
}
get.risk <- function(data, est)
{
  # est = p x niters 
  est.bias = apply(est, 2, function(colum) 
    log10(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))
  est.bias[length(est.bias)]
}

sample.data.sgd <- function(data, lambda0, plot=T)
{
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  #   lambda0 = 0.01
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (theta.old - ai * lpred * xi) + ai * data$Y[i] * xi    
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  if(plot) {
    plot.risk(data, theta.sgd)
  } else {
    return(theta.sgd)
  }
}
sample.data.asgd <- function(data, lambda0, plot=T)
{
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  #   lambda0 = 0.01
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)^(2/3)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (theta.old - ai * lpred * xi) + ai * data$Y[i] * xi
    theta.new = (theta.old*(i) + theta.new)/(i+1)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  if(plot) {
    plot.risk(data, theta.sgd)
  } else {
    return(theta.sgd)
  }
}
sample.data.implicit <- function(data, lambda0, plot=T)
{
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  #   lambda0 = 0.01
  
  for(i in 1:n) {
    xi = data$X[i, ]
    xi.norm = sum(xi^2)
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    # Implicit SGD
    theta.new = (theta.old  - ai * fi * lpred * xi) + (ai * data$Y[i] * xi - ai^2 * fi * data$Y[i] * xi.norm * xi)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  if(plot) {
    plot.risk(data, theta.sgd)
  } else {
    return(theta.sgd)
  }
}
#Batch
sample.data.batch <- function(data, plot=T)
{
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=99)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  #   lambda0 = 0.01
  
  for(i in 100:n)
  {
    #only calculate for every 100
{
  if(i %% 100 == 0)
  {
    xi = data$X[1:i, ]
    yi = data$Y[1:i]
    linear.fit = lm(yi~xi-1)
    # Batch
    theta.new = as.numeric(linear.fit$coef[1:100])  
  }
  else
  {
    #otherwise, set it to the last multiple of 100
    theta.new = theta.sgd[,floor(i/100)*100]
  }
}
theta.new
theta.sgd = cbind(theta.sgd, theta.new)
  }
colnames(theta.sgd)=NULL
if(plot) {
  plot.risk(data, theta.sgd)
} else {
  return(theta.sgd)
}
}

start_time=Sys.time()
A = generate.A(p=100)
data = sample.data(dim.n = 10000,dim.p = 100)
sgd = sample.data.sgd(data, lambda0=0.01,plot=F)
asgd = sample.data.asgd(data, lambda0=0.01,plot=F)
implicit = sample.data.implicit(data, lambda0=0.01,plot=F)
batch = sample.data.batch(data,plot=F)

descentvalues=c()
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))
task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(task.id==1)
{
  descent_values = sample.data.sgd(data,lambda = 0.01,plot=F)
  #   plot.risk(descent_values$data,descent_values$theta)
}
if(task.id==2)
{
  descent_values = sample.data.asgd(data,lambda = 0.01, plot=F)
  #   plot.risk(data,descent_values)
}
if(task.id==3)
{
  descent_values = sample.data.implicit(data,lambda = 0.01, plot=F)
  #   plot.risk(data,descent_values)
}
if(task.id==4)
{
  descent_values = sample.data.batch(data,plot=F)
  #   plot.risk(data,descent_values)
}
end_time = Sys.time()
total_time = as.numeric(end_time-start_time, units="mins")

save(list=c("descent_values","data","start_time","end_time","total_time"),file=sprintf("odyssey/pset3/2b_job%d_task%d.rda", job.id, task.id))









