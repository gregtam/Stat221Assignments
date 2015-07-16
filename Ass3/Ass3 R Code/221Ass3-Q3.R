rm(list=ls())
library(mvtnorm)
library(glmnet)

# genjerry, genx2 are functions taken from the above paper.
# These functions generate the simulation data.
# NOTE: use function sample.data(..) instead.
genjerry = function(x, snr)
{
  # generate data according to Friedman's setup
  n=nrow(x)
  p=ncol(x)
  b=((-1)^(1:p))*exp(-2*((1:p)-1)/20)
  # b=sample(c(-0.8, -0.45, 0.45, 0.9, 1.4), size=p, replace=T)
  # ((-1)^(1:p))*(1:p)^{-0.65}#exp(-2*((1:p)-1)/20)
  f=x%*%b
  e=rnorm(n)
  k=sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, beta=b))
}

genx2 = function(n,p,rho)
{
  #    generate x's multivariate normal with equal corr rho
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}
  
  return(x)
}

sample.data <- function(dim.n, dim.p, rho=0.0, snr=1)
{
  # Samples the dataset according to Friedman et. al.
  #
  # 1. sample covariates
  X = genx2(dim.n, dim.p, rho)
  # 2. ground truth params.
  theta = ((-1)^(1:dim.p))*exp(-2*((1:dim.p)-1)/20)
  
  f= X %*% theta
  e = rnorm(dim.n)
  k= sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, X=X, theta=theta))
}

dist <- function(x, y)
{
  if(length(x) != length(y))
    stop("MSE should compare vectors of same length")
  sqrt(mean((x-y)^2))
}

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

sample.data.sgd <- function(data, lambda0, plot=T)
{
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1/sum(diag(data$X %*% t(data$X)))  
  #   lambda0 = 0.01
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)^(2/3)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (theta.old - ai * lpred * xi) + ai * data$y[i] * xi    
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  if(plot) {
    plot.risk(data, theta.sgd)
  } else {
    return(theta.sgd)
  }
}

# Main function to run this experiment.
run.glmnet <- function(dim.n, dim.p,
                       rho.values=c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95),
                       nreps=3, 
                       verbose=F)
{
  ## Runs glmnet() for various param values.
  ##
  niters = 0
  cols = c("rho", "rep", "time", "mse")
  timings = matrix(nrow=0, ncol=length(cols))
  colnames(timings) <- cols
  rownames(timings) = NULL
  total.iters = nreps * length(rho.values)
  
  pb = txtProgressBar(style=3)
  
  seeds=sample(1:1e9, size=total.iters)
  for(i in 1:nreps)
  {
    for(rho in rho.values)
    {
      niters = niters + 1
      set.seed(seeds[niters])
      # 1. (For every repetition) Sample the dataset
      dataset = sample.data(dim.n=dim.n, dim.p=dim.p, rho=rho, snr=3.0)
      
      lambda0 = min(eigen(dataset$X %*% t(dataset$X))$values)
      
      true.theta = dataset$theta
      x = dataset$X
      y = dataset$y
      stopifnot(nrow(x) == dim.n, ncol(x) == dim.p)
      # 1b. Define metrics:
      #   dt = time for the method to finish
      #   mse = Distance (e.g. RMSE) of the estimates to the ground truth.
      #         (q1, q2, q3) representing the quartiles (since glmnet returns grid of estimates)
      #         Implicit has (x, x, x) i.e., the same value in all places.
      new.dt = 0
      new.mse = NA
      # 2. Run the method.
    
      new.dt = system.time({fit = sample.data.sgd(dataset,lambda0,plot=F)})[1]
      new.mse = median(apply(fit, 2, function(est) dist(est, true.theta)))
      
            
      # 3. Tabulate timings
      timings = rbind(timings, c(rho, i, 
                                 new.dt, 
                                 new.mse))
      setTxtProgressBar(pb, niters/total.iters)
    }
  }
  return(timings)
}


#3c
times = matrix(,nrow=0,ncol=6)
params = cbind(N=c(1000,5000,100,100,100,100), p=c(100,100,1000,5000,20000,50000))

start_time = Sys.time()
task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))
N = params[task.id,"N"]
p = params[task.id,"p"]
ans = run.glmnet(N, p)
times = rbind(times,tapply(ans[,"time"],c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)), mean))
end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")

save(list=c("times","start_time","end_time","total_time"),file=sprintf("odyssey/pset3/3c_job%d_task%d.rda", job.id, task.id))


#3d
times = matrix(,nrow=0,ncol=6)

start_time = Sys.time()
task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))
print(paste("task.id ", task.id,"   ---  job.id  ",job.id))
N = 50000
p = 10e3
ans = run.glmnet(N, p)
times = rbind(times,tapply(ans[,"time"],c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)), mean))
end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")

save(list=c("times","start_time","end_time","total_time"),file=sprintf("odyssey/pset3/3d_job%d_task%d.rda", job.id, task.id))
