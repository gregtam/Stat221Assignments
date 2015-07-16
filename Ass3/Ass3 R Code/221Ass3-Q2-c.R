library(mvtnorm)

bias = function(theta)
{
  #returns the bias assuming theta.star is a vector of 1s
  log(sqrt(sum((theta-1)^2)))
#   sum(abs(theta-1))
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

lr <- function(alpha, n)
{
  ## learning rate
  alpha / (alpha + n)
}

tr = function(A)
{
  sum=0
  {
    if(nrow(A)!=ncol(A))
      stop("Not a square matrix")
    else
    {
      for(i in 1:nrow(A))
        sum = sum + A[i,i]
    }
  }
  sum
}

sample.data <- function(dim.n, A, model="gaussian")
{
  # Samples the dataset. Returns a list with (Y, X, A ,true theta)
  dim.p = nrow(A)
  # This call will make the appropriate checks on A.
  X = rmvnorm(dim.n, mean=rep(0, dim.p), sigma=A)
  theta = matrix(1, ncol=1, nrow=dim.p)
  epsilon = rnorm(dim.n, mean=0, sd=1)
  # Data generation
  y = X %*% theta  + epsilon
  return(list(Y=y, X=X, A=A, theta=theta))
}
#SGD
sample.data.sgd <- function(data, alpha, plot=T)
{
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = alpha / (alpha * tr(A) + n)
    
    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$Y[i]
    # Standard SGD
    theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  return(list(data = data, theta = theta.sgd))
}
#ASGD
sample.data.asgd <- function(data, alpha, plot=T)
{
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = alpha / (alpha * tr(A) + n)
    
    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$Y[i]
    # ASGD
    theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi
    theta.new = (theta.old*(i-1) + theta.new)/i
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  colnames(theta.sgd)=NULL
  return(list(data = data, theta = theta.sgd))
}
#Implicit
sample.data.implicit <- function(data, alpha, plot=T)
{
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)
  
  for(i in 1:n)
  {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = lr(alpha, i)
    
    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$Y[i]
    # Implicit SGD
    theta.new = (theta.old  - ai * fi * lpred * xi) +  
      (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
    
    theta.sgd = cbind(theta.sgd, theta.new)
  }  
  colnames(theta.sgd)=NULL
  return(list(data = data, theta = theta.sgd))
}

start_time = Sys.time()
task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))
#Vectors to convert task.id to method and alpha value
method_vec = floor(0:29 / 10)+1
#generates different alpha values to sample over
alpha_index = rep(seq(50,500,by=50),3)

method = method_vec[task.id]
alpha = alpha_index[task.id]
args <- as.numeric(commandArgs(trailingOnly = TRUE))
m = args[1]
if(length(args) != 1)
{
  stop("Not correct no. of args")
}

theta_mat = matrix(,nrow=0,ncol=10000)
for(i in 1:m)
{
  A = generate.A(p=100)
  d = sample.data(dim.n=1e4-1, A)
  descent_values = c()
  if(method==1)
  {
    descent_values = sample.data.sgd(d,alpha=alpha)
  }
  if(method==2)
  {
    descent_values = sample.data.asgd(d,alpha=alpha)
  }
  if(method==3)
  {
    descent_values = sample.data.implicit(d,alpha=alpha)
  }
  theta_mat = rbind(theta_mat,descent_values$theta)
}

theta_avg = c()
get_bias = function(num, theta_mat)
{
  #take average of bias(theta_mat[1:100,num]), ..., bias(theta_mat[1:100+m*100,num])
  #bias for theta_num across all m replications
  mean(sapply(1:m,function(i) bias(theta_mat[(i-1)*100 + 1:100,num])))
}
biases = sapply(1:1e4, function(num) get_bias(num,theta_mat))

get_variance = function(num,theta_mat)
{
  #gets the variance-covariance matrix for a given theta_num across m replications
  temp_mat = matrix(,nrow=0,ncol=100)
  for(i in 1:m)
    temp_mat = rbind(temp_mat,theta_mat[(i-1)*100 + 1:100,num])
  tr(var(temp_mat))
}
variances = sapply(1:1e4, function(num) get_variance(num, theta_mat))

end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")
save(list = c("biases","variances","start_time","end_time","total_time"), file = sprintf("odyssey/pset3/2c/2c_job%d_task%d.rda",  job.id, task.id))















