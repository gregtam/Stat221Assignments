library(MASS) #ginv
library(mvtnorm)
dat = read.csv("1router_allcount.dat")
attach(dat)

c=2
w=11
h=(w-1)/2
T=length(nme)/25
I=16
J=8

task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))

x.indices = 1:16
y.indices = 17:24
nme[which(nme %in% nme[x.indices])]
x = matrix(value[which(nme %in% nme[x.indices])],nrow=16)
nme[which(nme %in% nme[y.indices])]
y = matrix(value[which(nme %in% nme[y.indices])],nrow=8)

#set A matrix
A = matrix(,nrow=J-1,ncol=I)
A[1,]=c(rep(1,4),rep(0,12))
A[2,]=c(rep(0,4),rep(1,4),rep(0,8))
A[3,]=c(rep(0,8),rep(1,4),rep(0,4))
A[4,]=c(rep(0,12),rep(1,4))
A[5,]=rep(c(1,0,0,0),4)
A[6,]=rep(c(0,1,0,0),4)
A[7,]=rep(c(0,0,1,0),4)

# eta.hat0 = 1
# Sigma0=diag(rep(1,16))
# V=matrix()
# 
# #Step 1.
# eta.hat=c()
# t=h+1
# eta.hat[t-1]=eta.hat0
# Sigma.hat=Sigma0
# Sigma.hatV = Sigma.hat
# 
# #Step 2.
# eta.hat[t] = rnorm(1,mu=)

c=2
lambda = rep(10,16)
phi = 20
theta = matrix(c(lambda,phi),ncol=1)
eta = log(theta)

t = h+1
eta.old=c()
# Sigma0 = diag(c(phi*lambda^c,))
Sigma0 = diag(rep(1,17))
V = diag(rep(1,17))
eta.hat.old = eta.hat0
Sigma.hat.old = Sigma0

repeat
{
  #Step 2
  Sigma.hat.tt = Sigma.hat.old + V
  
  #Step 3
  dmvnorm(c(eta),mean=eta.hat.old,sigma=Sigma.hat.tt,log=TRUE)
  dmvnorm(y[-8,1], mean=A%*%lambda, sigma=A%*%Sigma%*%t(A),log=TRUE)
  g = function(eta.t)
  {
    lambda.temp = exp(eta.t[-length(eta.t)])
    phi.temp = exp(eta.t[length(eta.t)])
    Sigma.temp = diag(phi.temp*lambda.temp^c)
    ans = dmvnorm(c(eta),mean=eta.t,sigma=Sigma.hat.tt,log=TRUE) +
        dmvnorm(y[-8,1], mean=A%*%lambda.temp, sigma=A%*%Sigma.temp%*%t(A),log=TRUE)
    ans
  }
  eta.hat.new = optim(eta.hat.old,g,control=list(fnscale=-1))$par
  Sigma.hat.new = -ginv(Sigma.hat.tt)
}














