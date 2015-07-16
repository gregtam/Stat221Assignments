impala = read.table("impala.txt", header = TRUE)
impala = impala[,1]

waterbuck = read.table("waterbuck.txt", header = TRUE)
waterbuck = waterbuck[,1]

log.lik <- function(Y, N, theta)
{
  # Log-likelihood of the data
  sum(dbinom(Y, N, theta, log=T))
}

log.prior <- function(N, theta)
{
  #1/N times the indicator that theta is in [0,1]
  log(1/N * ifelse(theta>=0 & theta<=1,1,0))
}

time = 1:1000

#impala
a = 0.0005
b = 3
gamma = 0.7

#waterbuck
a = 0.002
b = 5
gamma = 0.7

epsilon = a*(b+time)^(-gamma)
epsilon[1]
epsilon[length(epsilon)]

y = impala
y = waterbuck
S = sum(y)
theta.samp = 0.9
N.samp = max(y)

for(t in 1:100)
{
  theta.old = tail(theta.samp,1)
  N.old = tail(N.samp,1)

  #ensures N is large enough
  repeat
  {
    #geometric over entire thing
    N.new = rgeom(1,1/(1+N.old))
    #only keep sample if N is in the boundaries
    if(N.new >= max(y) & N.new <= 300)
      break
  }

  logprior.grad = -log(N.old)
  ll.grad = function(index)
  {
    y[index]/theta.old - (N.old - y[index])/(1-theta.old)
  }

  ll.grad.vec = sapply(1:5,ll.grad)
  eta = rnorm(1,mean=0,sd=sqrt(epsilon[t]))

  #proposed update
  theta.new = theta.old + epsilon[t]/2*(logprior.grad + sum(ll.grad.vec)) + eta
  #if theta.new is outside boundaries, project back to [0,1]
  #theta.new cannot be 0 or 1 otherwise ll.grad becomes infinite
  #so we add a bit of noise
  if(theta.new>1)
    theta.new=1-1e-6
  if(theta.new<0)
    theta.new=1e-6

  #Add theta to samples
  theta.samp = c(theta.samp,theta.new)
  theta.old = theta.new
  #Add N to samples
  N.samp = c(N.samp,N.new)
  N.old = N.new
}

N.samp
theta.samp



pdf("waterbuck-scatterplot.pdf",width=6,height=6)
plot(N.samp,theta.samp,type="b",ylim=range(0,1))
dev.off()


pdf("waterbuck-diagnostics.pdf",width=6,height=6)
par(mfrow=c(2,2))
acf(N.samp, main="ACF for N")
acf(theta.samp, main="ACF for theta")
plot(N.samp,type='l',main="Traceplot for N")
plot(theta.samp,type='l',main="Traceplot for theta")
dev.off()


pdf("impala-scatterplot.pdf",width=6,height=6)
plot(N.samp,theta.samp,type="b",ylim=range(0,1))
dev.off()


pdf("impala-diagnostics.pdf",width=6,height=6)
par(mfrow=c(2,2))
acf(N.samp, main="ACF for N")
acf(theta.samp, main="ACF for theta")
plot(N.samp,type='l',main="Traceplot for N")
plot(theta.samp,type='l',main="Traceplot for theta")
dev.off()


