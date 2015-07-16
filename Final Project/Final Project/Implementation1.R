sigma = c(sqrt(10),1)
sigma.x = sqrt(2)

theta = matrix(c(0,1),nrow=2,ncol=1)
# theta = matrix(c(1,-1),nrow=2,ncol=1)
eta = matrix(,nrow=2,ncol=0)

x = 1/2*rnorm(100, mean=theta[1,1], sd=sigma.x) +
  1/2*rnorm(100, mean=theta[1,1] + theta[2,1], sd = sigma.x)

N=length(x)
n=1

time = 1:10000
a = 0.015
b = 1
gamma = 0.55


epsilon = a*(b+time)^(-gamma)
epsilon[1]
epsilon[length(epsilon)]

t=1
for(t in 1:10000)
{
  #set most recent theta as previous theta
  theta.old = theta[,ncol(theta)]
  #gradient of the log of the prior
  logprior.grad = -theta.old/sigma
  #gradient of the loglikelihood
  ll.grad = function(xtemp)
  {
    const = 1/(exp(-(xtemp-theta.old[1])^2/(2*sigma.x^2)) +
         exp(-(xtemp-theta.old[1]-theta.old[2])^2/(2*sigma.x^2)))
    vec = c()
    vec[1] = exp(-(xtemp - theta.old[1])^2/(2*sigma.x^2)) * (xtemp - theta.old[1])/sigma.x^2 +
      exp(-(xtemp - theta.old[1] - theta.old[2])^2/(2*sigma.x^2)) * (xtemp - theta.old[1] - theta.old[2])/sigma.x^2
    vec[2] = exp(-(xtemp - theta.old[1] - theta.old[2])^2/(2*sigma.x^2)) *
      (xtemp - theta.old[1] - theta.old[2])/sigma.x^2
    const*vec
  }

  eta = cbind(eta,rnorm(2,mean=0,sd=sqrt(epsilon[t])))
  ll.grad.vec = sapply(x[1:n],ll.grad)

  #proposed update
  theta.new = theta.old + epsilon[t]/2*(logprior.grad + N/n*apply(ll.grad.vec,1,sum)) + eta[,ncol(eta)]
  theta = cbind(theta,theta.new)
  theta.old = theta.new
}

plot(theta[1,],theta[2,],xlim=range(-1.5,2.5),ylim=range(-3,3))

pdf("fig1right.pdf",width=6,height=6)
plot(theta[1,],theta[2,],xlim=range(-1.5,2.5),ylim=range(-3,3),col=rgb(0,0,0,alpha=0.1),main="Estimated Posterior")
dev.off()

# plot(c(),xlim=range(-1.5,2.5),ylim=range(-3,3))
# text(theta[1,],theta[2,],labels=c(1:10000))



#Contours of Posterior
library(lattice)
prior = function(theta1,theta2){dnorm(theta1,mean=0,sd=sigma[1])*dnorm(theta2,mean=0,sd=sigma[2])}
likelihood = function(x,theta1,theta2)
{
  1/2*dnorm(x,mean=theta1,sd=sigma.x) + 1/2*dnorm(x,mean=theta1+theta2,sd=sigma.x)
}
C = 0
for(theta1 in seq(-5,5,by=0.25))
{
  for(theta2 in seq(-5,5, by=0.25))
  {
    total_likelihood = 1
    for(i in 1:length(x))
      total_likelihood = total_likelihood * likelihood(x[i],theta1,theta2)
    C = C + prior(theta1, theta2) * total_likelihood
  }
}
C
posterior = function(theta1,theta2,x)
{
  post = prior(theta1,theta2)
  for(i in 1:length(x))
    post = post * likelihood(x[i],theta1,theta2)
  post/C
}

#check if constant works
final_likelihood = 0
posts = c()
for(theta1 in seq(-5,5,by=0.25))
{
  for(theta2 in seq(-5,5, by=0.25))
  {
    final_likelihood = final_likelihood + posterior(theta1, theta2, x)
    posts = c(posts,posterior(theta1, theta2, x))
  }
}
final_likelihood



theta1 = seq(-1.5,2.5,0.25)
theta2 = seq(-3,3,0.25)
ab.grid = expand.grid(theta1=theta1,theta2=theta2)
z = matrix(nrow=length(theta1),ncol=length(theta2))
for(i in 1:length(theta1))
  for(j in 1:length(theta2))
    z[i,j] = posterior(theta1[i],theta2[j],x)
ab.grid$z = c(z)

#Figure 1
pdf("fig1left.pdf",width=6,height=6)
print(contourplot(z ~ theta1*theta2, data=ab.grid, cuts=8, labels=FALSE,main="Posterior Contours"))
dev.off()

#Figure 2
pdf("fig2left.pdf",width=6,height=6)
plot(epsilon,type='l',log="xy",ylim=range(10e-6,10e0),col="red",xlab="iteration",ylab="noise variance")
lines(diff(theta[1,]),col="blue")
lines(diff(theta[2,]),col="green")
legend(1,0.0001,c("theta1 noise","theta2 noise","injected noise"),col=c("blue","green","red"),lwd=c("2","2","2"))
dev.off()






