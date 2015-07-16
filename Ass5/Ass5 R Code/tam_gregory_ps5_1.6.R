library(numDeriv) #hessian function
library(mvtnorm)
library(psych) #trace function
library(MASS) #ginnv

dlogprior = function(eta,eta.hat,Sigma.hat)
{
  dmvnorm(eta, mean=eta.hat, sigma=Sigma.hat, log=TRUE)
}
dloglikelihood = function(Y,A,lambda,Sigma)
{
  dmvnorm(Y, mean=A %*% lambda, sigma= A %*% Sigma %*% t(A), log=TRUE)
}

library(MASS) #ginv
dat = read.csv("1router_allcount.dat")
attach(dat)

start_time=Sys.time()

task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))

x.indices = 1:16
y.indices = 17:24
# nme[which(nme %in% nme[x.indices])]
x = matrix(value[which(nme %in% nme[x.indices])],nrow=16)
# nme[which(nme %in% nme[y.indices])]
y = matrix(value[which(nme %in% nme[y.indices])],nrow=8)




I=16
J=8
c=2
T=length(nme)/25
A = matrix(,nrow=J-1,ncol=I)
rownames(A) = c("f","s","l","c","f","s","l")
A[1,]=c(rep(1,4),rep(0,12))
A[2,]=c(rep(0,4),rep(1,4),rep(0,8))
A[3,]=c(rep(0,8),rep(1,4),rep(0,4))
A[4,]=c(rep(0,12),rep(1,4))
A[5,]=rep(c(1,0,0,0),4)
A[6,]=rep(c(0,1,0,0),4)
A[7,]=rep(c(0,0,1,0),4)



#Step 2
smoothed_EM = function(data,c,A,V.init,eta.init)
{
  I=ncol(A)
  J=nrow(A)+1
  c=2
  
  #Step 1 (Set initial parameters)
  h=5
  w=2*h+1
  t=h+1
  eta = matrix(rep(eta.init,17),nrow=17,ncol=1)
  V = diag(rep(V.init,17))
  f.array = matrix(,ncol=0,nrow=I+1)


  #Step 2
  for(t in 6:(T-5))
  {
    eta.hat.old = eta[,ncol(eta)]
    eta.hat.new = c()
    
    #E-step
    #set lambda and phi based on eta
    lambda = exp(head(eta.hat.old,length(eta.hat.old)-1))
    phi = exp(tail(eta.hat.old,1))
    
    Sigma.full = diag(c(phi*lambda^c, phi))
    Sigma = diag(phi*lambda^c)
    min.t = t-h
    max.t = t+h
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-8,i] - A %*% lambda))
    R = Sigma - Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% A %*% Sigma
    
    Q.local = -w/2*(log(det(Sigma))) + tr(ginv(Sigma) %*% R) - 1/2*sum(sapply(1:11,function(i) t(m[,i]-lambda) %*% ginv(Sigma) %*% (m[,i]-lambda)))  
    Q.prior = function(eta.val)
    {
      eta.hat.old
      Sigma.hat = Sigma.full + V
      dmvnorm(eta.val, mean=eta.hat.old, sigma=Sigma.hat, log=TRUE)
    }
    g = function(eta.val){Q.local + Q.prior(eta.val)}
    
    #Find mode of g
    #Optimizing over all eta's at once doesn't work, so we do it termwise
    g.marginal = function(eta.change, eta.val,index)
    {
      #eta.change: eta coordinate we change
      #index: index for eta.change
      #eta.val: vector of etas
      eta.val[index]=eta.change
      g(eta.val)
    }
    
    eta.hat.new=eta.hat.old
    eta.hat.new[1] = optim(1, g.marginal, eta.val=eta.hat.old, index=1, control=list(fnscale=-1),method="Brent",lower=1.0001,upper=1e7)$par
    for(reps in 1:20)
    {
      for(i in 1:17)
      {
        eta.hat.new[i] = optim(eta.hat.new[i], g.marginal, eta.val=eta.hat.new, index=i, control=list(fnscale=-1),method="Brent",lower=1.0001,upper=1e7)$par
      }
    }
    eta.hat.old = eta.hat.new
    eta = cbind(eta,eta.hat.new)
    f.array = cbind(f.array,exp(dlogprior(head(eta.hat.new,length(eta.hat.new)-1),head(eta.hat.old,length(eta.hat.old)-1),Sigma) + dloglikelihood(t(data[-8,t]),A,lambda,Sigma)))
    print(paste("Step",t))
  }
  list("f"=f.array,"eta"=eta, "theta"=exp(eta))
}

f_em = smoothed_EM(y,c,A,V.init=100,eta.init=11)


#T-10 columns since we ignore the first 5 and last 5 (window size)
#ignore the first entry of theta (the -1 term) since that is our initial value.
#To find this we look at the matrix A and see which plots correspond to columns of A
theta.MA = matrix(,nrow=25,ncol=T-10)
theta.MA[6:9,] = f_em$theta[13:16,-1]
theta.MA[11:14,] = f_em$theta[9:12,-1]
theta.MA[16:19,] = f_em$theta[5:8,-1]
theta.MA[21:24,] = f_em$theta[1:4,-1]

#fill in other boxes with appropriate sums
theta.MA[1,] = apply(theta.MA[1+seq(5,20,5),],2,sum)
theta.MA[2,] = apply(theta.MA[2+seq(5,20,5),],2,sum)
theta.MA[3,] = apply(theta.MA[3+seq(5,20,5),],2,sum)
theta.MA[4,] = apply(theta.MA[4+seq(5,20,5),],2,sum)

theta.MA[10,] = apply(theta.MA[10-4:1,],2,sum)
theta.MA[15,] = apply(theta.MA[15-4:1,],2,sum)
theta.MA[20,] = apply(theta.MA[20-4:1,],2,sum)
theta.MA[25,] = apply(theta.MA[25-4:1,],2,sum)

theta.MA[5,] = apply(theta.MA[c(6:9,11:14,16:19,21:24),],2,sum)




dev.off()
par(mfrow=c(5,5),oma=c(1.5,1.5,1.5,1.5))
par(mar=c(2.5,2.5,1,0))
name.list = c("dst fddi", "dst switch", "dst local", "dst corp", "total", 
              "corp->fddi", "corp->switch", "corp->local", "corp->corp", "src corp",
              "local->fddi", "local->switch", "local->local", "local->corp", "src local",
              "switch->fddi", "switch->switch", "switch->local", "switch->corp", "src switch",
              "fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp", "src fddi")
plotted.values = matrix(,nrow=25,ncol=T)
for(i in 1:25)
{
  indices = which(nme == name.list[i])
  plotted.values[i,]=value[indices]
  if(i==21) #bottom left corner
  {
    plot(plotted.values[i,],type='l',ylim=range(0,1e6),main=name.list[i],lwd=2,col="cyan")
  }else if(i>20) #bottom row
  {
    plot(plotted.values[i,],type='l',ylim=range(0,1e6),yaxt="n",main=name.list[i],lwd=2,col="cyan")
  }else if(i%%5==1) #left side
  {
    plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",main=name.list[i],lwd=2,col="cyan")
  }else #everything else (no axis labels)
  {
    plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",yaxt="n",main=name.list[i],lwd=2,col="cyan")
  }
  #Plot fitted data
  lines(1:277+5,theta.MA[i,],lwd=2)
}

mtext("Smoothed Observation",outer=TRUE)
mtext("Bytes/sec",outer=TRUE,side=2)
mtext("Hour of Day",outer=TRUE,side=1)












