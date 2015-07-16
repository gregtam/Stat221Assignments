library(MASS) #ginv
dat = read.csv("2router_linkcount.dat")
attach(dat)

start_time=Sys.time()

task.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
job.id = as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID"))

y = matrix(value[which(nme %in% nme[1:16])],nrow=16)



#set A matrix
I=64
J=16
c=2
T=length(nme)/16
A = matrix(,nrow=J-1,ncol=I)
rownames(A) = c("router5","r4-local","switch","r4-others","gw1","gw2","gw3","gw-others","router5","r4-local","switch","r4-others","gw1","gw2","gw3")
A[1,]=c(rep(1,8),rep(0,56))
A[2,]=c(rep(0,8),rep(1,8),rep(0,48))
A[3,]=c(rep(0,16),rep(1,8),rep(0,40))
A[4,]=c(rep(0,24),rep(1,8),rep(0,32))
A[5,]=c(rep(0,32),rep(1,8),rep(0,24))
A[6,]=c(rep(0,40),rep(1,8),rep(0,16))
A[7,]=c(rep(0,48),rep(1,8),rep(0,8))
A[8,]=c(rep(0,56),rep(1,8))
A[9,]=rep(c(1,0,0,0,0,0,0,0),4)
A[10,]=rep(c(0,1,0,0,0,0,0,0),4)
A[11,]=rep(c(0,0,1,0,0,0,0,0),4)
A[12,]=rep(c(0,0,0,1,0,0,0,0),4)
A[13,]=rep(c(0,0,0,0,1,0,0,0),4)
A[14,]=rep(c(0,0,0,0,0,1,0,0),4)
A[15,]=rep(c(0,0,0,0,0,0,1,0),4)

ftest = function(theta,a,b,c)
{
  #the function f which we try to set to 0
  #takes theta as input, returns f(theta) vector
  ans = c()
  lambda = head(theta,length(theta)-1)
  phi = tail(theta,1)
  for(i in 1:length(theta))
  {
    if(i<=I)
      ans[i] = c*phi*lambda[i]^c + (2-c)*lambda[i]^2 - 2*(1-c)*lambda[i]*b[i] - c*a[i]
    if(i==I+1)
      ans[i] = sum(lambda^(-c+1)*(lambda - b))
  }
  ans
}

#############
#LOCAL MODEL#
#############

locally_iid_EM = function(data, c, A, lambda.init = 1e6, phi.init=0.5)
{
  w=11
  h=(w-1)/2
  T=length(nme)/16
  #initialize the theta vector for simplicity, keep all lambdas the same. 
  #choice of these parameters are important. I initially had small lambdas, and
  #the EM would not work.
  lambda = rep(lambda.init,64)
  phi = phi.init
  theta = matrix(c(lambda,phi),ncol=1)
  
  f.array = matrix(,ncol=0,nrow=I+1)
  #replicate for each t
  for(t in 6:(T-5))
  {
    #previous theta vector (lambda, phi)
    theta.old = theta[,ncol(theta)]
    theta.new=c()
    
    #E-step
    Sigma = diag(phi*lambda^c)
    min.t = t-h
    max.t = t+h
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-16,i] - A %*% lambda))
    R = Sigma - Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% A %*% Sigma
    
    #finds f, which is the derivative of Q(theta,theta^{(t)})
    #Newton-Raphson is done on f to find the maximum of Q   
    a=c()
    b=c()
    f=c()
    for(i in 1:I)
    {
      a[i] = R[i,i] + mean(m[i,]^2)
      b[i] = mean(m[i,])
      f[i] = c*phi*lambda[i]^c + (2-c)*lambda[i]^2 - 2*(1-c)*lambda[i]*b[i] - c*a[i]
    }
    f[I+1] = sum(lambda^(-c+1)*(lambda - b))
    
    
    #analytical solution: when c=2, it sets the lambda such that f()=0
    lambda = (-b + sqrt(b^2 + 4*a*phi))/(2*phi)
    #replace lambdas with analytical solutions
    theta.new[1:I] = lambda
    #check that f is close to 0
    ftest(c(lambda,phi),a,b,c)
    
    #M-step
    
    #Create Fdot matrix (derivative)
    Fdot = matrix(,nrow=I+1,ncol=I+1)
    for(i in 1:I)
      for(j in 1:I)
        Fdot[i,j] = (i==j)*(phi*c^2*lambda[i]^(c-1) + 2*(2-c)*lambda[i] - 2*(1-c)*b[i])
    for(j in 1:I)
      Fdot[I+1,j] = (2-c)*lambda[j]^(1-c) - (1-c)*lambda[j]^(-c)*b[j]
    for(i in 1:I)
      Fdot[i,I+1] = c*lambda[i]^c
    Fdot[I+1,I+1] = 0
    
    mult=1
    tempInv = ginv(Fdot) %*% f
    repeat
    {
      #Fractional Newton-Raphson: Shrinks amount it takes off to ensure phi remains greater than 0
      phi = theta.old[I+1] - mult*tempInv[I+1,]
      if(phi>0)
        break
      mult = mult/2
    }
    theta.new[I+1] = phi
    ftest(c(lambda,phi),a,b,c)
    
    theta.old = theta.new
    theta = cbind(theta,theta.new)
    f.array = cbind(f.array,ftest(theta.new,a,b,c))
  }
  list("f"=f.array,"theta"=theta,"a"=a,"b"=b)
}

#aggregated ftest: Since ftest is a vector, we take the sum of square since we want to minimize it
optim_ftest = function(theta,a,b,c){sum(ftest(theta,a,b,c)^2)}
Sys.time()
f_em = locally_iid_EM(y,c,A,lambda.init=temp)
f_em$f


end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")
save(list = c("f_em","start_time","end_time","total_time","last_par"), file=sprintf("odyssey/pset5/converge_task%d_job%d.rda", task.id, job.id))


#T-10 columns since we ignore the first 5 and last 5 (window size)
#ignore the first entry of theta (the -1 term) since that is our initial value.
#To find this we look at the matrix A and see which plots correspond to columns of A
theta.MA = matrix(,nrow=81,ncol=T-10)

theta.MA[10:17,] = f_em$theta[13:16,-1]
theta.MA[9+10:17,] = f_em$theta[9:12,-1]
theta.MA[9*2+10:17,] = f_em$theta[5:8,-1]
theta.MA[9*3+10:17,] = f_em$theta[1:4,-1]
theta.MA[9*4+10:17,] = f_em$theta[13:16,-1]
theta.MA[9*5+10:17,] = f_em$theta[9:12,-1]
theta.MA[9*6+10:17,] = f_em$theta[5:8,-1]
theta.MA[9*7+10:17,] = f_em$theta[1:4,-1]

#fill in other boxes with appropriate sums
theta.MA[1,] = apply(theta.MA[1+seq(9,72,9),],2,sum)
theta.MA[2,] = apply(theta.MA[2+seq(9,72,9),],2,sum)
theta.MA[3,] = apply(theta.MA[3+seq(9,72,9),],2,sum)
theta.MA[4,] = apply(theta.MA[4+seq(9,72,9),],2,sum)
theta.MA[5,] = apply(theta.MA[5+seq(9,72,9),],2,sum)
theta.MA[6,] = apply(theta.MA[6+seq(9,72,9),],2,sum)
theta.MA[7,] = apply(theta.MA[7+seq(9,72,9),],2,sum)
theta.MA[8,] = apply(theta.MA[8+seq(9,72,9),],2,sum)

theta.MA[18,] = apply(theta.MA[18-8:1,],2,sum)
theta.MA[27,] = apply(theta.MA[27-8:1,],2,sum)
theta.MA[36,] = apply(theta.MA[36-8:1,],2,sum)
theta.MA[45,] = apply(theta.MA[45-8:1,],2,sum)
theta.MA[54,] = apply(theta.MA[54-8:1,],2,sum)
theta.MA[63,] = apply(theta.MA[63-8:1,],2,sum)
theta.MA[72,] = apply(theta.MA[72-8:1,],2,sum)
theta.MA[81,] = apply(theta.MA[80-8:1,],2,sum)

theta.MA[9,] = apply(theta.MA[c(10:17,19:26,28:35,37:44,46:53,55:62,64:71,73:79),],2,sum)


#Plot actual data
dev.off()
pdf("tam_gregory_fig6_2routerlocal.pdf",width=9,height=9)
par(mfrow=c(9,9),oma=c(1.5,1.5,1.5,1.5))
par(mar=c(2.5,2.5,1,0))
name.list = c("dst router5", "dst r4-local", "dst switch", "dst r4-others","dst gw1", "dst gw2", "dst gw3", "dst gw-others", 
              "ori router5", "ori r4-local", "ori switch", "ori r4-others","ori gw1", "ori gw2", "ori gw3", "ori gw-others")
plotted.values = matrix(,nrow=81,ncol=T)

#stores indices of destimation and origin values so we know for which values of i to plot
index_vec = c(1:8,seq(18,81,9))
for(i in 1:81)
{
  if(i %in% index_vec)
  {
    indices = which(nme == name.list[which(index_vec==i)])
    plotted.values[i,]=value[indices]
    if(i==1)
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }else if(i==81)
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),yaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }else
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",yaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }
  }else if(i==73) #bottom left corner
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),main=name.list[i],lwd=2,col="cyan")
  }else if(i>73) #bottom row
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),yaxt="n",lwd=2,col="cyan")
  }else if(i%%9==1) #left side
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),xaxt="n",lwd=2,col="cyan")
  }else #everything else (no axis labels)
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),xaxt="n",yaxt="n",lwd=2,col="cyan")
  }
  #Plot fitted data
  lines(1:278+5,theta.MA[i,],lwd=2)
}
mtext("Local Observation",outer=TRUE)
mtext("Bytes/sec",outer=TRUE,side=2)
mtext("Hour of Day",outer=TRUE,side=1)
dev.off()


















################
#SMOOTHED MODEL#
################

smoothed_EM = function(data,c,A,V.init,eta.init)
{
  I=ncol(A)
  J=nrow(A)+1
  c=2
  
  #Step 1 (Set initial parameters)
  h=5
  w=2*h+1
  t=h+1
  eta = matrix(rep(eta.init,65),nrow=65,ncol=1)
  V = diag(rep(V.init,65))
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
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-16,i] - A %*% lambda))
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





f_em2 = smoothed_EM(y,c,A,V.init=5,eta.init=3)
f_em2$f


end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")
save(list = c("f_em","start_time","end_time","total_time","last_par"), file=sprintf("odyssey/pset5/converge_task%d_job%d.rda", task.id, job.id))


#T-10 columns since we ignore the first 5 and last 5 (window size)
#ignore the first entry of theta (the -1 term) since that is our initial value.
#To find this we look at the matrix A and see which plots correspond to columns of A
theta.MA = matrix(,nrow=81,ncol=T-10)

theta.MA[10:17,] = f_em2$theta[13:16,-1]
theta.MA[9+10:17,] = f_em2$theta[9:12,-1]
theta.MA[9*2+10:17,] = f_em2$theta[5:8,-1]
theta.MA[9*3+10:17,] = f_em2$theta[1:4,-1]
theta.MA[9*4+10:17,] = f_em2$theta[13:16,-1]
theta.MA[9*5+10:17,] = f_em2$theta[9:12,-1]
theta.MA[9*6+10:17,] = f_em2$theta[5:8,-1]
theta.MA[9*7+10:17,] = f_em2$theta[1:4,-1]

#fill in other boxes with appropriate sums
theta.MA[1,] = apply(theta.MA[1+seq(9,72,9),],2,sum)
theta.MA[2,] = apply(theta.MA[2+seq(9,72,9),],2,sum)
theta.MA[3,] = apply(theta.MA[3+seq(9,72,9),],2,sum)
theta.MA[4,] = apply(theta.MA[4+seq(9,72,9),],2,sum)
theta.MA[5,] = apply(theta.MA[5+seq(9,72,9),],2,sum)
theta.MA[6,] = apply(theta.MA[6+seq(9,72,9),],2,sum)
theta.MA[7,] = apply(theta.MA[7+seq(9,72,9),],2,sum)
theta.MA[8,] = apply(theta.MA[8+seq(9,72,9),],2,sum)

theta.MA[18,] = apply(theta.MA[18-8:1,],2,sum)
theta.MA[27,] = apply(theta.MA[27-8:1,],2,sum)
theta.MA[36,] = apply(theta.MA[36-8:1,],2,sum)
theta.MA[45,] = apply(theta.MA[45-8:1,],2,sum)
theta.MA[54,] = apply(theta.MA[54-8:1,],2,sum)
theta.MA[63,] = apply(theta.MA[63-8:1,],2,sum)
theta.MA[72,] = apply(theta.MA[72-8:1,],2,sum)
theta.MA[81,] = apply(theta.MA[80-8:1,],2,sum)

theta.MA[9,] = apply(theta.MA[c(10:17,19:26,28:35,37:44,46:53,55:62,64:71,73:79),],2,sum)


#Plot actual data
dev.off()
pdf("tam_gregory_fig6_2router.pdf",width=9,height=9)
par(mfrow=c(9,9),oma=c(1.5,1.5,1.5,1.5))
par(mar=c(2.5,2.5,1,0))
name.list = c("dst router5", "dst r4-local", "dst switch", "dst r4-others","dst gw1", "dst gw2", "dst gw3", "dst gw-others", 
              "ori router5", "ori r4-local", "ori switch", "ori r4-others","ori gw1", "ori gw2", "ori gw3", "ori gw-others")
plotted.values = matrix(,nrow=81,ncol=T)

#stores indices of destimation and origin values so we know for which values of i to plot
index_vec = c(1:8,seq(18,81,9))
for(i in 1:81)
{
  if(i %in% index_vec)
  {
    indices = which(nme == name.list[which(index_vec==i)])
    plotted.values[i,]=value[indices]
    if(i==1)
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }else if(i==81)
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),yaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }else
    {
      plot(plotted.values[i,],type='l',ylim=range(0,1e6),xaxt="n",yaxt="n",main=name.list[which(index_vec==i)],lwd=2,col="cyan")
    }
  }else if(i==73) #bottom left corner
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),main=name.list[i],lwd=2,col="cyan")
  }else if(i>73) #bottom row
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),yaxt="n",lwd=2,col="cyan")
  }else if(i%%9==1) #left side
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),xaxt="n",lwd=2,col="cyan")
  }else #everything else (no axis labels)
  {
    plot(-1,type='l',xlim=range(1,288),ylim=range(0,1e6),xaxt="n",yaxt="n",lwd=2,col="cyan")
  }
  #Plot fitted data
  lines(1:278+5,theta.MA[i,],lwd=2)
}
mtext("Smoothed Observation",outer=TRUE)
mtext("Bytes/sec",outer=TRUE,side=2)
mtext("Hour of Day",outer=TRUE,side=1)
dev.off()

