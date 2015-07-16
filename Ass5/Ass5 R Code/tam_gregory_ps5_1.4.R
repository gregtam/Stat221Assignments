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

#set A matrix
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
  T=length(nme)/25
  #initialize the theta vector for simplicity, keep all lambdas the same. 
  #choice of these parameters are important. I initially had small lambdas, and
  #the EM would not work.
  lambda = rep(lambda.init,16)
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
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-8,i] - A %*% lambda))
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
f_em = locally_iid_EM(y,c,A,lambda.init=1e4)


end_time = Sys.time()
total_time = as.numeric(end_time-start_time,units="mins")
save(list = c("f_em","start_time","end_time","total_time","last_par"), file=sprintf("odyssey/pset5/converge_task%d_job%d.rda", task.id, job.id))


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


#Plot actual data
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
mtext("Local Observation",outer=TRUE)
mtext("Bytes/sec",outer=TRUE,side=2)
mtext("Hour of Day",outer=TRUE,side=1)

