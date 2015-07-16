dlogprior = function(eta,eta.hat,Sigma.hat)
{
  dmvnorm(eta, mean=eta.hat, sigma=Sigma.hat, log=TRUE)
}
dloglikelihood = function(Y,A,lambda,Sigma)
{
  dmvnorm(Y, mean=A %*% lambda, sigma= A %*% Sigma %*% t(A), log=TRUE)
}

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
  T=ncol(y)
  #initialize the theta vector for simplicity, keep all lambdas the same. 
  #choice of these parameters are important. I initially had small lambdas, and
  #the EM would not work.
  lambda = rep(lambda.init,ncol(A))
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
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-nrow(data),i] - A %*% lambda))
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


################
#SMOOTHED MODEL#
################

smoothed_EM = function(data,c,A,V.init=5,eta.init=3)
{
  I=ncol(A)
  J=nrow(A)+1
  c=2
  
  #Step 1 (Set initial parameters)
  h=5
  w=2*h+1
  t=h+1
  eta = matrix(rep(eta.init,ncol(A)+1),nrow=ncol(A)+1,ncol=1)
  V = diag(rep(V.init, ncol(A)+1))
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
    m = sapply(min.t:max.t, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (data[-nrow(data),i] - A %*% lambda))
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
    f.array = cbind(f.array,exp(dlogprior(head(eta.hat.new,length(eta.hat.new)-1),head(eta.hat.old,length(eta.hat.old)-1),Sigma) + dloglikelihood(t(data[-nrow(data),t]),A,lambda,Sigma)))
    print(paste("Step",t))
  }
  list("f"=f.array,"eta"=eta, "theta"=exp(eta))
}





