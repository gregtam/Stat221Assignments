#initialize theta
lambda = rep(10,16)
phi = 10
theta = matrix(c(lambda,phi),ncol=1)
f.array = matrix(,ncol=0,nrow=I+1)

#set A matrix
A = matrix(,nrow=J-1,ncol=I)
A[1,]=c(rep(1,4),rep(0,12))
A[2,]=c(rep(0,4),rep(1,4),rep(0,8))
A[3,]=c(rep(0,8),rep(1,4),rep(0,4))
A[4,]=c(rep(0,12),rep(1,4))
A[5,]=rep(c(1,0,0,0),4)
A[6,]=rep(c(0,1,0,0),4)
A[7,]=rep(c(0,0,1,0),4)

root = matrix(,nrow=16,ncol=0)
#replicate
for(reps in 1:100)
{
  theta.old = theta[,ncol(theta)]
  Sigma = diag(phi*lambda^c)
  
  A %*% Sigma %*% t(A)
  solve(A %*% Sigma %*% t(A))
  # lambda + Sigma %*% t(A) %*% solve(A %*% Sigma %*% t(A)) * (y - A %*% lambda)
  m = sapply(1:T, function(i) lambda + Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% (y[-8,i] - A %*% lambda))
  R = Sigma - Sigma %*% t(A) %*% ginv(A %*% Sigma %*% t(A)) %*% A %*% Sigma
  
  a=c()
  b=c()
  f=c()
  for(i in 1:I)
  {
    a[i] = R[i,i] + mean(m[i,]^2)
    b[i] = mean(m[i,])
    f[i] = c*phi*lambda[i]^c + (2-c)*lambda[i]^2 - 2*(1-c)*lambda[i]*b[i] - c*a[i]
  }
  f[i+1] = sum(lambda^(-c+1)*(lambda - b))
  
  #analytical solution for lambda (c=2)
  lambda = (-b + sqrt(b^2 + 4*a*phi))/(2*phi)
  #check analytical solution is a root
  root = cbind(root,2*phi*lambda^2 + 2*lambda*b-2*a)
  
  Fdot = matrix(,nrow=I+1,ncol=I+1)
  for(i in 1:I)
    for(j in 1:I)
      Fdot[i,j] = (i==j)*(phi*c^2*lambda[i]^(c-1) + 2*(2-c)*lambda[i] - 2*(1-c)*b[i])
  for(j in 1:I)
    Fdot[I+1,j] = (2-c)*lambda[j]^(1-c) - (1-c)*lambda[j]^(-c)*b[j]
  for(i in 1:I)
    Fdot[i,I+1] = c*lambda[i]^c
  Fdot[I+1,I+1] = 0
  
  #   theta.new = theta.old - ginv(Fdot) %*% f
  theta.new[1:I] = lambda
  mult=1
  tempInv = ginv(Fdot) %*% f
  repeat
  {
    #Fractional Newton-Raphson: Shrinks amount it takes off so parameters are greater than 0
    phi = (theta.old - 1/mult*tempInv)[I+1,]
    if(phi>0)
      break
    mult = mult+1
  }
  
  theta.new[I+1] = phi
  #replace lambdas with analytical solutions
  theta = cbind(theta,theta.new)
  f.array = cbind(f.array,ftest(theta.new))
}

par(mfrow=c(2,1))
plot(f.array[1,50:ncol(theta)-1])
f.array[1,]
plot(f.array[17,-1])
f.array[17,]
theta[,100]
ftest(theta[,100])
ftest(theta.new)

