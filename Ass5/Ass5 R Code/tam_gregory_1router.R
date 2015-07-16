library(MASS) #ginv
dat = read.csv("1router_allcount.dat")
attach(dat)
source("tam_gregory_functions.R")

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


f_em = locally_iid_EM(y,c,A,lambda.init=1e4)



library(numDeriv) #hessian function
library(mvtnorm)
library(psych) #trace function
library(MASS) #ginnv


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


f_em = smoothed_EM(y,c,A,V.init=100,eta.init=11)









