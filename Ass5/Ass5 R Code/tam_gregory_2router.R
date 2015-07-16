library(MASS) #ginv
dat = read.csv("2router_linkcount.dat")
attach(dat)
source("tam_gregory_functions.R")

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

#aggregated ftest: Since ftest is a vector, we take the sum of square since we want to minimize it
f_em = locally_iid_EM(y,c,A,lambda.init=1e4)
f_em$f


f_em2 = smoothed_EM(y,c,A,V.init=5,eta.init=3)
f_em2$f



