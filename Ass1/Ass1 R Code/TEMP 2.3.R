# ll.G = function(par, G, index, theta, X)
# {
#   N = length(G)/2
#   P = length(theta)
#   log_likelihood_sum = 0
#   G[index]=par
#   for(n in 1:N) #for each observation
#   {
#     for(p in 1:P) #for each feature
#     {
#       val_H = theta[[p]]$high[X[-1,-1][n,p]]
#       val_L = theta[[p]]$low[X[-1,-1][n,p]]
#       if(G[n] * val_H + G[n + N/2] * val_L == 0)
#         log_likelihood_sum = log_likelihood_sum + 1e-16
#       else
#         log_likelihood_sum = log_likelihood_sum + log(G[n] * val_H + G[n + N/2] * val_L)
#     }
#   }
#   log_likelihood_sum
# }
# ll.G(1,G,1,theta0List,X)
# 
# opt[i] = optim(0.5, ll.G, control=list(fnscale=-1), method = "L-BFGS-B", theta = theta0List, X=X[-1,-1], G=G, index=1)
# 
# getHighVal = function(p,n)
# {
#   theta0List[[p]]$high[X[-1,-1][n,p]]
# }
# getLowVal = function(p,n)
# {
#   theta0List[[p]]$low[X[-1,-1][n,p]]
# }
# sapply(1:49, getHighVal,1:69)
# sapply(1:49, getLowVal,1:69)



ll.G = function(param, G, index, theta, X)
{
  N = length(G)/2
  P = length(theta)
  log_likelihood_sum = 0
  G[index] = param
  
  getHighVal = function(p,n)
  {
    theta[[p]]$high[Xmin[n,p]]
  }
  getLowVal = function(p,n)
  {
    theta[[p]]$low[Xmin[n,p]]
  }
  val_H = sapply(1:49, getHighVal,n=1:69)
  val_L = sapply(1:49, getLowVal,n=1:69)
  
  val_G = matrix(G[1:N], nrow=N,ncol=dim(val_H)[2])
  val_G2 = 1 - val_G
  
  mix = val_G * val_H + val_G2 * val_L
  sum(log(ifelse(mix>0, mix, mix+exp(1e-16))))
}
ll.G(1,G,1,theta0List,X)


ll.thetaL = function(param, G, index, theta, X)
{
  N = length(G)/2
  P = length(theta)
  log_likelihood_sum = 0
  #u_i
  scaled_param = exp(param)/(sum(exp(param)))
  theta[[index]]$low = scaled_param
  
  getHighVal = function(p,n)
  {
    theta[[p]]$high[Xmin[n,p]]
  }
  getLowVal = function(p,n)
  {
    theta[[p]]$low[Xmin[n,p]]
  }
  val_H = sapply(1:49, getHighVal,n=1:69)
  val_L = sapply(1:49, getLowVal,n=1:69)
  
  val_G = matrix(G[1:N], nrow=N,ncol=dim(val_H)[2])
  val_G2 = 1 - val_G
  
  mix = val_G * val_H + val_G2 * val_L
  sum(log(ifelse(mix>0, mix, mix+exp(1e-16))))
}
ll.thetaL(0.2,G,3,theta0List,X)



ll.thetaH = function(param, G, index, theta, X)
{
  N = length(G)/2
  P = length(theta)
  log_likelihood_sum = 0
  #u_i
  scaled_param = exp(param)/(sum(exp(param)))
  theta[[index]]$high = scaled_param
  
  getHighVal = function(p,n)
  {
    theta[[p]]$high[Xmin[n,p]]
  }
  getLowVal = function(p,n)
  {
    theta[[p]]$low[Xmin[n,p]]
  }
  val_H = sapply(1:49, getHighVal,n=1:69)
  val_L = sapply(1:49, getLowVal,n=1:69)
  
  val_G = matrix(G[1:N], nrow=N,ncol=dim(val_H)[2])
  val_G2 = 1 - val_G
  
  mix = val_G * val_H + val_G2 * val_L
  sum(log(ifelse(mix>0, mix, mix+exp(1e-16))))
}
ll.thetaH(0.3,G,3,theta0List,X)



Sys.time()
optTheta = theta0List
optG = G
likeli=c()
for(counter in 1:11)
{
  print(paste("Doing Count", counter, Sys.time()))
  #optimize G
  for(i in 1:69)
    optG[i] = optim(0.5, ll.G, control=list(fnscale=-1), method = "L-BFGS-B", lower=0, upper=1, theta = optTheta, X=X, G=optG, index=i)$par
  optG[,2] = 1-optG[,1]
  
  
  #optimize thetaL
  for(index in 1:length(theta))
  {
    #optimize thetaL
    init = rep(0,length(theta[[index]]$low))
    thetaL = optim(init, ll.thetaL, control=list(fnscale=-1), method = "L-BFGS-B", theta = optTheta,X=X,G=optG,index=index)$par
    optL = exp(thetaL)/(sum(exp(thetaL)))
    optL
    optTheta[[index]]$level = c(1: length(optL)-1)
    optTheta[[index]]$low = optL
  }
  for(index in 1:length(theta))
  {
    #optimize thetaH
    init = rep(0,length(theta[[index]]$high))
    thetaH = optim(init, ll.thetaH, control=list(fnscale=-1), method = "L-BFGS-B", theta = optTheta,X=X,G=optG,index=index)$par
    optH = exp(thetaH)/(sum(exp(thetaH)))
    optH
    optTheta[[index]]$high = optH
  }
  likeli[counter] = ll(optG, optTheta,X[-1,-1])
}
Sys.time()

ll(optG, optTheta, X[-1,-1])









#TEST
x = seq(0,1,0.01)
y = c()
for(i in x)
  y = c(y, ll.G(i, G, 1, theta0List, X))
plot(x,y,type='l')






