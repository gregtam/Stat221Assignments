data = read.csv("a9a",sep=" ",header=FALSE)
data = data[,-16] #removes last column which is glitchy (all NA)

x = matrix(,nrow=32561,ncol=14)

for(cols in 1:14)
{
  #Extracts the number x from the string "x:1"
  x[,cols] = sapply(as.character(data[,cols+1]), function(i) as.numeric(strsplit(i,":")[[1]][1]))
}
y = data[,1]

#Split x and y into training and test sets
train.indices = sample(1:nrow(x),0.8*nrow(x))
x.train = x[train.indices,]
x.test = x[-train.indices,]
y.train = y[train.indices]
y.test = y[-train.indices]

#Creates design matrix from the matrix of indices x.train
x.mat.train = matrix(0,nrow=nrow(x.train), ncol=123)
for(i in 1:nrow(x.mat.train))
  x.mat.train[i,x.train[i,]]=1
x.mat.train = cbind(1,x.mat.train)

x.mat.test = matrix(0,nrow=nrow(x.test), ncol=123)
for(i in 1:nrow(x.mat.test))
  x.mat.test[i,x.test[i,]]=1
x.mat.test = cbind(1,x.mat.test)


#Set Parameters
N = nrow(x)
n = 10

time = 1:50
a = 0.015
b = 1
gamma = 0.55

epsilon = a*(b+time)^(-gamma)
epsilon[1]
epsilon[length(epsilon)]

expit = function(x){1/(1+exp(-x))}
beta.mat = matrix(0,nrow=124,ncol=1)

for(t in 1:50)
{
  for(batch.rep in seq(1,124,10))
  {
    n=10
    if(batch.rep==121)
      n=4
    beta.old = beta.mat[,ncol(beta.mat)]
    logprior.grad = sign(beta.old)
    ll.grad = function(index)
    {
      expit(y.train[index] * sum(beta.old * x.mat.train[index,])) * y.train[index] * x.mat.train[index,]
    }

    eta = rnorm(124,mean=0,sd=sqrt(epsilon[t]))
    ll.grad.vec = sapply(1:n,ll.grad)

    #proposed update
    beta.new = beta.old + epsilon[t]/2*(logprior.grad + N/n*apply(ll.grad.vec,1,sum)) + eta
    beta.mat = cbind(beta.mat,beta.new)
    beta.old = beta.new
  }
}

length(y.train)
dim(x.mat.train)
dim(beta.mat)

#indices when entire sweep of dataset has finished
beta.index = seq(14,13*50+1,13)
beta.final = beta.mat[,beta.index]

#Average log joint probability per data item vs number of sweeps
avg.log.prob = c()
for(i in 1:ncol(beta.final))
{
  temp = mean(expit(y.train * x.mat.train %*% beta.final[,i]))
  avg.log.prob = c(avg.log.prob, log(temp))
}

pdf("avglogprob.pdf",weight=6,height=6)
dev.off()


#Accuracy on test set vs number of swepes
length(y.test)
dim(x.mat.test)
dim(beta.mat)

#vector of accuracies
acc = c()
for(i in 1:ncol(beta.mat))
{
  probs = expit(y.test * x.mat.test %*% beta.mat[,i])
  temp = length(which(probs>0.5))/length(probs)
  acc = c(acc,temp)
}
pdf("fig3.pdf",width=8,height=4)
par(mfrow=c(1,2))
plot(avg.log.prob,type='l',main="Average Log Probability")
plot((1:28-1)/13, acc[1:28],type='l',xlab="Number of Sweeps",ylab="Accuracy",main="Accuracy vs Number of Sweeps")
dev.off()








