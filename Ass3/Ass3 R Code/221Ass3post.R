#2a
load("2a.rda")
plot(data.sgd$Xplot, data.sgd$Yplot,type='l',ylim=range(-6,2),xlab = "training size t", ylab = "excess risk",axes=FALSE)
lines(data.asgd.good$Xplot, data.asgd.good$Yplot,col="red")
lines(data.asgd.bad$Xplot, data.asgd.bad$Yplot, col="blue")
lines(data.implicit$Xplot, data.implicit$Yplot, col="green")

png("221Ass3-Q2-a.png",width=600,height=400,res=200)
plot(log10(1:length(data.sgd$Yplot)), data.sgd$Yplot,type='l',ylim=range(-6,2),xlab = "training size t", ylab = "excess risk",axes=FALSE)
lines(log10(1:length(data.asgd.good$Yplot)), data.asgd.good$Yplot,col="red")
lines(log10(1:length(data.asgd.bad$Yplot)), data.asgd.bad$Yplot, col="blue")
lines(log10(1:length(data.implicit$Yplot)), data.implicit$Yplot, col="green")



legend("bottomleft", c("sgd","asgd.good","asgd.bad","implicit"), 
       lwd = rep(2,4), col=c("black","red","blue","maroon"), title = "Legend",cex=0.9)
#x axis
axis(1,0:6,c("10e0","10e1","10e2","10e3","10e4","10e5","10e6"))
#y axis
axis(2,-6:2,c("10e-6","10e-5","10e-4","10e-3","10e-2","10e-1","10e0","10e1","10e2"))
dev.off()


#2b
plot.risk <- function(data, est,col)
{
  # est = p x niters 
  est.bias = apply(est, 2, function(colum) 
    log10(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))
  plot(log10(1:length(est.bias)),est.bias, type="l", lty=3,xlim=range(2,5),ylim=range(-4,4),col=col,axes=FALSE,xlab = "training size t", ylab="excess risk")
}
lines.risk <- function(data, est,col)
{
  # est = p x niters 
  est.bias = apply(est, 2, function(colum) 
    log10(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))
  lines(log10(1:length(est.bias)),est.bias, type="l", lty=3,xlim=range(2,5),ylim=range(-4,4),col=col)
}
for(task.num in 1:4)
{
  job.id = 22729420
  job.id = 22805866
  job.id = 22814597
  load(paste("2b_job",job.id,"_task",task.num,".rda",sep=""))
  colours = rainbow(4)
  
  A = data$A
  theta = descent_values
  dim.t = length(theta[1,])
  excess_risk = sapply(1:dim.t, function(i) theta[,i] %*% A %*% theta[,i])
  if(task.num==1)
  {
    plot.risk(data,descent_values, col=colours[task.num])
  }else
  {
    lines.risk(data,descent_values,col = colours[task.num])
  }
}

legend("bottomleft", c("sgd","asgd","implicit","batch"), 
       lwd = rep(2,4), col=colours, title = "Legend",cex=0.9)
#x axis
axis(1,0:5,c("10e0","10e1","10e2","10e3","10e4","10e5"))
#y axis
axis(2,-4:4,c("10e-4","10e-3","10e-2","10e-1","10e0","10e1","10e2","10e3","10e4"))
dev.off()







#2c
job.id = 22738468
bias_mat_sgd = matrix(,nrow=0,ncol=10000)
var_mat_sgd = matrix(,nrow=0,ncol=10000)
for(task.num in 1:10)
{
  load(paste("2c_job",job.id,"_task",task.num,".rda",sep=""))
  bias_mat_sgd = rbind(bias_mat_sgd, biases)
  var_mat_sgd = rbind(var_mat_sgd,variances)
}
bias_mat_asgd = matrix(,nrow=0,ncol=10000)
var_mat_asgd = matrix(,nrow=0,ncol=10000)
for(task.num in 11:20)
{
  load(paste("2c_job",job.id,"_task",task.num,".rda",sep=""))
  bias_mat_asgd = rbind(bias_mat_asgd, biases)
  var_mat_asgd = rbind(var_mat_asgd,variances)
}
bias_mat_implicit = matrix(,nrow=0,ncol=10000)
var_mat_implicit = matrix(,nrow=0,ncol=10000)
for(task.num in 21:30)
{
  load(paste("2c_job",job.id,"_task",task.num,".rda",sep=""))
  bias_mat_implicit = rbind(bias_mat_implicit, biases)
  var_mat_implicit = rbind(var_mat_implicit,variances)
}

pdf("2cbias1.pdf",height=4,width=6)
plot(log10(1:10000), log10(bias_mat_sgd[1,]),type='l',ylim=range(0,0.5),col=rainbow(10)[1],main="SGD bias")
for(i in 2:10)
{
  lines(log10(1:10000), log10(bias_mat_sgd[i,]),col="red")
}
dev.off()

pdf("2cbias2.pdf",height=4,width=6)
plot(log10(1:10000), log10(bias_mat_asgd[1,]),type='l',ylim=range(0,0.5),col=rainbow(10)[1],main="ASGD bias")
for(i in 2:10)
{
  lines(log10(1:10000), log10(bias_mat_asgd[i,]),col=rainbow(10)[i])
}
dev.off()

pdf("2cbias3.pdf",height=4,width=6)
plot(log10(1:10000), log10(bias_mat_implicit[1,]),type='l',ylim=range(0,0.5),col=rainbow(10)[1],main="Implicit bias")
for(i in 2:10)
{
  lines(log10(1:10000), log10(bias_mat_implicit[i,]),col=rainbow(10)[i])
}
dev.off()

pdf("2cvar1.pdf",height=4,width=6)
plot(log10(1:10000), log10(var_mat_sgd[1,]),type='l',ylim=range(0,max(log10(var_mat_sgd))),col=rainbow(10)[0],main="SGD variance")
for(i in 2:10)
{
  lines(log10(1:10000), log10(var_mat_sgd[i,]),col=rainbow(10)[i])
}
dev.off()

pdf("2cvar2.pdf",height=4,width=6)
plot(log10(1:10000), log10(var_mat_asgd[1,]),type='l',ylim=range(-3,max(log10(var_mat_asgd))), col = rainbow(10)[1],main="ASGD variance")
for(i in 2:10)
{
  lines(log10(1:10000), log10(var_mat_asgd[i,]),col=rainbow(10)[i])
}
dev.off()

pdf("2cvar3.pdf",height=4,width=6)
plot(log10(1:10000), log10(var_mat_implicit[1,]),type='l',ylim=range(-1.5,max(log10(var_mat_implicit))),col=rainbow(10)[1],main="Implicit variance")
for(i in 2:10)
{
  lines(log10(1:10000), log10(var_mat_implicit[i,]),col=rainbow(10)[i])
}
dev.off()



#3c

timesmat = matrix(,nrow=0,ncol=6)
job.id = 22853657
for(task.num in 1:6)
{
  load(paste("3c_job",job.id,"_task",task.num,".rda",sep=""))
  timesmat = rbind(timesmat,times)
}


#3d?
job.id = 2
load(paste("3c_job",job.id,"_task",2,".rda",sep=""))










