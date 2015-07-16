dat = read.csv("1router_allcount.dat")
attach(dat)

nme
levels(nme)
x.corp=c("corp->corp","fddi->corp","local->corp","switch->corp")
x.fddi=c("corp->fddi","fddi->fddi","local->fddi","switch->fddi")
x.local=c("corp->local","fddi->local","local->local","switch->local")
x.switch=c("corp->switch","fddi->switch","local->switch","switch->switch")
x.src=c("src corp", "src fddi", "src local", "src switch")


x.corp=c("corp->corp","corp->fddi","corp->local","corp->switch")
x.fddi=c("fddi->corp","fddi->fddi","fddi->local","fddi->switch")
x.local=c("local->corp","local->fddi","local->local","local->switch")
x.switch=c("switch->corp","switch->fddi","switch->local","switch->switch")
x.src=c("src corp", "src fddi", "src local", "src switch")

x.corp=c("src corp", "dst corp")
x.fddi=c("src fddi", "dst fddi")
x.local=c("src local", "dst local")
x.switch=c("src switch", "dst switch")


corp.index = which(nme %in% x.corp)
fddi.index = which(nme %in% x.fddi)
local.index = which(nme %in% x.local)
switch.index = which(nme %in% x.switch)

par(mfrow=c(4,1),oma=c(1.5,1.5,1.5,1.5))
par(mar=c(3,3,1,0))
corp.index = t(sapply(1:2,function(i) which(nme %in% x.corp[i])))
corp.val = matrix(,nrow=nrow(corp.index),ncol=ncol(corp.index))
corp.val[1,] = value[corp.index[1,]]
corp.val[2,] = value[corp.index[2,]]
plot(corp.index[1,],corp.val[1,],type='l',ylim=range(min(corp.val),max(corp.val)),main="corp")
lines(corp.index[2,],corp.val[2,])

local.index = t(sapply(1:2,function(i) which(nme %in% x.local[i])))
local.val = matrix(,nrow=nrow(local.index),ncol=ncol(local.index))
local.val[1,] = value[local.index[1,]]
local.val[2,] = value[local.index[2,]]
plot(local.index[1,],local.val[1,],type='l',ylim=range(min(local.val),max(local.val)),main="local")
lines(local.index[2,],local.val[2,])

switch.index = t(sapply(1:2,function(i) which(nme %in% x.switch[i])))
switch.val = matrix(,nrow=nrow(switch.index),ncol=ncol(switch.index))
switch.val[1,] = value[switch.index[1,]]
switch.val[2,] = value[switch.index[2,]]
plot(switch.index[1,],switch.val[1,],type='l',ylim=range(min(switch.val),max(switch.val)),main="switch")
lines(switch.index[2,],switch.val[2,])

fddi.index = t(sapply(1:2,function(i) which(nme %in% x.fddi[i])))
fddi.val = matrix(,nrow=nrow(fddi.index),ncol=ncol(fddi.index))
fddi.val[1,] = value[fddi.index[1,]]
fddi.val[2,] = value[fddi.index[2,]]
plot(fddi.index[1,],fddi.val[1,],type='l',ylim=range(min(fddi.val),max(fddi.val)),main="fddi")
lines(fddi.index[2,],fddi.val[2,])

mtext("Origin Destination",outer=TRUE)
mtext("Bytes/sec",outer=TRUE,side=2)
mtext("Hour of Day",outer=TRUE,side=1)































