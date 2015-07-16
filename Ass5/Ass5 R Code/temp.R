dat = read.csv("1router_allcount.dat")
attach(dat)

nme
levels(nme)
x.corp=c("corp->corp","fddi->corp","local->corp","switch->corp")
x.fddi=c("corp->fddi","fddi->fddi","local->fddi","switch->fddi")
x.local=c("corp->local","fddi->local","local->local","switch->local")
x.switch=c("corp->switch","fddi->switch","local->switch","switch->switch")
x.src=c("src corp", "src fddi", "src local", "src switch")

corp.index = which(nme %in% x.corp)
fddi.index = which(nme %in% x.fddi)
local.index = which(nme %in% x.local)
switch.index = which(nme %in% x.switch)
src.index = which(nme %in% x.src)

plot(corp.index,value[corp.index],type='l')

plot(time[corp.index],value[corp.index],type='l')






















