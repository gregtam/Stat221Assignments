#1.5
impala = read.table("impala.txt", header = TRUE)
impala = impala[,1]

waterbuck = read.table("waterbuck.txt", header = TRUE)
waterbuck = waterbuck[,1]

unscaled.marg.post.N = function(N,y)
{
  S = sum(y)
  n = length(y)
  1/choose(n*N,S) * 1/(n*N+1) * prod(choose(N,y)/N)
}

y = waterbuck
N.seq = seq(max(y),max(y)+200,1)
marg.post.seq = sapply(N.seq, function(temp) unscaled.marg.post.N(temp,y))
waterbuck.constant = 1/sum(marg.post.seq)

y = impala
N.seq = seq(max(y),max(y)+1000,1)
marg.post.seq = sapply(N.seq, function(temp) unscaled.marg.post.N(temp,y))
impala.constant = 1/sum(marg.post.seq)

#1.6
marg.post.N.waterbuck = function(N,y)
{
  S = sum(y)
  n = length(y)
  val = 1/choose(n*N,S) * 1/(n*N+1) * prod(choose(N,y)/N)*waterbuck.constant
  if(is.na(val))
    return(0)
  val
}
marg.post.N.impala = function(N,y)
{
  S = sum(y)
  n = length(y)
  val = 1/choose(n*N,S) * 1/(n*N+1) * prod(choose(N,y)/N)*impala.constant
  if(is.na(val))
    return(0)
  val
}

y = waterbuck
p.waterbuck = sum(sapply(101:7000,function(i) marg.post.N.waterbuck(i,y)))
p.waterbuck

y = impala
p.impala = sum(sapply(101:7000, function(i) marg.post.N.impala(i,y)))
p.impala










