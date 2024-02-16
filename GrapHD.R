GrapHD=function(G){
  require(igraph)
  require(future.apply)
  v=V(G)
  #aM=as_adjacency_matrix(G)
  plan("multisession", workers=15)
  d=future_sapply(v, function(vi){
    vv=c(1)
    i=1
    repeat {
      ev=unique(unlist(names(ego(G, i, nodes=vi)[[1]])))
      vv=c(vv,length(ev))
      i=i+1
      if (length(ev)>(0.8*vcount(G))) break
    }
    x=log(seq_along(vv))
    y=log(vv)
    k=sapply(3:length(vv),function(j) lm(y[1:j]~x[1:j])$coefficients[2])
    k=max(k)
    return(k)
  })
  closeAllConnections()
  hist(log10(d), breaks=length(d)/10)
  return(c(mean=10^(mean(log10(d))), max=max(d)))
}
