# packages ----
require(parallel)
require(matrixStats)
require(reshape2)
require(future.apply)
require(stringi)
require(imputeR)
require(FactoMineR)
require(factoextra)
require(rgl)
require(igraph)
require(uwot)
require(limma)

# Load data -----
pth=".//Peptide arrays"
fs=list.files(path=pth)
fs=grep(".csv",fs, value=T)
pepmap=read.csv(paste(pth,fs[1],sep="//"))
datal=lapply(fs[fs!="Mapptd.csv"],function(ff){
  read.csv(paste(pth,ff,sep="//"))
})
peps=table(unlist(lapply(datal, function(ff) ff[,1])))
table(peps)
peps=names(peps)

D=sapply(datal,function(f){
  x=as.numeric(f[,2])
  names(x)=f[,1]
  x[peps]
})
colnames(D)=stri_extract(fs[fs!="Mapptd.csv"], regex="\\w+")
rownames(D)=peps
iNA=which(is.na(D), arr.ind = T)
Di=impute(D, lmFun = "lassoR")$imp
dimnames(Di)=dimnames(D)
i0=which(rowSums(Di)==0)
Di=Di[-i0,]
Dil=log10(Di+0.001)

groups=substr(colnames(D),2,4)
gi=cbind(BIV=rowMeans(Dil[,groups=="BIV"]),
         BOC=rowMeans(Dil[,groups=="BOC"]),
         CIV=Dil[,groups=="CIV"],
         MIV=rowMeans(Dil[,groups=="MIV"]),
         uma=Dil[,groups=="uma"])

gipos=gi-min(gi)



pcaD=PCA(t(Dil))
plot3d(cmdscale(dist(t(Dil)), k=3), size=15, col=as.numeric(as.factor(groups)))
plot3d(pcaD$ind$coord[,1:3], size=15, col=as.numeric(as.factor(groups)))

# Rearaph ----

G=ReGr(Dil)
compG=components(G)
G1=induced.subgraph(G, vids=V(G)[compG$membership==1])
Lm=embed_laplacian_matrix(G1,no=vcount(G1)-1)
j=GrapHD(G1)
j=ceiling(j["max"])
j=(vcount(G1)-j):(vcount(G1)-1)
Mx=Lm$X[,j]
uMx=umap(Mx, n_neighbors = 100, min_dist = 0.005, spread=10,  n_epochs = 1500, 
         nn_method = "fnn", verbose=F, n_components = 3, init = "normlaplacian")

uMx2=umap(Mx, n_neighbors = 100, min_dist = 0.005, spread=10,  n_epochs = 1500, 
         nn_method = "fnn", verbose=F, n_components = 2, init = "normlaplacian")
G1=set.vertex.attribute(G1, name="x", value=uMx2[,1])
G1=set.vertex.attribute(G1, name="y", value=uMx2[,2])
for (i in 1:5) G1=set.vertex.attribute(G1, name=colnames(gi)[i], value=gipos[names(V(G1)),i])
write.graph(G1, format = "graphml", file="G1.graphml")
rglplot(G1, layout=uMx, vertex.size=0.5, edge.width=0.2, 
        vertex.label=NA, vbertex.color=rank(vertex_attr(G1)$uma))

clq=max_cliques(G1, min=30)
