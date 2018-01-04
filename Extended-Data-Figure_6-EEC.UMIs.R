library(rsvd)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)
library(KernSmooth)
library(beeswarm)
library(stringr)
library(reshape2)
library(pvclust)
library(NMF)

source('Fxns.R')
data_name<-'GSE92332_EEC_UMIcounts.txt.gz'
EEC_UMIs<-load_data(data_name)

###  whether need this step??
gene_expr<-as.numeric(apply(EEC_UMIs,1,sum))
EEC_UMIs<-EEC_UMIs[which(gene_expr>0),]
###

#  select genes
v=get.variable.genes(EEC_UMIs,min.cv2=100)
var.genes=as.character(rownames(v)[v$p.adj<0.05])
EEC_tpm=data.frame(log2(1+tpm(EEC_UMIs)))

sample.names<-colnames(EEC_tpm)
region.names<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][3])))
cell.type.names<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][4])))
all.genes<-rownames(EEC_tpm)
# EEC.pv<-pvclust(as.data.frame(t(EEC_tpm[var.genes,])),method.hclust = 'ward.D2',method.dist = 'correlation',
#                 nboot = 10000,parallel = TRUE,quiet = FALSE)
# save(EEC.pv,file='EEC_pvclust.RData')
load('EEC_pvclust.RData')



###  Figure b
EEC.pca<-read.table('EEC_pca_scores.txt')[,1:11]

Heatmap_corr<-function(data,condition,all.condition){
  cat(sprintf('There ara %d conditions\n',length(condition)))
  tpm<-data.frame()
  for(i in 1:length(condition)){
    tpm<-rbind(tpm,t(data[,all.condition%in%condition[i]]))
  }
  #cat(sprintf('Whether creat data accurate %d \n',sum(dim(tpm.data)[1]==dim(tpm)[1])))
  tpm<-data.frame(t(tpm))
  cat(sprintf('Whether creat data accurate %d \n',sum(dim(data)[1]==dim(tpm)[2])))
  ### create Condition
  Condition<-c()
  for(i in 1:length(condition)){
    Condition<-c(Condition,rep(condition[i],sum(all.condition%in%condition[i])))
  }
  
  return(list(Condition,tpm))
}

EEC.pca.sort<-Heatmap_corr(data=t(EEC.pca),condition = unique(cell.type.names),all.condition = cell.type.names)
EEC.corr<-cor(EEC.pca.sort[[2]],method = 'pearson')

aheatmap(EEC.corr,Colv = NA,Rowv = NA,annCol = EEC.pca.sort[[1]],
          annRow = EEC.pca.sort[[1]])

###   Figure d
genes.d<-c("Vipr1","Ffar2","Gper1","Gpr119","Gpbar1","Ptger3","Galr1",
           "Cnr1","Ffar4","Adgrd1","Ffar1")

EEC.heatmap.d<-Heatmap_fun(genes = genes.d,tpm.data = EEC_tpm,
                           condition = unique(cell.type.names),all.condition = cell.type.names)

aheatmap(EEC.heatmap.d[[2]],Colv = NA,Rowv = NA,annCol = EEC.heatmap.d[[1]])

