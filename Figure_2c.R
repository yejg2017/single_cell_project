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
library(destiny)
library(RColorBrewer)
library(rgl)
palette(brewer.pal(6, 'Spectral'))

source('Fxns.R')

## Figure f
Regional_UMIs<-load_data("./data/GSE92332_Regional_UMIcounts.txt.gz")
Regional_UMIs<-Regional_UMIs[which(unlist(apply(Regional_UMIs,1,sum))>0),]
Regional_tpm<-data.frame(log2(1+tpm(Regional_UMIs)))
#region.var.genes<-get.variable.genes_cvdiff(Regional_tpm)

regional.v = get.variable.genes(Regional_UMIs, min.cv2 = 100)
regional.var.genes = as.character(rownames(regional.v)[regional.v$p.adj<0.05])
regional.genes<-rownames(Regional_tpm)
region_groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][2])))

cell.groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][4])))
region.sample.names<-colnames(Regional_tpm)



#   Figure 2.c
all.diffusion.matrix<-cbind(data.frame(Cell=region.sample.names[cell.groups%in%c('TA','Stem','EP','Enterocyte')]),
                            t(Regional_tpm[regional.var.genes,cell.groups%in%c('TA','Stem','EP','Enterocyte')]))
all.ct<-as.ExpressionSet(all.diffusion.matrix)
#phenoData(all.ct)$regions<-as.factor(region_groups[cell.groups%in%c('TA','Stem','EP','Enterocyte')])

all.dif<-DiffusionMap(all.ct,verbose = T)
save(all.dif,file='./data/Regional_UMIs_DiffusionMap_4.RData')
# load('./data/Regional_UMIs_DiffusionMap_4.RData')
plot.DiffusionMap(all.dif,dims = c(1,3,4),
                  col =as.factor(region_groups[cell.groups%in%c('TA','Stem','EP','Enterocyte')]),pch=20)


plot.DiffusionMap(all.dif,dims = c(1,3,4),
                  col =as.factor(cell.groups[cell.groups%in%c('TA','Stem','EP','Enterocyte')]),pch=20)


# rgl     :   failed  ??? 
# plot3d(eigenvectors(all.dif)[,c(1,4,3)],
#        col=as.integer(as.factor(region_groups[cell.groups%in%c('TA','Stem','EP','Enterocyte')])),pch=20)
# rgl.close()


DiffusionMap_Meantpm<-function(gene,DMap.obj,tpm.data,all.cells=cell.groups,
                               cells=c('TA','Stem','EP','Enterocyte'),do.plot=TRUE){
  
  #   gene : the gene use for caculate mean expression
  #   DMap.obj :  the object from DiffusionMap function
  #   tpm.data : the TPM data:gene must in the data gene.names
  #   all.cells:  the all cells of TPM data
  #   cells: the cells use to select 
  #   do.plot: whether plot
  if(!gene%in%rownames(tpm.data)){
    cat(sprintf('%s is not in TPM data\n',gene))
    do.plot=FALSE
    return(FALSE)
  }else{
  gene.mean.tpm<-tpm.data[gene,all.cells%in%cells]
  if(do.plot){
    plot.DiffusionMap(DMap.obj,dims = c(1,3,4),
                      col =as.numeric(gene.mean.tpm),pch=20,legend_main = gene)
    
  }else{
    return(gene.mean.tpm)
  }
  }
}

Gkn3<-DiffusionMap_Meantpm(gene = 'Gkn3',DMap.obj = all.dif,tpm.data = Regional_tpm)
Fabp1<-DiffusionMap_Meantpm(gene = 'Fabp1',DMap.obj = all.dif,tpm.data = Regional_tpm)
Bex1<-DiffusionMap_Meantpm(gene = 'Bex1',DMap.obj = all.dif,tpm.data = Regional_tpm)
Fabp6<-DiffusionMap_Meantpm(gene = 'Fabp6',DMap.obj = all.dif,tpm.data = Regional_tpm)
