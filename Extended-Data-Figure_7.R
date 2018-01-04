library(NMF)
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

source('Fxns.R')
FullAtlas_UMIs<-load_data('./data/GSE92332_AtlasFullLength_TPM.txt.gz')
FullAtlas_UMIs<-FullAtlas_UMIs[which(unlist(apply(FullAtlas_UMIs,1,sum))>0),]

FullAtlas_tpm<-data.frame(log2(1+tpm(FullAtlas_UMIs)))
FullAtlas.v<-get.variable.genes(FullAtlas_tpm,min.cv2 = 100)
FullAtlas.var.genes<-as.character(rownames(FullAtlas.v)[FullAtlas.v$p.adj<0.05])
all.genes<-rownames(FullAtlas_tpm)
all.cells<-colnames(FullAtlas_tpm)

cell.groups<-unlist(lapply(all.cells,function(x)return(str_split(x,'_')[[1]][4])))
Tuft.cells<-all.cells[cell.groups%in%'Tuft']

Tuft.tpm<-FullAtlas_tpm[FullAtlas.var.genes,Tuft.cells]
Tuft.umis<-FullAtlas_UMIs[FullAtlas.var.genes,Tuft.cells]


#   Figure a
Tuft.tsne.rot<-PCA_TSNE.scores(data.tpm=Tuft.tpm,data.umis=Tuft.umis,
                                 data_name = './data/FullAtlas.Tuft',is.var.genes = FALSE,var_genes = NULL,sig.pcs = FALSE)


colnames(Tuft.tsne.rot)<-c('tSNE_1','tSNE_2')
Tuft.pca<-read.table('./data/FullAtlas.Tuft_pca_scores.txt')
dm<-as.matrix(dist(Tuft.pca[,1:15]))
# build nearest neighbor graph
knn = build_knn_graph(dm, k = 40)
clustering = cluster_graph(knn)$partition   # only 1 cluster,why ???


## Try kmeans method for subcluster,may be not right
Tuft.kmeans<-kmeans(Tuft.tsne.rot,2)
Tuft.clustering<-as.character(Tuft.kmeans$cluster)
Tuft.clustering<-paste('Tuft-',Tuft.clustering,sep='')

ggplot(data = Tuft.tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=Tuft.clustering))+geom_point()+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())



###    Figure b


# library(DESeq2)
# condition<-as.factor(Tuft.clustering)
# countData<-as.matrix(Tuft.umis)
# dds<-DESeqDataSetFromMatrix(countData =countData,DataFrame(condition),~condition)

Tuft_1_marker_genes<-as.character(read.table('./data/Tuft_1_marker_genes.txt',header = FALSE)$V1)
Tuft_2_marker_genes<-as.character(read.table('./data/Tuft_2_marker_genes.txt',header = FALSE)$V1)
Tuft_12_marker.genes<-c(Tuft_1_marker_genes[1:25],Tuft_2_marker_genes[1:25])

Tuft.1.tpm<-Tuft.tpm[Tuft_12_marker.genes,which(Tuft.clustering=='Tuft-1')]
Tuft.2.tpm<-Tuft.tpm[Tuft_12_marker.genes,which(Tuft.clustering=='Tuft-2')]

Tuft.marker.12.genes.tpm<-cbind(Tuft.1.tpm,Tuft.2.tpm)
annCol<-c(rep('Tuft-1',dim(Tuft.1.tpm)[2]),rep('Tuft-2',dim(Tuft.2.tpm)[2]))
annRow<-rep(c('Tuft-1','Tuft-2'),each=25)
# NMF::aheatmap(Tuft.marker.12.genes.tpm,Rowv = NA,Colv = NA,
#               annCol = Tuft.clustering,annRow = annRow,
#               scale = 'row')#,filename = 'Heatmap.Tuft.jpeg')
NMF::aheatmap(Tuft.marker.12.genes.tpm,Rowv = NA,Colv = NA,
              annCol = annCol,annRow = annRow,
              scale = 'row')#,filename = 'Heatmap.Tuft.jpeg'

##  Figure d

library(bayesboot)
 

S8.table.name<-c('Gene','p.max','fdr.max','Pfisher','FDR.fisher',
                 'FC.min','FC.lower','FC.mean','group.me','other.me','Fac.non.me')

genes<-c('Ninj1','Nradd','Nrep','Plekhg5','Lyn','Rhog','Il17rb','Irf7','Rac2')

Tuft1_FC<-read.csv('Tuft1_FC.csv',header = FALSE,sep = '\t')
Tuft2_FC<-read.csv('Tuft2_FC.csv',header = FALSE,sep='\t')
colnames(Tuft1_FC)<-S8.table.name
colnames(Tuft2_FC)<-S8.table.name

Tuft_FC<-rbind(Tuft1_FC,Tuft2_FC)

# 
# Tuft.FC.genes<-as.character(Tuft_FC$Gene)
# for(g in genes){
#   if(!g%in%Tuft.FC.genes){
#     cat(sprintf('%s\n',g))
#   }
# }
rownames(Tuft_FC)<-Tuft_FC$Gene

Tuft.Genes.FC<-Tuft_FC[genes,c('FC.min','FC.lower','FC.mean')]

Caculate_SD<-function(x,N=10){
  zscore=qnorm(0.95)
  sd<-(x[3]-x[2])/zscore
  return(rnorm(n=N,mean = x[3],sd=sd))
}

Tuft.Genes.NormFc<-t(as.data.frame(apply(Tuft.Genes.FC,1,Caculate_SD)))
Tuft.Genes.bayesboot.FC<-t(data.frame(apply(Tuft.Genes.NormFc,1,bayesboot,mean,R=100)))
rownames(Tuft.Genes.bayesboot.FC)<-genes
Tuft.Genes.bayesboot.FC<-as.data.frame(Tuft.Genes.bayesboot.FC)
cols<-rep(c('red','blue'),time=c(4,5))
boxplot(log2(t(Tuft.Genes.bayesboot.FC)),outline=FALSE,horizontal=TRUE,col=cols)
legend('bottomright',legend = c('Inflammation-related','Neuron-related'),fill= c('blue','red'))


##  Figure e
Il33.Enterocyte<-FullAtlas_tpm['Il33',cell.groups%in%'Enterocyte']
Il33.Tuft<-data.frame(Il33=rep(0,dim(Tuft.tpm)[2]))
tpm.e<-rbind(t(Il33.Enterocyte),Il33.Tuft)
tpm.e$Cell<-c(rep('Enterocyte',dim(Il33.Enterocyte)[2]),Tuft.clustering)

ggplot(data=tpm.e,aes(x=Cell,y=Il33,fill=Cell))+geom_violin(aes(colour=Cell))+
  geom_jitter(shape=16, position=position_jitter())+xlab(NULL)+ylab('Il33 Expression\nLog2(TPM+2)')

## Figure i

### left top
Fig.I.tpm<-FullAtlas_tpm[,cell.groups%in%c('Enterocyte','Paneth','Stem','Tuft')]
Fig.I.umis<-FullAtlas_UMIs[,cell.groups%in%c('Enterocyte','Paneth','Stem','Tuft')]
Fig.I.tsne.rot<-PCA_TSNE.scores(data.tpm=Fig.I.tpm,data.umis=Fig.I.umis,
                               data_name = 'FullAtlas.Fig.I',is.var.genes = TRUE,var_genes =FullAtlas.var.genes,sig.pcs = FALSE)

Fig.I.tsne.rot<-data.frame(Fig.I.tsne.rot)
colnames(Fig.I.tsne.rot)<-c('tSNE_1','tSNE_2')
Fig.I.cells<-cell.groups[cell.groups%in%c('Enterocyte','Paneth','Stem','Tuft')]
ggplot(data = Fig.I.tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=Fig.I.cells))+geom_point()+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())

### Figure i top left

# Fig.I.tpm.tosort<-unlist(apply(Fig.I.tpm,2,mean))
# Fig.I.sort.cell<-colnames(Fig.I.tpm)[order(Fig.I.tpm.tosort,decreasing = T)]
# 
# Fig.I.sort.cell.tpm<-Fig.I.tpm[,Fig.I.sort.cell[1:332]]
# Fig.I.sort.cell.umis<-Fig.I.umis[,Fig.I.sort.cell[1:332]]
# Fig.I.sort.cells<-Fig.I.cells[order(Fig.I.tpm.tosort,decreasing = T)][1:332]
# 
# Fig.I.tsne.rot.sort<-PCA_TSNE.scores(data.tpm=Fig.I.sort.cell.tpm,data.umis=Fig.I.sort.cell.umis,
#                                 data_name = 'FullAtlas.Fig.I.sort.cell',is.var.genes = TRUE,var_genes =FullAtlas.var.genes,sig.pcs = FALSE)
# 
# Fig.I.tsne.rot.sort<-data.frame(Fig.I.tsne.rot.sort)
# colnames(Fig.I.tsne.rot.sort)<-c('tSNE_1','tSNE_2')
# ggplot(data = Fig.I.tsne.rot.sort,aes(x=tSNE_1,y=tSNE_2,color=Fig.I.sort.cells))+geom_point()+
#   scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())

###   Figure i top right
Dclk1<-Genes_mean_tpm(genes = 'Dclk1',tpm_data = Fig.I.tpm,tsne_data = Fig.I.tsne.rot,title =NULL,doplot = T)


###  Figure i bottom left,right
Tuft1_Meanexpr<-unlist(apply(FullAtlas_tpm[Tuft1_FC$Gene,cell.groups%in%c('Enterocyte','Paneth','Stem','Tuft')],2,mean))
Tuft2_Meanexpr<-unlist(apply(FullAtlas_tpm[Tuft2_FC$Gene,cell.groups%in%c('Enterocyte','Paneth','Stem','Tuft')],2,mean))


ggplot(Fig.I.tsne.rot, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Tuft1_Meanexpr))+theme(legend.title = element_text(size=8,color='blue',face='bold'),
                                                                                legend.position = 'right') +ggtitle('Tuft-1 score')+
                                                                                scale_color_gradient2(low='lightblue',mid='green',high='red',name='Log2\nTPM+1')

ggplot(Fig.I.tsne.rot, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Tuft2_Meanexpr))+theme(legend.title = element_text(size=8,color='blue',face='bold'),
                                                                                            legend.position = 'right') +ggtitle('Tuft-2 score')+
                                                                                scale_color_gradient2(low='lightblue',mid='green',high='red',name='Log2\nTPM+1')

