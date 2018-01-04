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


# data
LargeCellSort_UMIs<-read.delim('./data/GSE92332_LargeCellSort_UMIcounts.txt.gz')
info(sprintf("Data dimensions: %s" , paste(dim(LargeCellSort_UMIs), collapse = "x")))

###  whether need this step??
gene_expr<-as.numeric(apply(LargeCellSort_UMIs,1,sum))
LargeCellSort_UMIs<-LargeCellSort_UMIs[which(gene_expr>0),]
###  

v = get.variable.genes(LargeCellSort_UMIs, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])
LargeCellSort_tpm<-data.frame(log2(1+tpm(LargeCellSort_UMIs)))  # transform for tpm value



Cell_groups<-unlist(lapply(colnames(LargeCellSort_UMIs),function(x) return(str_split(x,'_')[[1]][4])))
batches<-unlist(lapply(colnames(LargeCellSort_UMIs),function(x) return(str_split(x,'_')[[1]][3])))
mices<-unlist(lapply(colnames(LargeCellSort_UMIs),function(x) return(str_split(x,'_')[[1]][2])))

table(batches)   

# check batch effect
batch_mean_tpm = group.means(counts = LargeCellSort_tpm, groups = batches)
x = batch_mean_tpm[, 1]
y = batch_mean_tpm[,2]
expr.cor = round(cor(x,y),2)
smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("Before batch correction, correlation between \ntwo illustrative batches is %s", expr.cor), xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")
#   result show that there is not batch effect



# LargeCellsort_tpm_norm = batch.normalise.comBat(counts = as.matrix(LargeCellSort_tpm), batch.groups = batches)
# batch_mean_tpm_norm = group.means(counts = LargeCellsort_tpm_norm, groups = batch.labels)
# x = batch_mean_tpm_norm[, 1]
# y = batch_mean_tpm_norm[,2]
# expr.cor = round(cor(x,y),2)
# smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("After batch correction, correlation between \ntwo illustrative batches is %s", expr.cor), xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")


# Dim Redection
#if(file.exists('LargeCellSort_pca_scores.txt')){
#  pca<-read.table('LargeCellSort_pca_scores.txt')
#}else{
#  pca = rpca(t(LargeCellSort_tpm[var.genes,]), center=T, scale=T, retx=T, k=100)$x
#  write.table(pca,'LargeCellSort_pca_scores.txt',quote = F)
#}
#y = sig.pcs.perm(dat=t(LargeCellSort_tpm[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)  

# if(file.exists('LargeCellSort_tsne.txt')){
#  tsne.rot = read.table("LargeCellSort_tsne.txt")
# }else{
#  barnes_hut_tsne = Rtsne(pca[, 1:13], check_duplicates=T,pca=FALSE, #dont run PCA again
#                          initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
#  tsne.rot = barnes_hut_tsne$Y
#  write.table(tsne.rot,'LargeCellSort_tsne.txt',quote = F)
# }


tsne.rot<-PCA_TSNE.scores(data.tpm=LargeCellSort_tpm,data.umis=LargeCellSort_UMIs,var_genes = var.genes,data_name = './data/LargeCellSort')
colnames(tsne.rot)<-c('tSNE_1','tSNE_2')

#   Figure a

###   It is terrible,the picture shows.And  why  ??? ,nedd to check
ggplot(data = tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=Cell_groups))+geom_point()+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())


#  Figure b

genes.names<-rownames(LargeCellSort_tpm)
cells<-c("Paneth.2","Paneth.1")
genes<-c('Defa20','Gm15308','Defa22','Defa21','Gm15315','Gm21002','Nupr1',
         'Gm15293','Pnliprp2','Itln1','Rnase1','AY761184','Gm15284','Defa23',
         'Clps','Defa17','Gm15299','Gm15292','Guca2b')




heatmap.1<-Create_plot_data(genes = genes,origin.data = LargeCellSort_tpm,
                            cell_groups = Cell_groups,cells = cells,do.plot=TRUE)



# library(scales)
# numColors <- length(levels(heatmap.1$Groups)) # How many colors you need
# getColors <- brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package)
# myPalette <- getColors(numColors)
# names(myPalette) <- levels(heatmap.1$Groups)

p<-ggplot(heatmap.1,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+scale_fill_gradient2(mid ='grey')
  #scale_fill_continuous(low='lightblue',high='red')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color=c("black","red"), size=10,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(size=8),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))



#   Figure c
Regional_UMIs<-load_data("./data/GSE92332_Regional_UMIcounts.txt.gz")
Regional_tpm<-data.frame(log2(1+tpm(Regional_UMIs)))
regional.genes<-rownames(Regional_UMIs)
region_groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][2])))

cell.groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][4])))



#genes.regions.tpm<-Regional_tpm[regional.genes%in%genes,]
heatmap.2<-Create_plot_data(genes = genes,origin.data = Regional_tpm[,cell.groups%in%'Paneth'],cell_groups = region_groups[cell.groups%in%'Paneth'])
p<-ggplot(heatmap.2,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
    scale_fill_continuous(low='lightblue',high='red')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", size=10,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(size=8),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))


#  Figure e
genes.e<-c('Eef1g','Gstm1','Rgcc','Wfdc10','Amica1','Ifitm3','Gkn3','Uba52','Acot1',
           'Nrtn','Eef1b2','Hmgcs2','Hspa8','Pdgfa','Btf3','2210407C18Rik','Clca3b','Sord',
           'Tfpi','Bex4','Bex1','Gas6','Bspry','Cyba')
heatmap.3<-Create_plot_data(genes = genes.e,origin.data = Regional_tpm[,cell.groups%in%'Stem'],cell_groups = region_groups[cell.groups%in%'Stem'])  # use Stem data

p<-ggplot(heatmap.3,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
  scale_fill_continuous(low='lightblue',high='red')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", size=10,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(size=8),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))

