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

gene.names<-rownames(EEC_tpm)
sample.names<-colnames(EEC_tpm)

# Extract messages from samples
region.names<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][3])))
mices<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][2])))
cell.types<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][4])))


# do not known whether it has batch effect.  Ignore here.   If need to be checked,and will do it
tsne.rot<-PCA_TSNE.scores(data.tpm = EEC_tpm,data.umis=EEC_UMIs,var_genes = var.genes,data_name = 'EEC',is.var.genes = TRUE)
tsne.rot<-as.data.frame(tsne.rot)
colnames(tsne.rot)<-c('tSNE_1','tSNE_2')

## Figure a
ggplot(data = tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=region.names))+geom_point()+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())+ggtitle('Regions')

## Figure b  
## do not how to plot???


## Figure c
###   1.
ggplot(data = tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=cell.types))+geom_point()+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())+ggtitle('Cells Type')


###   2.

genes.p<-c('Neurog3','Sct','Tac1','Sst','Cck','Gcg','Ghrl','Gip','Nts','Reg4','Pyy')
for(gene in genes.p){
  if(!gene%in%gene.names){
    cat(sprintf('%s is not exist',gene))
  }
}


for(gene in genes.p){
  if(!gene%in%gene.names){
    next
  }else{
  print(Genes_mean_tpm(genes = gene,tpm_data = EEC_tpm,tsne_data = tsne.rot,title = NULL,doplot = TRUE))
  }
}


##  Figure d

#  (left)
heatmap.1<-Create_plot_data(genes = genes.p,origin.data = EEC_tpm,var_genes = var.genes,if.use.var.genes = TRUE,
                            cell_groups = cell.types)

p<-ggplot(heatmap.1,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+scale_fill_gradient2(mid ='grey',high='red')
#scale_fill_continuous(low='lightblue',high='red')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color=c("red"), size=10,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(size=8),
                                                             panel.grid.major.y  = element_blank(),
                                                             panel.grid.minor.y = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))



EEC_UMIs_1<-as.data.frame(t(EEC_UMIs[gene.names%in%genes.p,]))
Count<-lapply(unique(cell.types),function(x){
  x.1<-EEC_UMIs_1[cell.types%in%x,]
  x.2<-unlist(apply(x.1,2,sum))
  return(data.frame(x=round(x.2/sum(x.2),3)))
})
EEC_UMIs_prop<-as.data.frame(Count)*100
colnames(EEC_UMIs_prop)<-unique(cell.types)
#EEC_UMIs_prop<-data.frame(t(EEC_UMIs_prop)*100)


aheatmap(EEC_UMIs_prop, color = cubehelix1.16, border_color = list("cell"="white"), txt=EEC_UMIs_prop, Colv = NA, Rowv = NA )



