library(KernSmooth)
library(NMF)
library(rsvd)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)
library(destiny)
library(stringr)
library(reshape2)
library(vioplot)

source('Fxns.R')

#   load data
Salmonellalnfect_UMIs<-load_data('/data/GSE92332_SalmonellaInfect_UMIcounts.txt.gz')
Salmonellalnfect_UMIs<-Salmonellalnfect_UMIs[which(unlist(apply(Salmonellalnfect_UMIs,1,sum))>0),]


#  caculate the log2(1+tpm)
Salmonellalnfect_tpm<-data.frame(log2(1+tpm(Salmonellalnfect_UMIs)))

# get variables
Sal.v = get.variable.genes(Salmonellalnfect_UMIs, min.cv2 = 100)
Sal.var.genes = as.character(rownames(Sal.v)[Sal.v$p.adj<0.05])  # select genes

#  message from sample
Salmonellalnfect.all.genes<-rownames(Salmonellalnfect_UMIs)
Salmonellalnfect.all.cells<-colnames(Salmonellalnfect_UMIs)


Salmonellalnfect.cells.group<-unlist(lapply(Salmonellalnfect.all.cells,function(x)return(str_split(x,'_')[[1]][4])))
Salmonellalnfect.condition<-unlist(lapply(Salmonellalnfect.all.cells,function(x)return(str_split(x,'_')[[1]][2])))
Salmonellalnfect.mices<-unlist(lapply(Salmonellalnfect.all.cells,function(x)return(str_split(x,'_')[[1]][3])))


#   
#  Dimensionality reduction
Salmonellalnfect.tsne.rot<-PCA_TSNE.scores(data.tpm = Salmonellalnfect_tpm,data.umis = Salmonellalnfect_UMIs,var_genes = Sal.var.genes,
                                   data_name = './data/Salmonellalnfect',sig.pcs = TRUE,is.var.genes = TRUE)

Salmonellalnfect.tsne.rot<-data.frame(Salmonellalnfect.tsne.rot)
colnames(Salmonellalnfect.tsne.rot)<-c('tSNE_1','tSNE_2')


#  Figure b
Sal.genes.b<-c('Reg3b','Reg3g','Clec2h','Fabp1','Anpep')
# for(gene in Sal.genes.b){
#   Genes_mean_tpm(genes = gene,tpm_data = Salmonellalnfect_tpm[,Salmonellalnfect.condition%in%'Control'],
#                tsne_data = Salmonellalnfect.tsne.rot[Salmonellalnfect.condition%in%'Control',],title ='Control')
# }
# 
# for(gene in Sal.genes.b){
#   Genes_mean_tpm(genes = gene,tpm_data = Salmonellalnfect_tpm[,Salmonellalnfect.condition%in%'Salm'],
#                  tsne_data = Salmonellalnfect.tsne.rot[Salmonellalnfect.condition%in%'Salm',],title ='Salm')
# }


# ggplot(data=tsne,aes(x=tSNE_1,y=tSNE_2,colour=Gene.Mp))+geom_point()+
#   facet_wrap(~ Condition, ncol = 2)+
#   theme(legend.title = element_text(size=10,color='blue',face='bold'),legend.position = 'right')+
#   scale_color_gradient2(low='lightblue',mid='green',high='red',name='Log2\nTPM+1')

### function for facet_wrap
Facet_wrap_fun<-function(gene,tpm.data,tsne.data,
                         condition=c('Control','Salm'),all.condition=Salmonellalnfect.condition){
  cat(sprintf('There ara %d conditions\n',length(condition)))
  tsne<-data.frame()
  for(i in 1:length(condition)){
    tsne<-rbind(tsne,tsne.data[all.condition%in%condition[i],])
  }
  cat(sprintf('Whether creat data accurate %d \n',sum(dim(tsne.data)[1]==dim(tsne)[1])))
  
  ###  create  gene expression TPM data
  gene.mp<-c()
  for(i in 1:length(condition)){
    gene.mp<-c(gene.mp,as.numeric(tpm.data[gene,all.condition%in%condition[i]]))
  }
  
  ### create Condition
  Condition<-c()
  for(i in 1:length(condition)){
    Condition<-c(Condition,rep(condition[i],sum(all.condition%in%condition[i])))
  }
  tsne$Gene.Mp<-gene.mp
  tsne$Condition<-Condition
  tsne$Gene<-rep(gene,dim(tsne)[1])
  
  return(tsne)
}

All.Facet.tsne<-data.frame()
for(gene in Sal.genes.b){
  All.Facet.tsne<-rbind(All.Facet.tsne,
                        Facet_wrap_fun(gene=gene,tpm.data = Salmonellalnfect_tpm,tsne.data = Salmonellalnfect.tsne.rot))
}


ggplot(data=All.Facet.tsne,aes(x=tSNE_1,y=tSNE_2,colour=Gene.Mp))+geom_point()+
  facet_wrap(~ Gene+Condition, nrow = 5,ncol = 2)+
  theme(legend.title = element_text(size=10,color='blue',face='bold'),legend.position = 'right')+
  scale_color_gradient2(low='lightblue',mid='green',high='red',name='Log2\nTPM+1')




#  Figure d
heatmap.d1.genes<-c("Asns",'Calm2','Dnase1','Frk','Gstp1','Nlrp6','Mthfd2','Pfkfb3',
                    'Tgm2','Tnfsf10')


# Prepare data
###   d1
Control.tpm<-Salmonellalnfect_tpm[heatmap.d1.genes,Salmonellalnfect.condition%in%'Control']
Salm.tpm<-Salmonellalnfect_tpm[heatmap.d1.genes,Salmonellalnfect.condition%in%'Salm']

Control.Salm.tpm<-cbind(Control.tpm,Salm.tpm)
annCol.d1<-c(rep('Control',dim(Control.tpm)[2]),rep('Salm',dim(Salm.tpm)[2]))
NMF::aheatmap(Control.Salm.tpm,Rowv = NA,Colv = NA,
              annCol = annCol.d1,
              scale = 'row')

### d2

#### load  SalmHelm_UMIs data

heatmap.d2.genes<-c('Calm2','Gpr160','Frk','Il4ra','Nars','Mthfd2','Pfkfb3','Tgm2','Tnfsf10')
SalmHelm_UMIs<-load_data(".data/GSE92332_SalmHelm_UMIcounts.txt.gz")
SalmHelm_UMIs<-SalmHelm_UMIs[which(unlist(apply(SalmHelm_UMIs,1,sum))>0),]
v = get.variable.genes(SalmHelm_UMIs, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])  # select genes

SalmHelm_tpm<-data.frame(log2(1+tpm(SalmHelm_UMIs)))
SalmHelm.all.cells<-colnames(SalmHelm_UMIs)
SalmHelm.all.genes<-rownames(SalmHelm_UMIs)

SalmHelm.condition<-unlist(lapply(SalmHelm.all.cells,function(x)return(str_split(x,'_')[[1]][3])))    # cell by conditon
SalmHelm.cell.groups<-unlist(lapply(SalmHelm.all.cells,function(x)return(str_split(x,'_')[[1]][4])))  # cell groups
SalmHelm.batches<-unlist(lapply(SalmHelm.all.cells,function(x)return(str_split(x,'_')[[1]][1])))  #  batches


#### remove batch effect
SalmHelm_tpm_norm = batch.normalise.comBat(counts = as.matrix(SalmHelm_tpm), batch.groups = SalmHelm.batches)

Control.SalmHelm.tpm<-SalmHelm_tpm_norm[heatmap.d2.genes,SalmHelm.condition%in%'Control']
HpolyD10.SalmHelm.tpm<-SalmHelm_tpm_norm[heatmap.d2.genes,SalmHelm.condition%in%'Hpoly.Day10']
HpolyD3.SalmHelm.tpm<-SalmHelm_tpm_norm[heatmap.d2.genes,SalmHelm.condition%in%'Hpoly.Day3']
Salmonella.SalmHelm.tpm<-SalmHelm_tpm_norm[heatmap.d2.genes,SalmHelm.condition%in%'Salmonella']

CHHS.tpm<-cbind(Control.SalmHelm.tpm,HpolyD10.SalmHelm.tpm,HpolyD3.SalmHelm.tpm,Salmonella.SalmHelm.tpm)
annCol.d2<-c(rep('Control',dim(Control.SalmHelm.tpm)[2]),
            rep("Hpoly.Day10",dim(HpolyD10.SalmHelm.tpm)[2]),
             rep("Hpoly.Day3",dim(HpolyD3.SalmHelm.tpm)[2]),
             rep("Salmonella",dim(Salmonella.SalmHelm.tpm)[2]))
NMF::aheatmap(CHHS.tpm,Rowv = NA,Colv = NA,
              annCol = annCol.d2,
              scale = 'row')



#    Figure e

Cells<-c("Enterocyte.dist","Enterocyte.prox","Tuft","Goblet","Enteroendocrine","Paneth")
Sal.Saa1.tpm<-data.frame(t(Salmonellalnfect_tpm[c('Saa1'),Salmonellalnfect.cells.group%in%Cells]))
Sal.Saa1.tpm$Cell<-Salmonellalnfect.cells.group[Salmonellalnfect.cells.group%in%Cells]
Sal.Saa1.tpm$Condition<-Salmonellalnfect.condition[Salmonellalnfect.cells.group%in%Cells]
colnames(Sal.Saa1.tpm)<-c('Value','Cell','Condition')
Sal.Saa1.tpm$Gene<-rep('Saa1',dim(Sal.Saa1.tpm)[1])

Sal.Saa2.tpm<-data.frame(Gene=t(Salmonellalnfect_tpm[c('Saa2'),Salmonellalnfect.cells.group%in%Cells]))
Sal.Saa2.tpm$Cell<-Salmonellalnfect.cells.group[Salmonellalnfect.cells.group%in%Cells]
Sal.Saa2.tpm$Condition<-Salmonellalnfect.condition[Salmonellalnfect.cells.group%in%Cells]
colnames(Sal.Saa2.tpm)<-c('Value','Cell','Condition')
Sal.Saa2.tpm$Gene<-rep('Saa2',dim(Sal.Saa2.tpm)[1])

Sal.Saa.tpm<-rbind(Sal.Saa1.tpm,Sal.Saa2.tpm)

ggplot(data=Sal.Saa.tpm,aes(x=Cell,y=Value))+
  geom_violin(aes(fill=Value,colour=Condition))+
  geom_jitter(aes(group=Condition),shape=16, position=position_jitter())+
  facet_wrap(~ Gene, nrow = 2)+
  xlab(NULL)+ylab('Expression\n(Log2(TPM+1)')


# ggplot(data=Sal.Saa.tpm,aes(x=Cell,y=Saa2))+
#   geom_violin(aes(fill=Saa2,colour=Condition))+
#   geom_jitter(aes(group=Condition),shape=16, position=position_jitter())+
#   xlab(NULL)+ylab('Saa2 Expression\n(Log2(TPM+1)')


# Figure  f

Paneth.f1.tpm<-data.frame(t(Salmonellalnfect_tpm[c('Ang4','Defa23','Defa24'),Salmonellalnfect.cells.group%in%'Paneth']))
Paneth.f1.tpm$Condition<-Salmonellalnfect.condition[Salmonellalnfect.cells.group%in%'Paneth']
Paneth.f1.tpm.melt<-melt(Paneth.f1.tpm)

ggplot(data=Paneth.f1.tpm.melt,aes(x=variable,y=value))+
  geom_violin(aes(fill=value,colour=Condition))+
  scale_x_discrete()+
  geom_boxplot(aes(colour=Condition),width=0.1)+
  geom_jitter(aes(group=Condition),shape=16, position=position_jitter())+
  xlab(NULL)+ylab('Expression Level\nLog2(TPM+1)')




#  Figure g
Paneth.Condition.mices.tab<-as.matrix(table(Salmonellalnfect.condition[Salmonellalnfect.cells.group%in%'Paneth'],
                            Salmonellalnfect.mices[Salmonellalnfect.cells.group%in%'Paneth']))
attr(Paneth.Condition.mices.tab,'class')<-NULL
Paneth.Condition.mices.tab<-data.frame(Paneth.Condition.mices.tab)/as.vector(table(Salmonellalnfect.condition))
# Paneth.Condition.mices.tab$Condition<-c('Control','Salm')
# Paneth.Condition.mices.tab.melt<-melt(Paneth.Condition.mices.tab)


All.tab<-as.matrix(table(Salmonellalnfect.cells.group,Salmonellalnfect.condition))
attr(All.tab,'class')<-NULL
All.tab<-data.frame(All.tab)
All.tab<-All.tab/unlist(apply(All.tab,2,sum))
All.Paneth.tab<-All.tab['Paneth',]

barplot(as.matrix(All.Paneth.tab),beside = TRUE,width = 0.1,ylim = c(0.0,0.03),col = c('grey','pink'))


Mice.prop<-data.frame(Mice.Prop=unlist(apply(Paneth.Condition.mices.tab,2,sum)))
Mice.prop$Condition<-c(rep('Control',4),rep('Salm',4))
Mice.prop$Total<-c(rep(All.Paneth.tab[,'Control'],4),rep(All.Paneth.tab[,'Salm'],4))

p<-ggplot(data = Mice.prop,aes(x=Condition))+
  geom_bar(aes(y=Total,fill=Condition),stat="identity",position = 'dodge')

p<-p+geom_point(aes(y=Mice.Prop),stat="identity",position=position_dodge(width = 0.9),alpha=.8,
                size=3)+ylim(0.0,0.030)+ylab('Fraction of Paneth cells')


