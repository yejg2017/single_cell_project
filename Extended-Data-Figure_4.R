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

source('Fxns.R')
#  load data
atlas_umis<-load_data(data_name = "GSE92332_atlas_UMIcounts.txt.gz")
atlas_umis<-atlas_umis[which(unlist(apply(atlas_umis,1,sum))>0),]
v = get.variable.genes(atlas_umis, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])  # select genes

get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
batch.labels = factor(unlist(lapply(colnames(atlas_umis), get_field, 1,"_")))
table(batch.labels)
atlas_tpm = data.frame(log2(1+tpm(atlas_umis)))

#  based on earlier analysis,wo knowed that this data has batch effect.
atlas_tpm_norm = batch.normalise.comBat(counts = as.matrix(atlas_tpm), batch.groups = batch.labels)

sample.names<-colnames(atlas_tpm_norm)
cell.types<-unlist(lapply(sample.names,function(x)return(str_split(x,'_')[[1]][3])))

###
diffusion_matrix<-cbind(data.frame(Cell=sample.names),t(atlas_tpm_norm[var.genes,])) 
rownames(diffusion_matrix)<-1:dim(diffusion_matrix)[1]


cells.1<-c('Stem','Enterocyte.Progenitor',"Enterocyte.Progenitor.Early","Enterocyte.Progenitor.Late" ,
           "Enterocyte.Immature.Distal","Enterocyte.Immature.Proximal")
ct.1 <- as.ExpressionSet(diffusion_matrix[cell.types%in%cells.1,])
dif.1<-DiffusionMap(ct.1,verbose = T,vars = NULL)  #   奇怪,这次样本没有减少 ???

# save(dif,file='Atlas_UMIs_DiffusionMap.RData')
DC<-as.data.frame(eigenvectors(dif.1))

## Figure a
cell.types.2<-unlist(lapply(cell.types[cell.types%in%cells.1],function(x){
  if(str_detect(x,'\\.')){
    return(str_split(x,'\\.')[[1]][2])
  }else{
    return(str_split(x,'\\.')[[1]][1])
  }
}))


ggplot(DC,aes(x=DC1,y=DC3))+geom_point(aes(color=cell.types[cell.types%in%cells.1]))+
 scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())+ggtitle('Enterocyte maturation')
# annotate(c('Immature','Progenitor','Stem'),x=c(-0.050,0.000,0.025),y=c(-0.025,0.050,-0.025))


plot(x=DC$DC1,y=DC$DC3,col=as.integer(as.factor(cell.types[cell.types%in%cells.1])),pch=20,cex=2,xlab='DC1',ylab='DC3',
     main='Enterocyte maturation')
text(x=c(-0.050,0.000,0.025),y=c(-0.025,0.050,-0.025),labels = c('Immature','Progenitor','Stem'),cex=2)
# legend(x=-0.06,y=0.04,legend =unique(cell.types[cell.types%in%cells.1]),
#        col = unique(as.integer(as.factor(cell.types[cell.types%in%cells.1]))))

gene.names<-rownames(atlas_tpm_norm)

## Figure c
plot_tpm<-function(gene,all.genes,cells,celltype,tpm.data,title=NULL,DC.data=DC,DC.F=c(1,2)){
  #  gene: the gene to caculate mean TPM expression
  #  all.genes: all genes of data
  #  cells: the cells of sample to select
  #  celltype : all the sample cells type
  #  tpm.data: the TPM data to use
  #  DC.data: data from DiffusionMap function
  #  DC.F: the component to plot  of DC.data
  #stopifnot(gene%in%rownames(tpm.data))
  gene.tpm<-tpm.data[all.genes%in%gene,celltype%in%cells]
  Logtpm<-as.numeric((apply(gene.tpm,2,mean)))
  xlabel<-paste('DC',DC.F[1],sep='-')
  ylabel<-paste('DC',DC.F[2],sep='-')
  print(ggplot(DC.data, aes(x=DC[,DC.F[1]], y=DC[,DC.F[2]]))+geom_point(aes(color=Logtpm))+theme(legend.title = element_text(size=8,color='blue',face='bold'),
                                                                                  legend.position = 'right') +ggtitle(title)+theme_bw()+labs(x=xlabel,y=ylabel)+
                                                                                  scale_color_gradient2(low='green',mid='blue',high='red',name='Log2\nTPM+1'))
  
  return(Logtpm)
}

Sox4<-plot_tpm(gene = 'Sox4',all.genes = gene.names,cells=cells.1,
           celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Sox4(stem/TA)',DC.F = c(1,3))

Foxm1<-plot_tpm(gene = 'Foxm1',all.genes = gene.names,cells=cells.1,
                celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Foxm1(progeniters)',DC.F = c(1,3))

Mxd3<-plot_tpm(gene = 'Mxd3',all.genes = gene.names,cells=cells.1,
               celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Mxd3(progeniters)',DC.F = c(1,3))

Batf2<-plot_tpm(gene = 'Batf2',all.genes = gene.names,cells=cells.1,
                celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Batf2(Immature Enterocyte)',DC.F = c(1,3))


##  Figure b
cells.2<-c("Enterocyte.Immature.Distal","Enterocyte.Immature.Proximal", "Enterocyte.Progenitor.Early",
           "Enterocyte.Progenitor.Late",'Stem',"TA","TA.G1","TA.G2")


ct.2<- as.ExpressionSet(diffusion_matrix[cell.types%in%cells.2,])
dif.2<-DiffusionMap(ct.2,verbose = T,vars = NULL)  # not reduce samples  # had saved:Atlas_UMIs_DiffusionMap_2.RData

DC<-as.data.frame(eigenvectors(dif.2))

ggplot(DC,aes(x=DC1,y=DC2))+geom_point(aes(color=cell.types[cell.types%in%cells.2]))+
  scale_color_manual(values=brewer16)+scale_fill_discrete()+theme(legend.title=element_blank())+ggtitle('Enterocyte maturation')



##  Figure d
Creb3l3<-plot_tpm(gene = 'Creb3l3',all.genes = gene.names,cells=cells.2,
               celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Creb3l3(proximal)')


Gata4<-plot_tpm(gene = 'Gata4',all.genes = gene.names,cells=cells.2,
                  celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Gata4(proximal)')


Nr1i3<-plot_tpm(gene = 'Nr1i3',all.genes = gene.names,cells=cells.2,
                  celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Nr1i3(proximal)')


Osr2<-plot_tpm(gene = 'Osr2',all.genes = gene.names,cells=cells.2,
               celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Osr2(distal)')

Jund<-plot_tpm(gene = 'Jund',all.genes = gene.names,cells=cells.2,
         celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Jund(distal)')

Nr1h4<-plot_tpm(gene = 'Nr1h4',all.genes = gene.names,cells=cells.2,
                celltype = cell.types,tpm.data = atlas_tpm_norm,title = 'Nr1h4(distal)')



## Figure f
Regional_UMIs<-load_data("GSE92332_Regional_UMIcounts.txt.gz")
Regional_UMIs<-Regional_UMIs[which(unlist(apply(Regional_UMIs,1,sum))>0),]
Regional_tpm<-data.frame(log2(1+tpm(Regional_UMIs)))
#region.var.genes<-get.variable.genes_cvdiff(Regional_tpm)

regional.v = get.variable.genes(Regional_UMIs, min.cv2 = 100)
regional.var.genes = as.character(rownames(regional.v)[regional.v$p.adj<0.05])
regional.genes<-rownames(Regional_tpm)
region_groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][2])))

cell.groups<-unlist(lapply(colnames(Regional_tpm),function(x)return(str_split(x,'_')[[1]][4])))
region.sample.names<-colnames(Regional_tpm)

#  only use Region Stem  Cells 
region.diffusion.matrix<-cbind(data.frame(Cell=region.sample.names[cell.groups%in%'Stem']),t(Regional_tpm[regional.var.genes,cell.groups%in%'Stem']))
region.ct<-as.ExpressionSet(region.diffusion.matrix)
region.dif<-DiffusionMap(region.ct)
#save(region.dif,file = 'Regional_UMIs_DiffusionMap_3.RData')

DC<-as.data.frame(eigenvectors(region.dif))
Lgr5<-plot_tpm(gene = 'Lgr5',all.genes = regional.genes,cells=c('Stem'),
                celltype =cell.groups ,tpm.data = Regional_tpm,DC.data = DC,DC.F = c(1,2),title = 'Lgr5(Stem)')


Gkn3<-plot_tpm(gene = 'Gkn3',all.genes = regional.genes,cells=c('Stem'),
         celltype =cell.groups ,tpm.data = Regional_tpm,DC.data = DC,DC.F = c(1,2),title = 'Gkn3(Stem)')

Bex1<-plot_tpm(gene = 'Bex1',all.genes = regional.genes,cells=c('Stem'),
               celltype =cell.groups ,tpm.data = Regional_tpm,DC.data = DC,DC.F = c(1,2),title = 'Bex1(Stem)')





