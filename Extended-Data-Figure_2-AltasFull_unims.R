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

get.variable.genes_cvdiff <- function(ed,do.plot=T){
  
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv <- sqrt(vars)/means
  #minMeanForFit <- unname( quantile( means[ which( cv > min.cv ) ], .95 ) )
  #useForFit <- means >= minMeanForFit # & spikeins
  #info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
  fit <- glm.fit( cbind( a0 = 1, a1tilde = 1/means),cv) #[useForFit] ), cv[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv))}
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
  
  df <- ncol(ed) - 1
  # add confidence interval
  if(do.plot){
    lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
    lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
  }
  afit <- a1/means+a0  # fited value
  CVdiff<-cv-afit
  
  CVdiff_mean<-mean(CVdiff)
  CVdiff_std<-sqrt(var(CVdiff))
  delta<-CVdiff_mean+1.67*CVdiff_std
  
  r<-data.frame(rownames(ed),CVdiff)
  colnames(r)<-c('Genes','CVdiff')
  r<-r[r$CVdiff>delta,]
  return(r)
}

source('Fxns.R')




# clsuetr
atlas_full_tpm<-read.delim('./data/GSE92332_AtlasFullLength_TPM.txt.gz')
info(sprintf("Data dimensions: %s" , paste(dim(atlas_full_tpm), collapse = "x")))
v = get.variable.genes(atlas_full_tpm, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])
atlas_full_tpm<-data.frame(log2(1+atlas_full_tpm))



atlas_full_tpm_sample<-colnames(atlas_full_tpm)
cell_types<-unlist(lapply(atlas_full_tpm_sample,function(x)return(str_split(x,"_")[[1]][4])))


if(file.exists('./data/AtlasFull_pca_scores.txt')){
  pca = read.table("./data/AtlasFull_pca_scores.txt")
}else{
  pca = rpca(t(atlas_full_tpm[var.genes,]), center=T, scale=T, retx=T, k=100)$x   #use selected variables for PCA
  write.table(pca,'./data/AtlasFull_pca_scores.txt',quote = F)
}
#pca = read.delim(file="regions_pca_scores.txt")



#   Test for sinificant principle compnents
y = sig.pcs.perm(dat=t(atlas_full_tpm[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)  # 20 PC


if(file.exists('./data/AtlasFull_tsne.txt')){
  tsne.rot = read.table("./data/AtlasFull_tsne.txt")
}else{
  barnes_hut_tsne = Rtsne(pca[, 1:y$r], check_duplicates=T,pca=FALSE, #dont run PCA again
                          initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
  tsne.rot = barnes_hut_tsne$Y
  write.table(tsne.rot,'./data/AtlasFull_tsne.txt',quote = F)
}




dm = as.matrix(dist(pca[, 1:14]))
# build nearest neighbor graph
knn = build_knn_graph(dm, k = 80)
clustering = cluster_graph(knn)$partition
clustering<-merge_clusters(clustering,c(11,12))
clustering<-merge_clusters(clustering,c(1,11))



#   Figure a

colnames(tsne.rot)<-c('tSNE_1','tSNE_2')
tsne.rot<-as.data.frame(tsne.rot)
ggplot(tsne.rot,aes(x=tSNE_1,y=tSNE_2,color=cell_types))+geom_point()+scale_color_manual(values=brewer16)
ggsave('AtlasFull_TSNE.jpeg')

Genes_mean_tpm<-function(genes,tpm_data,tsne_data,title,fun=mean,doplot=TRUE){
  Log2TPM<-as.numeric(apply(tpm_data[genes,],2,fun))
  if(doplot){
    title_1<-paste(genes,collapse = ',')
    title_2<-paste(title,'(',title_1,')',sep='')
    ggplot(tsne_data, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Log2TPM))+theme(legend.title = element_text(size=8,color='blue',face='bold'),
                                                                            legend.position = 'right') +ggtitle(title_2)+
                                                                           scale_color_gradient2(low='lightblue',mid='green',high='red')
    }
  else{
    return(Log2TPM)
  }
}
# stem mark genes
stem_mark_genes<-c('Lgr5','Ascl2','Slc12a2','Axin2','Olfm4','Gkn3')
Genes_mean_tpm(stem_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Stem')



# Cell cycle
cell_cycle_mark_genes<-c('Mki67','Cdk4','Mcm5','Mcm6','Pcna')
Genes_mean_tpm(cell_cycle_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Cell cycle')


# Enterocyte
Enterocyte_mark_genes<-c('Alpi','Apoa1','Apoa4','Fabp1')
Genes_mean_tpm(Enterocyte_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Enterocyte')


# Globlet
Globlet_mark_genes<-c('Muc2','Clca1','Tff3','Agr2')
Genes_mean_tpm(Globlet_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Globet')

# Paneth
Paneth_mark_genes<-c('Lyz1','Defa17','Defa22','Defa24','Ang4')
Genes_mean_tpm(Paneth_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Paneth')


#  Enteroendocrine
Enteroendocrine_mark_genes<-c('Chga','Chgb','Tac1','Tph1','Neurog3')
Genes_mean_tpm(Enteroendocrine_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Enteroendocrine')



# Tuft
Tuft_mark_genes<-c('Dclk1','Trpm5','Gfi1b','Il25')
Genes_mean_tpm(Tuft_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Tuft')





# number of detected genes
Count_genes<-function(x){
  count<-0
  for(c in x){
    if(c!=0){
      count<-count+1
    }
  }
  return(count)
}

Genes_per_cell<-as.numeric(apply(atlas_full_tpm,2,Count_genes))
ggplot(x, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Genes_per_cell))+ggtitle('Genes/Cell')+
  theme(legend.title = element_text("Genes/Cell",size=8,color='blue',face='bold'),
                  legend.position = 'right') +scale_color_gradient2(low='lightblue',mid='green',high='red')




#      Heatmap

Create_plot_data<-function(genes,fun=scale,origin.data=atlas_full_tpm,var_genes=var.genes,cell_groups=cell_types,
                           if.use.var.genes=FALSE,cells=NULL){

  #origin.data: [genes,cells]
  
  if(if.use.var.genes){
    origin.data_1<-origin.data[var_genes,]
  }else{
    origin.data_1<-origin.data
  }
  
  origin.data_2<-as.data.frame(t(origin.data_1))
  
  if(!is.null(cells)){
    heatp_1<-origin.data_2[cell_groups%in%cells,]
  
    heatp_2<-heatp_1[,colnames(heatp_1)%in%genes]
    heatp_3<-as.data.frame(apply(heatp_2,2,fun))
    heatp_3[,'Groups']<-cell_groups[cell_groups%in%cells] 
  }
  if(is.null(cells)){
    heatp_2<-origin.data_2[,colnames(origin.data_2)%in%genes]
    heatp_3<-as.data.frame(apply(heatp_2,2,fun))
    heatp_3[,'Groups']<-cell_groups
  }

  heatp_4<-melt(heatp_3)
  heatp_4<-heatp_4[order(heatp_4$Groups),]
  
  return(heatp_4)
}


## Extended data Figure 2  
###  Figure b


figure_1c_genes<-read.table('./data/Figure_1c_genes.txt')
figure_1c_genes<-as.character(figure_1c_genes$V1)

figure_b<-Create_plot_data(figure_1c_genes,cells =c("Paneth","Endocrine","Goblet","Tuft","Enterocyte"),
                           if.use.var.genes = TRUE)



p<-ggplot(figure_b,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
  scale_fill_continuous(low='lightblue',high='red')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color="blue", size=8,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(color='blue',size=4),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))


# Figure c
Mptx2_mark_genes<-c('Mptx2')
Genes_mean_tpm(Mptx2_mark_genes,tpm_data = atlas_full_tpm,tsne_data =tsne.rot,title = 'Mptx2')


atlas_umis = read.delim("./data/GSE92332_atlas_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(atlas_umis), collapse = "x")))
atlas_tpm = data.frame(log2(1+tpm(atlas_umis)))

tsne.rot = read.table("./data/atlas_tsne.txt")
colnames(tsne.rot)<-c('tSNE_1','tSNE_2')
Genes_mean_tpm(Mptx2_mark_genes,tpm_data = atlas_tpm,tsne_data =tsne.rot,title = 'Mptx2')



### Figure d
GPCRs<-read.table('./data/GPCRS.txt')
GPCRs<-as.character(GPCRs$V1)


#atlas_tpm_GPCRs<-atlas_full_tpm   #[var.genes,]
#atlas_tpm_GPCRs<-as.data.frame(t(atlas_tpm_GPCRs))
#atlas_GPCRs<-atlas_tpm_GPCRs[,colnames(atlas_tpm_GPCRs)%in%GPCRs]
#atlas_GPCRs<-as.data.frame(apply(atlas_GPCRs,2,scale))


#atlas_GPCRs[,'Groups']<-cell_types
#atlas_GPCRs_melt<-melt(atlas_GPCRs)
#atlas_GPCRs_melt<-atlas_GPCRs_melt[order(atlas_GPCRs_melt$Groups),]
atlas_gcprs<-Create_plot_data(genes = GPCRs)

p<-ggplot(atlas_gcprs,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
  scale_fill_continuous(low='lightblue',high='red')+ggtitle('GPCRS')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color="blue", size=8,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(color='blue',size=6),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))




### Figure e
LRRS<-read.table('./data/LRRS.txt')
LRRS<-as.character(LRRS$V1)

atlas_lrrs<-Create_plot_data(genes = LRRS)
p<-ggplot(atlas_lrrs,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
  scale_fill_continuous(low='lightblue',high='red')+ggtitle('LRRS')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color="blue", size=8,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(color='blue',size=6),
                                                             panel.border = element_blank(),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             axis.line = element_line(colour = "black"))

### Figure f
###  f1
transition_factors<-read.table('./data/Transition-factors.txt')
transition_factors<-as.character(transition_factors$V1)

atlas_transitions_f1<-Create_plot_data(genes=transition_factors)
p<-ggplot(atlas_transitions_f1,aes(x=Groups,y=variable))+geom_tile(aes(fill=value))+
  scale_fill_continuous(low='lightblue',high='red')+ggtitle('TFs')+theme_bw()


p+labs(x='',y='')+scale_x_discrete(expand = c(0,0),position = 'bottom')+
  scale_y_discrete(expand = c(0,0),position = 'right')+theme(axis.text.x = element_text(face="bold", color="blue", size=8,hjust = 1,angle = 90),
                                                             axis.text.y = element_text(color='blue',size=4),
                                                             panel.grid.major.y = element_blank(),
                                                             panel.grid.minor.y = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line(colour = "black"))



