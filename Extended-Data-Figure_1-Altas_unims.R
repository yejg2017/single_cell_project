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

atlas_umis = read.delim("./data/GSE92332_atlas_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(atlas_umis), collapse = "x")))
v = get.variable.genes(atlas_umis, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])




get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
batch.labels = factor(unlist(lapply(colnames(atlas_umis), get_field, 1,"_")))
table(batch.labels)
atlas_tpm = data.frame(log2(1+tpm(atlas_umis)))

#  Extended Data  Figure 1 
##  Figure e
cell_groups<-unlist(lapply(colnames(atlas_umis),function(x){
  return(str_split(x,'_')[[1]][3])
}))

cell_types<-unique(cell_groups)

Stem_cells<-atlas_tpm[var.genes,cell_groups=='Stem']
stem_pearson<-cor(Stem_cells,method = 'pearson')

stem_pearson_matrix<-as.matrix(stem_pearson)
stem_pearson_vec<-stem_pearson_matrix[upper.tri(stem_pearson_matrix)]


get_cell_tpm<-function(cell){
  cell_matrix<-as.matrix(atlas_tpm[var.genes,cell_groups==cell])
  cell_pearson<-cor(cell_matrix,method = 'pearson')
  cell_pearson_vec<-cell_pearson[upper.tri(cell_pearson)]
  return(cell_pearson_vec)
}

cell_pearson_list<-lapply(cell_types,get_cell_tpm)
jpeg(file='./data/Extend_data_1/Figure_e.jpeg')
boxplot(cell_pearson_list,names=cell_types,outline=F,horizontal=T,las=2,
        xlab=c('Pearson correlation between mice,n=6'),ylab.cex=0.3)
dev.off()



## Figure c
colors = c('#483D8B','#00FFFF','#EEE8AA','#CD5C5C','#CD853F','#B22222','#CDC9C9','#7CFC00',
           '#668B8B','#008B45','#FF6A6A','#8B4726','#FF3030','#8B0A50','#4F4F4F')

cell_batch<-as.matrix(table(cell_groups,batch.labels))
cell_names<-rownames(cell_batch)
cell_tables<-table(cell_groups)

jpeg(file='./data/Extend_data_1/Figure_b.jpeg')
par(mfrow=c(5,4))
for(i in 1:length(cell_types)){
  #pie(cell_batch[i,],labels = names(cell_batch[i,]),main=str_c(cell_names[i],'(n=',cell_tables[i],')'))
  if(i==length(cell_types)){
      pie(cell_batch[i,],labels=colnames(cell_batch[i,]),col=colors,main=str_c(cell_names[i],'(n=',cell_tables[i],')'))
  }else{
     pie(cell_batch[i,],labels=NA,col=colors,main=str_c(cell_names[i],'(n=',cell_tables[i],')'))
  }
}
dev.off()
par(mfrow=c(1,1))
#legend('bottomright',legend =colnames(cell_batch),fill = colors,horiz = TRUE)  # not good ,will try more


###  Figure d

XcellTypes<-c("Endocrine","Enterocyte.Mature.Distal","Enterocyte.Mature.Proximal",
               "Goblet","Paneth","Stem",  "Tuft" )


Xcell<-atlas_tpm[,colnames(atlas_tpm)[cell_groups%in%XcellTypes]]
Xcell_types<-unlist(lapply(colnames(Xcell),function(x){
  return(str_split(x,'_')[[1]][3])
}))

Xcell_batches<-unlist(lapply(colnames(Xcell),function(x){
  return(str_split(x,'_')[[1]][1])
}))

Xcell_table<-as.matrix(table(Xcell_types,Xcell_batches))
Xcell_table_fac<-Xcell_table/apply(Xcell_table,1,sum)
Xcell_table_fac<-as.data.frame.matrix(Xcell_table_fac)

jpeg(file='./data/Extend_data_1/Figure_d.jpeg')
beeswarm(as.data.frame(t(Xcell_table_fac)),las=2,col='blue',pch=20,ylab=c('Factions of cells'))  #,col = as.numeric(as.factor(rownames(Xcell_table_fac))))
bxplot(as.data.frame(t(Xcell_table_fac)),probs=0.5,add=T)
dev.off()
#Boxplot(t(Xcell_table_fac),id.method='none',outline=F,notch=F,las=2)



### Figure g


## take mean tpm across batches to show batch effect
batch_mean_tpm = group.means(counts = atlas_tpm, groups = batch.labels)
x = batch_mean_tpm[, 1]
y = batch_mean_tpm[,2]
expr.cor = round(cor(x,y),2)
#smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("Before batch correction, correlation between \ntwo illustrative batches is %s", expr.cor),
#              xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")


print(head(atlas_tpm[1:10,1:10]))

# Compensate for batch effect using ComBat
# Takes a few minutes
atlas_tpm_norm = batch.normalise.comBat(counts = as.matrix(atlas_tpm), batch.groups = batch.labels)
batch_mean_tpm_norm = group.means(counts = atlas_tpm_norm, groups = batch.labels)
x = batch_mean_tpm_norm[, 1]
y = batch_mean_tpm_norm[,2]
expr.cor = round(cor(x,y),2)
#smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("After batch correction, correlation between \ntwo illustrative batches is %s", expr.cor),
#              xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")


#  Dimensionality reduction
# Run (randomized) PCA, t-SNE

##  Run the rpca really take long time
#pca = rpca(t(atlas_tpm_norm[var.genes,]), center=T, scale=T, retx=T, k=100)$x   #use selected variables for PCA

#barnes_hut_tsne = Rtsne(pca[, 1:13], check_duplicates=T,       pca=FALSE, #dont run PCA again
#                                                                                       initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)

#tsne.rot = barnes_hut_tsne$Y


# Test for significant PCs. To avoid very long runtimes, run on a high memory server with lots of cores (n.cores).
y = sig.pcs.perm(dat=t(atlas_umis[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)

if(file.exists('./data/atlas_pca_scores.txt')){
  pca = read.table("./data/atlas_pca_scores.txt")
}else{
  pca = rpca(t(atlas_tpm_norm[var.genes,]), center=T, scale=T, retx=T, k=100)$x   #use selected variables for PCA
  write.table(pca,'./data/atlas_pca_scores.txt',quote = F)
}
#pca = read.delim(file="regions_pca_scores.txt")



#   Test for sinificant principle compnents
y = sig.pcs.perm(dat=t(atlas_tpm_norm[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)  # 20 PC


if(file.exists('./data/atlas_tsne.txt')){
  tsne.rot = read.table("./data/atlas_tsne.txt")
}else{
  barnes_hut_tsne = Rtsne(pca[, 1:y$r], check_duplicates=T,pca=FALSE, #dont run PCA again
                          initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
  tsne.rot = barnes_hut_tsne$Y
  write.table(tsne.rot,'./data/atlas_tsne.txt',quote = F)
}


# Unsupervised clustering
## Run kNN-graph clustering
# build cell-cell euclidean distance matrix using significant PC scores

dm = as.matrix(dist(pca[, 1:13]))
# build nearest neighbor graph
knn = build_knn_graph(dm, k = 200)
clustering = cluster_graph(knn)$partition

# merge a spurious cluster (cluster 16 is only a single cell) into the most similar cluster
clustering = merge_clusters(clustering, c(8, 16))

## confirm that clusters are extremely similar to those in the paper (infomap is a random-walk based alg, so there may begetwd minor differences)
clusters_from_paper = factor(unlist(lapply(colnames(atlas_umis), get_field, 3,"_")))

overlap = as.data.frame.matrix(table(clusters_from_paper, clustering))
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)

jpeg(file='./data/Extend_data_1/aheatmap.jpeg')
aheatmap(overlap, color = cubehelix1.16, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
dev.off()





#  Visualize the clustering overlaid onto the t-SNE (Figure 1b)
x = data.frame(tsne.rot, clustering)
colnames(tsne.rot)<-c('tSNE_1','tSNE_2')
#ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=clustering)) + geom_point() + scale_color_manual(values=brewer16)

ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=cell_groups)) + geom_point() + scale_color_manual(values=brewer16)
ggsave(file='./data/Extend_data_1/figure_g.jpeg')


# Figure g
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
Genes_mean_tpm(stem_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Stem')


# Cell cycle
cell_cycle_mark_genes<-c('Mki67','Cdk4','Mcm5','Mcm6','Pcna')
Genes_mean_tpm(cell_cycle_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Cell cycle')


# Enterocyte
Enterocyte_mark_genes<-c('Alpi','Apoa1','Apoa4','Fabp1')
Genes_mean_tpm(Enterocyte_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Enterocyte')


# Globlet
Globlet_mark_genes<-c('Muc2','Clca1','Tff3','Agr2')
Genes_mean_tpm(Globlet_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Globet')

# Paneth
Paneth_mark_genes<-c('Lyz1','Defa17','Defa22','Defa24','Ang4')
Genes_mean_tpm(Paneth_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Paneth')


#  Enteroendocrine
Enteroendocrine_mark_genes<-c('Chga','Chgb','Tac1','Tph1','Neurog3')
Genes_mean_tpm(Enteroendocrine_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Enteroendocrine')



# Tuft
Tuft_mark_genes<-c('Dclk1','Trpm5','Gfi1b','Il25')
Genes_mean_tpm(Tuft_mark_genes,tpm_data = atlas_tpm_norm,tsne_data =tsne.rot,title = 'Tuft')


# reads,umis
per_cell_umis<-as.numeric(apply(atlas_umis[var.genes,],2,sum))
ggplot(tsne.rot, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=per_cell_umis))+theme(legend.title = element_text("UMIS/Cell",size=8,color='blue',face='bold'),
                                                                        legend.position = 'right') +scale_color_gradient2(low='lightblue',mid='green',high='red')



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

Genes_per_cell<-as.numeric(apply(atlas_tpm_norm,2,Count_genes))
ggplot(x, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Genes_per_cell))+theme(legend.title = element_text("Genes/Cell",size=8,color='blue',face='bold'),
                                                                              legend.position = 'right') +scale_color_gradient2(low='lightblue',mid='green',high='red')




