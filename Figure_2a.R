library(NMF)
library(rsvd)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)


#  Figure 1
### Load all the required functions for this analysis
source("Fxns.R")

## Downloading UMI count data
if(file.exists('./data/GSE92332_Regional_UMIcounts.txt.gz')){
  print('File exiists!')
}else{
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_Regional_UMIcounts.txt.gz", destfile="GSE92332_Regional_UMIcounts.txt.gz")

}
## Reading UMI count data from file
regions_umis = read.delim("./data/GSE92332_Regional_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(regions_umis), collapse = "x")))


#  Get variable genes(select variables)
#v = get.variable.genes_cvdiff(regions_umis, min.cv =0.15)


regions_umis<-regions_umis[which(rowSums(regions_umis)>0),]  # select expression genes

regions_tpm<-data.frame(log2(1+tpm(regions_umis)))
#regions_expr_tmp<-data.frame(log2(1+tpm(regions_umis_expr)))


# selected highly variable genes
v<-get.variable.genes(regions_umis,min.cv2 =100)

var.genes = as.character(rownames(v)[v$p.adj<0.05])


#  Batch correction (ComBat)
#Check whether there is a batch effect

get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
#batch.labels = factor(unlist(lapply(colnames(regions_umis), get_field,3,"_")))
#table(batch.labels)

## take mean tpm across batches to show batch effect
#batch_mean_tpm = group.means(counts = regions_tpm, groups = batch.labels)
#x = batch_mean_tpm[, 1]
#y = batch_mean_tpm[,2]
#expr.cor = round(cor(x,y),2)
#smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("Before batch correction, correlation between \ntwo illustrative batches is %s", expr.cor), 
#              xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")




# Compensate for batch effect using ComBat
# Takes a few minutes
#regions_tpm_norm = batch.normalise.comBat(counts = as.matrix(regions_tpm), batch.groups = batch.labels)
#batch_mean_tpm_norm = group.means(counts = regions_tpm_norm, groups = batch.labels)
#x = batch_mean_tpm_norm[, 1]
#y = batch_mean_tpm_norm[,2]
#expr.cor = round(cor(x,y),2)
#smoothScatter(x, y, nrpoints=Inf, pch=16, cex=0.25, main=sprintf("After batch correction, correlation between \ntwo illustrative batches is %s", expr.cor), 
#              xlab="All genes Batch 2, mean log2(TPM+1)", ylab="All genes Batch  1, mean log2(TPM+1)")


#  Dimensionality reduction
# Run (randomized) PCA, t-SNE

##  Run the rpca really take long time
if(file.exists('./data/regions_pca_scores.txt')){
   pca = read.table(".data/regions_pca_scores.txt")
}else{
   pca = rpca(t(regions_tpm[var.genes,]), center=T, scale=T, retx=T, k=100)$x   #use selected variables for PCA
   write.table(pca,'./data/regions_pca_scores.txt',quote = F)
}
#pca = read.delim(file="regions_pca_scores.txt")


if(file.exists('./data/regions_tsne.txt')){
   tsne.rot = read.table("./data/regions_tsne.txt")
}else{
   barnes_hut_tsne = Rtsne(pca[, 1:13], check_duplicates=T,pca=FALSE, #dont run PCA again
                     								initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
   tsne.rot = barnes_hut_tsne$Y
   write.table(tsne.rot,'./data/regions_tsne.txt',quote = F)
}



colnames(tsne.rot)<-c('t_SNE_1','t_SNE_2')
library(stringr)
regions<-unlist(lapply(rownames(pca),function(x)return(str_split(x,"_")[[1]][2])))
cell_types<-unlist(lapply(rownames(pca),function(x)return(str_split(x,"_")[[1]][4])))

tsne<-cbind(tsne.rot,data.frame(regions=regions,Cell_Types=cell_types))
#tsne<-cbind(tsne.rot,data.frame(clustering=regions_cluster))                                                   

ggplot(tsne,aes(x=t_SNE_1,y=t_SNE_2,color=regions))+geom_point()+ scale_color_manual(values=brewer16)+ggtitle('Regions')
ggsave('Regions_TSNE.jpeg')



ggplot(tsne,aes(x=t_SNE_1,y=t_SNE_2,color=cell_types))+geom_point()+ scale_color_manual(values=brewer16)+ggtitle('Cell_types')
ggsave('Cell_Type_SNE.jpeg')
 
#tsne.rot = read.delim("regions_tsne.txt")

# Test for significant PCs. To avoid very long runtimes, run on a high memory server with lots of cores (n.cores).
#y = sig.pcs.perm(dat=t(regions_umis[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)



# Unsupervised clustering
## Run kNN-graph clustering
# build cell-cell euclidean distance matrix using significant PC scores

#dm = as.matrix(dist(pca[, 1:13]))
# build nearest neighbor graph
dm = as.matrix(dist(pca[, 1:13]))
# build nearest neighbor graph
knn = build_knn_graph(dm, k = 200)
clustering = cluster_graph(knn)$partition

# merge a spurious cluster (cluster 16 is only a single cell) into the most similar cluster
clustering = merge_clusters(clustering, c(21))

## confirm that clusters are extremely similar to those in the paper (infomap is a random-walk based alg, so there may begetwd minor differences)

### cluster for regions
clusters_from_paper = factor(unlist(lapply(colnames(regions_umis), get_field, 2,"_")))

overlap = as.data.frame.matrix(table(clusters_from_paper, clustering))
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]

jpeg('./data/Regions_aheatmap.jpeg')
aheatmap(overlap, color = cubehelix1.16, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
#jpeg('./data/Regions_aheatmap.jpeg')
dev.off()

### cluster for cell types
clusters_from_paper = factor(unlist(lapply(colnames(regions_umis), get_field, 4,"_")))

overlap = as.data.frame.matrix(table(clusters_from_paper, clustering))
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]

jpeg('./data/Cell_type_aheatmap.jpeg')
aheatmap(overlap, color = cubehelix1.16, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
#jpeg('./data/Cell_type_aheatmap.jpeg')
dev.off()
