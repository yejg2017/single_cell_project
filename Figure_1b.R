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
if(file.exists('./data/GSE92332_atlas_UMIcounts.txt.gz')){
  print('File exiists!')
}else{
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz", destfile="GSE92332_atlas_UMIcounts.txt.gz")

}
## Reading UMI count data from file
atlas_umis = read.delim("./data/GSE92332_atlas_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(atlas_umis), collapse = "x")))


#  Get variable genes(select variables)
v = get.variable.genes(atlas_umis, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])


#  Batch correction (ComBat)
#Check whether there is a batch effect

get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
batch.labels = factor(unlist(lapply(colnames(atlas_umis), get_field, 1,"_")))
print(table(batch.labels))

atlas_tpm = data.frame(log2(1+tpm(atlas_umis)))
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
pca = read.delim(file="atlas_pca_scores.txt")

#barnes_hut_tsne = Rtsne(pca[, 1:13], check_duplicates=T,	pca=FALSE, #dont run PCA again
#                        								initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)

#tsne.rot = barnes_hut_tsne$Y

tsne.rot = read.delim("atlas_tsne.txt")

# Test for significant PCs. To avoid very long runtimes, run on a high memory server with lots of cores (n.cores).
y = sig.pcs.perm(dat=t(atlas_umis[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)



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
overlap = overlap[,apply(overlap, 1, FUN=which.max)]
aheatmap(overlap, color = cubehelix1.16, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)





#  Visualize the clustering overlaid onto the t-SNE (Figure 1b)
x = data.frame(tsne.rot, clustering)
ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=clustering)) + geom_point() + scale_color_manual(values=brewer16)





