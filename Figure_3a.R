library(NMF)
library(rsvd)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)
library(stringr)

#  Figure 1
### Load all the required functions for this analysis
source("Fxns.R")

## Downloading UMI count data
if(file.exists('./data/GSE92332_EEC_UMIcounts.txt.gz')){
  print('File exiists!')
}else{
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_EEC_UMIcounts.txt.gz", destfile="GSE92332_EEC_UMIcounts.txt.gz")

}
## Reading UMI count data from file
EEC_umis = read.delim("./data/GSE92332_EEC_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(EEC_umis), collapse = "x")))
EEC_tpm<-data.frame(log2(1+tpm(EEC_umis)))



#  Get variable genes(select variables)
v = get.variable.genes(EEC_umis, min.cv2 = 100)
var.genes = as.character(rownames(v)[v$p.adj<0.05])


get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
##  Run the rpca really take long time
if(file.exists('./data/EEC_pca_scores.txt')){
   pca = read.table("./data/EEC_pca_scores.txt")
}else{
   pca = rpca(t(EEC_tpm[var.genes,]), center=T, scale=T, retx=T, k=100)$x   #use selected variables for PCA
   write.table(pca,'./data/EEC_pca_scores.txt',quote = F)
}
#pca = read.delim(file="regions_pca_scores.txt")



#   Test for sinificant principle compnents
y = sig.pcs.perm(dat=t(EEC_umis[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20, randomized=T)


if(file.exists('./data/EEC_tsne.txt')){
   tsne.rot = read.table("./data/EEC_tsne.txt")
}else{
   barnes_hut_tsne = Rtsne(pca[, 1:13], check_duplicates=T,pca=FALSE, #dont run PCA again
                                                                                initial_dims = 13, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
   tsne.rot = barnes_hut_tsne$Y
   write.table(tsne.rot,'./data/EEC_tsne.txt',quote = F)
}


sample_names<-colnames(EEC_umis)
clusters<-unlist(lapply(sample_names,function(x){
  return(str_split(x,"_",n=4)[[1]][4])
}))


colnames(tsne.rot)<-c('t_SNE_1','t_SNE_2')

tsne<-cbind(tsne.rot,data.frame(clusters=clusters))

ggplot(tsne,aes(x=t_SNE_1,y=t_SNE_2,color=clusters))+geom_point()+ scale_color_manual(values=brewer16)+ggtitle('EEC')
ggsave('./data/EEC_TSNE.jpeg')



# build nearest neighbor graph
dm = as.matrix(dist(pca[, 1:13]))
# build nearest neighbor graph
knn = build_knn_graph(dm, k = 200)
clustering = cluster_graph(knn)$partition

# merge a spurious cluster (cluster 16 is only a single cell) into the most similar cluster
clustering = merge_clusters(clustering, c(21))

## confirm that clusters are extremely similar to those in the paper (infomap is a random-walk based alg, so there may begetwd minor differences)

### cluster for regions
#clusters_from_paper = factor(unlist(lapply(colnames(regions_umis), get_field, 2,"_")))

overlap = as.data.frame.matrix(table(clusters, clustering))
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]

jpeg('./data/EEC_aheatmap.jpeg')
aheatmap(overlap, color = cubehelix1.16, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
dev.off()




