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
## Downloading UMI count data
if(file.exists("GSE92332_Regional_UMIcounts.txt.gz")){
  print('File exiists!')
}else{
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_Regional_UMIcounts.txt.gz", destfile="GSE92332_Regional_UMIcounts.txt.gz")
  
}
## Reading UMI count data from file
regional_umis = read.delim("GSE92332_Regional_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(regional_umis), collapse = "x")))


regional_umis_expr<-regional_umis[which(rowSums(regional_umis)>0),]  # select expression genes

regional_expr_tpm<-data.frame(log2(1+tpm(regional_umis_expr)))
#print(regional_expr_tpm[1:10,1:10])


cv_gene<-get.variable.genes_cvdiff(regional_expr_tpm,T)
gene.vars<-as.character(cv_gene$Genes)  # select genes


select_regional_tpm<-regional_expr_tpm[gene.vars,]
diffusion_matrix<-cbind(data.frame(Cell=colnames(select_regional_tpm)),t(select_regional_tpm))
rownames(diffusion_matrix)<-1:dim(diffusion_matrix)[1]


library(destiny)

library(stringr)
regions<-unlist(lapply(colnames(select_regional_tpm),function(x)return(str_split(x,"_")[[1]][2])))
cell_types<-unlist(lapply(colnames(select_regional_tpm),function(x)return(str_split(x,"_")[[1]][4])))

ct <- as.ExpressionSet(diffusion_matrix)
dif<-DiffusionMap(ct,verbose = T,rotate = T)
