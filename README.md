# single_cell_preoject
A single-cell survey of the small intestinal epithelium.  Realized lots of plot job in the parper

A single-cell survey of the small intestinal epithelium
            文章重现总结

一．环境配置：
Linux环境;R软件，需要用各种依赖包：
ggplot2 ;MAST ;NMF ;rsvd ;Rtsne ;cowplot ;igraph ;xlsx ;destiny ;DESeq2 ;devtools ;fpc ;goseq;doParallel;foreign;data.table;cccd;car;jackstraw ;KernelSmooth ;markdown ;FNN ;
sva;densityCluster;pvclust;foreach;lfa ;knitr ;RMySQL ;stringr ;stringi ;difussionMap;
optparse

安装命令:
source('http://bioconductor.org/biocLite.R')
例如：biocLite('destiny')

或者 install.packages(…)

在安装这些依赖包过程中,由于是在linux环境下配置，会出现各种环境配置，权限等各种问题，所以使用docker容器配置：bioconductor/devel_core2(docker pull bioconductor/devel_core2可安装);提供了Rstudio界面,方便敲代码。

二．图像实现说明：

参考了作者提供的一部分代码:https://github.com/yejg2017/single_cell_intestine(仅仅是提供了如何实现tSNE和变量筛选的代码)

使用的数据集:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92332

文章中很多图像的实现都无法完全重现(结果不太一样),主要的原因是缺少数据处理的很多细节，以及数据集没有能够提供图像实现的足够信息;但是该项目基本都有提供了各种图像实现的方法以及详细的代码(PS:只是重现了文章RESEARCH后的Extended Data Figure1-10的图像)。由于在数据方面很多细节方法可能还需要仔细考虑，所以目前项目重现所展现的只是方法以及代码上实现过程，缺少很多结果以及各种详细的说明。所有代码以及实现的图像效果都以markdown(生成html文件格式）的形式展现。

所有实现的过程(代码，数据集，结果数据，markdown)都已经提交到个人Github
上:https://github.com/yejg2017/single_cell_preoject
