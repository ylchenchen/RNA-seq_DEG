rm(list = ls()) 
options(stringsAsFactors = F)
source('function.R')
load('../data/step01_rawdata.rdata')

#DEG_DESeq2,change_id
DESeq2(filter_count = filter_count,group_list = group_list)

load('../data/Step03-DESeq2_nrDEG.Rdata')
#pheatmap
pheatmap(DEG_DESeq2,express_cpm)
#volcano
volcano(DEG_DESeq2)
