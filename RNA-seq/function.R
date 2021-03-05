#sample distribution
draw_pic=function(express_cpm,group_list){
#样本表达总体分布-箱式图
library(ggplot2)
exprSet <- express_cpm
data <- data.frame(expression=c(exprSet),sample=rep(colnames(exprSet),each=nrow(exprSet)))
#head(data)
p <- ggplot(data = data,aes(x=sample,y=expression,fill=sample))
p1 <- p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) + xlab(NULL) + ylab("log2(CPM+1)")
p1
# 保存图片
pdf(file = "../pic/sample_boxplot.pdf",width = 6,height = 8)
print(p1)
dev.off()

# 样本表达总体分布-小提琴图
p2 <- p + geom_violin() + theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle = 90)) + xlab(NULL) + ylab("log2(CPM+1)")
p2
# 保存图片
pdf(file = "../pic/sample_violin.pdf",width = 6,height = 8)
print(p2)
dev.off()

# 样本表达总体分布-概率密度分布图
m <- ggplot(data=data, aes(x=expression))
p3 <- m +  geom_density(aes(fill=sample, colour=sample),alpha = 0.2) + xlab("log2(CPM+1)")
p3
# 保存图片
pdf(file = "../pic/sample_density.pdf",width = 7,height = 8)
print(p3)
dev.off()

## 1.样本之间的相关性-层次聚类树----
dat <- express_cpm
sampleTree <- hclust(dist(t(dat),method = "manhattan"), method = "average")
p=plot(sampleTree)
p
print(sampleTree)
dev.off()

## 2.样本之间的相关性-PCA----
# 第一步，数据预处理
dat <- as.data.frame(t(dat))
dat$group_list <- group_list
library(FactoMineR)
library(factoextra)

# 画图仅需要数值型数据，去掉最后一列的分组信息
dat_pca <- PCA(dat[,-ncol(dat)], graph = FALSE)

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "text", # 只显示点，不显示文字
                  col.ind = dat$group_list, # 用不同颜色表示分组
                  palette = c("#00AFBB", "#E7B800"),
                  addEllipses = T, # 是否圈起来
                  legend.title = "Groups")
ggsave(p,filename = '../pic/sample_pca.png',width = 10,height = 7)


## 3.样本之间的相关性-cor----
# 利用绝对中位差mad/标准差sd统计学方法进行数据异常值检测
# 将表达量的绝对中位差mad从大到小排列取前500的结果
dat <- express_cpm
tmp <- sort(apply(dat,1, mad),decreasing = T)[1:500]
exprSet <-dat[names(tmp),]
# 使用500个基因的表达量来做相关性图
library(corrplot)
#dim(exprSet)

# 计算相关性
M <- cor(exprSet)
g <- corrplot(M,order = "AOE",addCoef.col = "white")

#corrplot(M,order = "AOE",type="upper",tl.pos = "d")
#corrplot(M,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n")
# 绘制样本相关性的热图
anno <- data.frame(sampleType=group_list)
rownames(anno) <- colnames(exprSet)
#anno
p=pheatmap::pheatmap(M,display_numbers = T,annotation_col = anno,fontsize = 11,cellheight = 28,cellwidth = 28)
# 保存图片
pdf(file = "../pic/sample_pheatmap.pdf",width = 7,height = 8)
print(p)
dev.off()
}


#DEG_DESeq2
DESeq2=function(filter_count,group_list){
  exprSet <- filter_count
  dim(exprSet)
  exprSet[1:4,1:4]
  # 加载包
  library(DESeq2)
  # 第一步，构建DESeq2的DESeq对象
  colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ group_list)
  # 第二步，进行差异表达分析
  dds2 <- DESeq(dds)
  # 提取差异分析结果，trt组对untrt组的差异分析结果
  tmp <- results(dds2,contrast=c("group_list","AD","Ctrl"))
  DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
  # 去除差异分析结果中包含NA值的行
  DEG_DESeq2 = na.omit(DEG_DESeq2)
  # 筛选上下调，设定阈值
  fc_cutoff <- 1.5
  fdr <- 0.05
  DEG_DESeq2$regulated <- "normal"
  
  loc_up <- intersect(which(DEG_DESeq2$log2FoldChange>=log2(fc_cutoff)),which(DEG_DESeq2$padj<fdr))
  loc_down <- intersect(which(DEG_DESeq2$log2FoldChange< (-log2(fc_cutoff))),which(DEG_DESeq2$padj<fdr))
  
  DEG_DESeq2$regulated[loc_up] <- "up"
  DEG_DESeq2$regulated[loc_down] <- "down"
  
  print(table(DEG_DESeq2$regulated))
  # 注释,添加一列gene symbol
  #BiocManager::install('org.Hs.eg.db')
  library(org.Hs.eg.db)
  keytypes(org.Hs.eg.db)
  
  library(clusterProfiler)
  id2symbol <- bitr(rownames(DEG_DESeq2), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
  
  symbol <- rep("NA",time=nrow(DEG_DESeq2))
  symbol[match(id2symbol[,1],rownames(DEG_DESeq2))] <- id2symbol[,2]
  DEG_DESeq2 <- cbind(rownames(DEG_DESeq2),symbol,DEG_DESeq2)
  colnames(DEG_DESeq2)[1] <- "GeneID"
  #删除无注释的基因
  print(table(DEG_DESeq2$symbol!='NA'))
  DEG_DESeq2=DEG_DESeq2[DEG_DESeq2$symbol!='NA',]
  # 保存
  #write.table(DEG_DESeq2,"data/DEG_DESeq2_all.csv",row.names = F,sep="\t",quote = F)
  
  ## 取表达差异倍数和p值,矫正后的pvalue并保存
  DEG_DESeq2 <- DEG_DESeq2[,c(1,2,4,7,8,9)]
  print(dim(DEG_DESeq2))
  print(table(DEG_DESeq2$regulated))
  save(DEG_DESeq2, file = "../data/Step03-DESeq2_nrDEG.Rdata")
  }

#pheatmap
pheatmap=function(DEG_DESeq2,express_cpm){
  DESeq2_sigGene <- DEG_DESeq2[DEG_DESeq2$regulated!="normal",1]
dat <- express_cpm[match(DESeq2_sigGene,rownames(express_cpm)),]
group <- data.frame(group=colnames(dat))
rownames(group)=colnames(dat)
library(pheatmap)
p <- pheatmap::pheatmap(dat,scale = "row",show_colnames =T,show_rownames = T, cluster_cols = T,annotation_col=group)
p
ggsave(p,filename = '../pic/DEG_pheatmap.png')
}

#volcano
volcano=function(DEG_DESeq2){
  logFC_cutoff=1.5
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG_DESeq2[DEG_DESeq2$regulated =='up',]) ,
                    '\nThe number of down gene is ',nrow(DEG_DESeq2[DEG_DESeq2$regulated =='down',]))

library(ggplot2)
g <- ggplot(data=DEG_DESeq2, 
            aes(x=log2FoldChange, y=-log10(padj), 
                color=regulated)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g)
ggsave(g,filename = '../pic/volcano.png')}



