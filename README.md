# RNA-seq_DEG
## 1. 数据下载，检查
```R
# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)
setwd('/home/ylchen/AD/DEG/C1')
rawcount=read.table(file = 'counts19.txt',header = TRUE)
rawcount[1:4,1:4]
rownames(rawcount)=rawcount[,1]
rawcount=rawcount[,-1]

rawcount=subset(rawcount, select=-c(SRR1931812,SRR1931816))
colnames(rawcount)
#"SRR1931814" "SRR1931815" "SRR1931818" "SRR1931819"
# 过滤在至少在75%的样本中都有表达的基因
dim(rawcount)
keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)

filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

# 加载edgeR包计算counts per millio(cpm) 表达矩阵,并对结果取log2值
library(edgeR)
express_cpm <- log2(cpm(filter_count)+1)
express_cpm[1:4,1:4]

#group文献
group_list=c('AD','AD','Ctrl','Ctrl')
group_list=factor(group_list)

# 保存表达矩阵和分组结果
save(filter_count,express_cpm,group_list,file = "../data/step01_rawdata.rdata")

# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)
source('function.R')
load('../data/step01_rawdata.rdata')
#sample distribution
draw_pic(express_cpm = express_cpm,group_list = group_list)
```
![image](https://user-images.githubusercontent.com/78598588/110124166-bcefe100-7dfc-11eb-9c9b-bff4c932359a.png)
![image](https://user-images.githubusercontent.com/78598588/110124209-c5481c00-7dfc-11eb-8fdc-736fc94b3429.png)

## 2. 差异分析
### DESeq2
* 差异分析的function要根据实际修改**分组**及**阈值**
```R
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
```
![image](https://user-images.githubusercontent.com/78598588/110125616-7b603580-7dfe-11eb-810e-f614ac226cec.png)
