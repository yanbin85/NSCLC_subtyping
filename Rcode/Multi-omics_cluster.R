#安装与加载包
# devtools::install_local("/home/yanbin/R包/MOVICS-master.zip")
# devtools::install_github("xlucpu/MOVICS")


library(MOVICS)
library(dplyr)
library(readr)
library(tibble)
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata")

# 读取并处理各组学数据 --------------------------------------------------------------
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/GEXP.rda")
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/MUTA.rda")
MUTA <- MUTA_NSCLC
rm(MUTA_NSCLC)
gc()
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/MET.rda")
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/PROT.rda")
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/TME.rda")
TME <- NSCLC
rm(NSCLC)
gc()
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/inputdata/CN.rda")

#临床数据读取与合并
luad_clinical <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/PHENO/clinical/TCGA-LUAD.clinical.tsv",
                            header = TRUE,
                            sep = "\t",
                            row.names = 1,
                            quote = "",    # 禁用引号解析
                            fill = TRUE,   # 允许填充缺失列
                            check.names = FALSE)  # 保留列名中的特殊字符
lusc_clinical <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/PHENO/clinical/TCGA-LUSC.clinical.tsv",
                            header = TRUE,
                            sep = "\t",
                            row.names = 1,
                            quote = "",    # 禁用引号解析
                            fill = TRUE,   # 允许填充缺失列
                            check.names = FALSE)  # 保留列名中的特殊字符
# 获取两个数据框的列名
cols_luad <- colnames(luad_clinical)
cols_lusc <- colnames(lusc_clinical)
# 共有列（出现在两个数据集中的列）
common_cols <- intersect(cols_luad, cols_lusc)
cat("共有列数量:", length(common_cols), "\n")
print(common_cols)
# LUAD独有列（仅出现在LUAD中的列）
luad_only <- setdiff(cols_luad, cols_lusc)
cat("\nLUAD独有列:", length(luad_only), "\n")
print(luad_only)
# LUSC独有列（仅出现在LUSC中的列）
lusc_only <- setdiff(cols_lusc, cols_luad)
cat("\nLUSC独有列:", length(lusc_only), "\n")
print(lusc_only)
# 仅保留两个数据框共有的列
luad_common <- luad_clinical[, common_cols]
lusc_common <- lusc_clinical[, common_cols]
# 合并
clinical.data <- rbind(luad_common, lusc_common)

#生存数据
luad_survial <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/PHENO/survival/TCGA-LUAD.survival.tsv",
                            header = TRUE,
                            sep = "\t",
                            row.names = 1,
                            quote = "",    # 禁用引号解析
                            fill = TRUE,   # 允许填充缺失列
                            check.names = FALSE)  # 保留列名中的特殊字符
lusc_survial <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/PHENO/survival/TCGA-LUSC.survival.tsv",
                           header = TRUE,
                           sep = "\t",
                           row.names = 1,
                           quote = "",    # 禁用引号解析
                           fill = TRUE,   # 允许填充缺失列
                           check.names = FALSE)  # 保留列名中的特殊字符
#合并数据
survial.data <- rbind(luad_survial,lusc_survial)

# 示例数据:TCGA的乳腺癌数据
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))



# 分子筛选 --------------------------------------------------------------
#一、数据处理
#基因表达数据处理
z-score再次对每一行进行归一化
row_means <- apply(GEXP, 1, mean, na.rm = TRUE)
row_sds <- apply(GEXP, 1, sd, na.rm = TRUE)
row_sds[row_sds == 0] <- 1  # 若标准差为0，设为1（此时Z-Score为0）
GEXP_zscore <- sweep(GEXP, 1, row_means, "-")  # 减去均值
GEXP_zscore <- sweep(GEXP_zscore, 1, row_sds, "/")    # 除以标准差


#DNA甲基化数据处理
# original_colnames <- colnames(MET) #提取原始列名
# is_B_sample <- substr(original_colnames, start = 16, stop = 16) == "B" #检查第16位是否为"B"（注意：字符串索引从1开始）
# MET_filtered <- MET[, !is_B_sample] #去除第16位为B的样本列


#TME数据处理
#转换为数值型
TME_numeric <- data.frame(lapply(TME, function(x) as.numeric(as.character(x))),
                         row.names = rownames(TME))

#多组学共同样本处理
#统一各组学数据样本ID格式(截取列名前15位，.转换为-)
colnames(CN) <- substr(colnames(CN), start = 1, stop = 15)
colnames(TME_numeric) <- substr(colnames(TME_numeric), start = 1, stop = 15)
colnames(MET) <- substr(colnames(MET), start = 1, stop = 15)
colnames(PROT) <- substr(colnames(PROT), start = 1, stop = 15)
colnames(CN) <- gsub("\\.", "-", colnames(CN))
colnames(GEXP_zscore) <- gsub("\\.", "-", colnames(GEXP_zscore))
colnames(PROT) <- gsub("\\.", "-", colnames(PROT))
colnames(TME_numeric) <- gsub("\\.", "-", colnames(TME_numeric))

# 提取各数据集的样本ID
samples_GEXP <- colnames(GEXP_zscore)
samples_MUTA  <- colnames(MUTA)
samples_CNV  <- colnames(CN)
samples_TME <- colnames(TME_numeric)
samples_MET <- colnames(MET)
samples_PROT <- colnames(PROT)

# 将样本列表合并为一个列表
sample_list <- list(samples_GEXP, samples_MUTA, samples_CNV, 
                    samples_TME, samples_MET, samples_PROT)
#save(sample_list,file="sample_list.rda")

# 计算所有数据集的共有样本
common_samples <- Reduce(intersect, sample_list)
# 检查共有样本数量
length(common_samples)  
#[1] 551 #剩下551个样本
has_duplicates <- any(duplicated(common_samples))
if (has_duplicates) {
  print("存在重复样本ID！")
} else {
  print("所有样本ID唯一。")
}
#save(common_samples,file="common_saples.rda")
# 筛选各数据集并严格对齐列顺序
GEXP_common <- GEXP_zscore[, common_samples]
MUTA_common <- MUTA[, common_samples]
CNV_common  <- CN[, common_samples]
TME_common  <- TME_numeric[, common_samples]
MET_common  <- MET[, common_samples]
PROT_common <- PROT[, common_samples]


# 检查行数是否一致......
ncol(GEXP_common ) == ncol(MUTA_common)  # TRUE

# 检查样本顺序是否一致
all(colnames(GEXP_common) == colnames(MUTA_common))  # TRUE
all(colnames(GEXP_common) == colnames(CNV_common))  # TRUE
all(colnames(GEXP_common) == colnames(TME_common))  # TRUE
all(colnames(GEXP_common) == colnames(MET_common))  # TRUE
all(colnames(GEXP_common) == colnames(PROT_common))  # TRUE


elite.GEXP <- getElites(dat       = GEXP_common,
                         method    = "mad",
                         elite.num = 3000, # 保留MAD前3000的基因
                         na.action = "impute", # 使用knn进行插补
                         elite.pct = 0.2) #保留20%的数据（此时不起作用）
sample_num <- round(0.05*551) #至少在5%的样本中突变

elite.MUTA <- getElites(dat       = MUTA_common,
                       method    = "freq", # must set as 'freq'
                       elite.num = sample_num, #在sample_num个及以上样本中突变（619个基因）
                       elite.pct = 0.1) # 此时该参数不起作用

#方案一、拷贝数数据二值化处理以筛选分子
CNV_binary <- CNV_common
CNV_binary[CNV_binary == 2] <- 0   #拷贝数为2替换为0（无突变）
CNV_binary[CNV_binary != 0] <- 1   #拷贝数不为2替换为1（存在突变）
elite.CNV <- getElites(dat       = CNV_binary,
                          method    = "freq", # must set as 'freq'
                         na.action = "rm", #移除NA值
                         elite.pct = 0.73) #变异频率大于0.73（2854个）
# CNV.genes <- rownames(elite.CNV.binary$elite.dat)
# elite.CNV.data <- CNV_common[CNV.genes,]

#方案二、基于sd值筛选拷贝数数据
# elite.CNV <- getElites(dat       =CNV_common,
#                         method    = "sd",
#                         elite.num = 3000, # 保留sd前3000的基因
#                         na.action = "impute", # 使用knn进行插补
#                         ) 
#方案三、
# load("~/NSCLC_subtyping/MULTI-OMICS/data/CN_top3000.rda")
# colnames(CN_top3000) <- substr(colnames(CN_top3000), start = 1, stop = 15)
# CNV_common  <- CN_top3000[, common_samples]
# elite.CNV <- getElites(dat       =CNV_common,
#                          method    = "sd",
#                          elite.pct = 1, 
#                          na.action = "impute",
#                          ) 

#转换id
cnv.gene_mapping <- read_tsv("https://gdc-hub.s3.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap")
cnv.gene_mapping2 <- data.frame(id=cnv.gene_mapping$id,gene=cnv.gene_mapping$gene)
elite.CNV$elite.dat <- elite.CNV$elite.dat %>% 
tibble::rownames_to_column(var = "ensembl_id") %>%
left_join(cnv.gene_mapping2, by=c("ensembl_id"="id")) %>%
select(gene, everything()) 
#去除重复的gene并将其设为行名(剩余2796)
elite.CNV$elite.dat <- elite.CNV$elite.dat %>%
distinct(gene, .keep_all = TRUE) %>%  # 保留每个 gene 的第一行
column_to_rownames(var = "gene")
elite.CNV$elite.dat <- elite.CNV$elite.dat[,-1]

elite.TME <- getElites( dat       = TME_common,
                        method    = "mad",
                        elite.pct = 1, # 保留所有特征（96个）
                        na.action = "impute", # 使用knn进行插补
                       ) 

elite.MET <- getElites(dat       = MET_common,
                       method    = "mad",
                       elite.num = 3000, # 保留MAD前3000的位点
                       na.action = "impute", # 使用knn进行插补
                       ) 


elite.PROT <- getElites(dat       = PROT_common,
                        method    = "mad",
                        elite.pct = 0.5, # 保留50%蛋白（243个）
                        na.action = "impute", # 使用knn进行插补
                        ) 

# 提取每个列表中的 elite.dat 数据框
GEXP_dat <- elite.GEXP$elite.dat  
MUTA_dat <- elite.MUTA$elite.dat 
CNV_dat  <- elite.CNV$elite.dat  
TME_dat  <- elite.TME$elite.dat   
MET_dat  <- elite.MET$elite.dat  
PROT_dat <- elite.PROT$elite.dat  

# #甲基化数据位点转换为基因
# met.gene_mapping <- read_tsv("https://gdc-hub.s3.us-east-1.amazonaws.com/download/HM450.hg38.manifest.gencode.v36.probeMap")
# met.gene_mapping2 <- data.frame(id=met.gene_mapping$`#id`,gene=met.gene_mapping$gene)
# MET_dat<- as.data.frame(MET_dat) %>% 
#   tibble::rownames_to_column(var = "id") %>%
#   left_join(met.gene_mapping2, by=c("id"="id")) %>%
#   select(gene, everything()) 

multiomics.data <- list(
  GEXP = GEXP_dat,  # 基因表达
  MUTA = MUTA_dat,  # 突变
  CNV  = CNV_dat,   # 拷贝数
  TME  = TME_dat,   # 肿瘤微环境
  MET  = MET_dat,   # 甲基化
  PROT = PROT_dat   # 蛋白
)

# 检查结构
str(multiomics.data, max.level = 1)
# List of 6
# $ GEXP:'data.frame':	3000 obs. of  551 variables:
#   $ MUTA:'data.frame':	619 obs. of  551 variables:
#   $ CNV :'data.frame':	2796 obs. of  551 variables:
#   $ TME :'data.frame':	96 obs. of  551 variables:
#   $ MET :'data.frame':	3000 obs. of  551 variables:
#   $ PROT:'data.frame':	243 obs. of  551 variables:




# 获得共识分型 ------------------------------------------------------------------
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS")
#使用getClustNum()确定最佳亚型数量为5
optk.nsclc <- getClustNum(data        = multiomics.data, # 6种组学数据
                         is.binary   = c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE), #第二种数据是二分类的
                         try.N.clust = 2:8, # 尝试亚型数量，从2到8
                         fig.name    = "CLUSTER NUMBER OF NSCLC")#保存的文件名

save(optk.nsclc,file="optk.nsclc.rda")
#利用10种算法进行聚类
#iClusterBayes算法
iClusterBayes.res <- getiClusterBayes(data        = multiomics.data,
                                      N.clust     = 5,
                                      type        =  c("gaussian", "binomial","binomial","gaussian", "gaussian","gaussian" ),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5,0.5,0.5,0.5,0.5,0.5), #先验概率，控制每个子数据集特征选择的先验概率向量
                                      sdev        = 0.05,  #提议标准差,默认0.05
                                      thin        = 3)   #MCMC抽样间隔（默认3），用于减少自相关性。


#同时使用s剩下9种算法进行分型
moic.res.list <- getMOIC(data        = multiomics.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"), #9种算法
                         N.clust     = 5, #聚类数选择4
                         type        = c("gaussian", "binomial","binomial","gaussian", "gaussian","gaussian" ))

#整合iClusterBayes结果并保存
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/outputdata")
save(moic.res.list, file = "moic.res.list.rda")#更改


# 可视化 ---------------------------------------------------------------------
#绘制一致性热图
load(file = "moic.res.list.rda")
#整合10种算法分型结果获得共识亚型
cmoic.nsclc <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

save(cmoic.nsclc,file="cmoic.cluster.results.rda")#更改
#计算Silhouette判断分型质量
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS")
getSilhouette(sil      = cmoic.nsclc$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

#绘制多组学综合热图
# β值矩阵转换为M值矩阵
indata <- multiomics.data
indata$MET<- log2(indata$MET/ (1 - indata$MET))

# 对数据进行标准化
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,NA,NA,2,2,2), # no truncation for mutation
                     centerFlag = c(T,F,F,T,T,T), # no center for mutation
                     scaleFlag  = c(T,F,F,T,T,T)) # no scale for mutation
# 检查每个数据集中的缺失值比例
sapply(plotdata, function(x) sum(is.na(x)) / length(x))
#去除缺失值
plotdata$MET <- na.omit(plotdata$MET)
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/outputdata")
save(multiomics.data,file="multiomics.data.rda")
save(plotdata,file="normalized.multiomics.data.rda")#更改

# # 替换Inf/NaN为合理值（如表达数据）
# plotdata$GEXP[is.infinite(plotdata$GEXP)] <- max(plotdata$GEXP[is.finite(plotdata$GEXP)], na.rm=TRUE)
# plotdata$GEXP[is.nan(plotdata$GEXP)] <- median(plotdata$GEXP, na.rm=TRUE)
# 基于iClusterBayes的结果在每个组学中选择前10个分子进行标注
feat<- iClusterBayes.res$feat.res
feat.GEXP <- feat[which(feat$dataset == "GEXP"),][1:10,"feature"] 
feat.MUTA <- feat[which(feat$dataset=="MUTA"),][1:10,"feature"] 
feat.CNV <- feat[which(feat$dataset=="CNV"),][1:10,"feature"] 
feat.TME <- feat[which(feat$dataset=="TME"),][1:10,"feature"] 
feat.MET <- feat[which(feat$dataset=="MET"),][1:10,"feature"] 
feat.PROT <- feat[which(feat$dataset=="PROT"),][1:10,"feature"] 




annRow <- list(feat.GEXP, feat.MUTA, feat.CNV, feat.TME,feat.MET,feat.PROT)

# 为每个组学的热图自定义颜色
GEXP.col <- c("#008000","#00FF00", "white","#FF0000", "#800000" )
MUTA.col <- c("grey90" , "black")
CNV.col <- c( "#FFFFFF","#FF9999" )
TME.col <- c("#6699CC", "white"  , "#FF3C38")
MET.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
PROT.col <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")

# mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
# lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
# meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
# mut.col    <- c("grey90" , "black")
# col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

col.list   <- list(GEXP.col,MUTA.col,CNV.col,TME.col,MET.col,PROT.col)

# 绘制iClsuerBayes热图
 setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS")
getMoHeatmap(data          = plotdata,
             row.title     = c("gene expression","somatic mutation","copy number variation","tumor microenvironment","DNA methylation","protein expression"),
             is.binary     = c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE), 
             legend.name   = c("mRNA.counts","mutated","varied","tme","M value","rppa value"),
             show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             show.colnames = FALSE, # show no sample names
             show.col.dend = TRUE,
             annRow        = annRow, # mark selected features
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    =iClusterBayes.res$clust.dend, # show dendrogram for samples
             clust.col = c("#2CA02C",  "#FF9F1C", "#FFA5AB",  "#9D4EDD","#118AB2"),
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")



# 绘制共识热图
getMoHeatmap(data          = plotdata,
             row.title     = c("gene expression","somatic mutation","copy number variation","tumor microenvironment","DNA methylation","protein expression"),
             is.binary     = c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE), 
             legend.name   = c("mRNA.counts","mutated","varied","tme","M value","rppa value"),
             show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             show.colnames = FALSE, # show no sample names
             show.col.dend = FALSE,
             annRow        = NULL, # mark selected features
             clust.res     = cmoic.nsclc$clust.res, # cluster results
             clust.dend    = cmoic.nsclc$clust.dend, # show dendrogram for samples
             clust.col = c("#2CA02C",  "#FF9F1C", "#FFA5AB",  "#9D4EDD","#118AB2"),
             show.row.dend = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ConsensusCluster")


getMoHeatmap(data          = plotdata,
             row.title     = c("gene expression","somatic mutation","copy number variation","tumor microenvironment","DNA methylation","protein expression"),
             is.binary     = c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE), 
             legend.name   = c("mRNA.counts","mutated","varied","tme","M value","rppa value"),
             show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
             show.colnames = FALSE, # show no sample names
             show.col.dend = FALSE,
             annRow        = annRow, # mark selected features（iCluster Bayes）
             clust.res     = cmoic.nsclc$clust.res, # cluster results
             clust.dend    = cmoic.nsclc$clust.dend, # show dendrogram for samples
             clust.col = c("#2CA02C",  "#FF9F1C", "#FFA5AB",  "#9D4EDD","#118AB2"),
             show.row.dend = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ConsensusCluster(iCluster Bayes makers)")


#保存工作空间
setwd("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS")
save.image("NSCLC_MULTIOMICS_Clustering.RData")
