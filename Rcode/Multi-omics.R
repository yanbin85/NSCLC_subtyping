#安装与加载包
# devtools::install_local("/home/yanbin/R包/MOVICS-master.zip")
# devtools::install_github("xlucpu/MOVICS")
# # 或设置超时时间
# options(timeout = 6000)
# devtools::install_github("xlucpu/MOVICS")

library(MOVICS)

setwd("/home/yanbin/NSCLC_subtyping/MULTI-OMICS/")

# 读取并处理各组学数据 --------------------------------------------------------------
GEXP <- read.table("/home/data/yanbin/NSCLC_Dataset/LUNG_GEXP/TCGA.LUNG.sampleMap_HiSeqV2_PANCAN", header = TRUE, sep = "\t",row.names=1)
MUTA <- read.table("/home/yanbin/NSCLC_subtyping/MUTA/MUTA_NSCLC.csv", header = TRUE, sep = ",")
MET <- read.table("/home/yanbin/NSCLC_subtyping/DNA_MET/MET.csv", header = TRUE, sep = ",",row.names=1)
PROT<- read.table("/home/yanbin/NSCLC_subtyping/PROT/PROT.csv", header = TRUE, sep = ",",row.names=1)
TME <-  read.table("/home/yanbin/NSCLC_subtyping/TME/TME.csv", header = TRUE, sep = ",",row.names=1)
CN <- read.table("/home/yanbin/NSCLC_subtyping/CN/CN.csv", header = TRUE, sep=",",row.names=1)
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



# GET Module --------------------------------------------------------------
#基因表达数据处理
#z-score再次对每一行进行归一化
row_means <- apply(GEXP, 1, mean, na.rm = TRUE)
row_sds <- apply(GEXP, 1, sd, na.rm = TRUE)
row_sds[row_sds == 0] <- 1  # 若标准差为0，设为1（此时Z-Score为0）
GEXP_zscore <- sweep(GEXP, 1, row_means, "-")  # 减去均值
GEXP_zscore <- sweep(GEXP_zscore, 1, row_sds, "/")    # 除以标准差

#突变数据处理
anyNA(MUTA[,1])  # 返回TRUE表示存在NA
which(is.na(MUTA[,1]))
#[1] 7224
MUTA <- MUTA[!is.na(MUTA[,1]), ]
rownames(MUTA) <- MUTA[,1]
MUTA <- MUTA[,-1]

#DNA甲基化数据处理
# original_colnames <- colnames(MET) #提取原始列名
# is_B_sample <- substr(original_colnames, start = 16, stop = 16) == "B" #检查第16位是否为"B"（注意：字符串索引从1开始）
# MET_filtered <- MET[, !is_B_sample] #去除第16位为B的样本列
MET1 <- MET

#拷贝数数据二值化处理
# CN <- na.omit(CN) #去除缺失值
CN_binary <- CN
CN_binary[CN_binary == 2] <- 0   #拷贝数为2替换为0（无突变）
CN_binary[CN_binary != 0] <- 1   #拷贝数不为2替换为1（存在突变）

#TME数据处理
#转换为数值型
TME_numeric <- data.frame(lapply(TME, function(x) as.numeric(as.character(x))),
                         row.names = rownames(TME))
# original_colnames2 <- colnames(TME_numeric) #提取原始列名
# is_B_sample2 <- substr(original_colnames, start = 16, stop = 16) == "B" #检查第16位是否为"B"（注意：字符串索引从1开始）
# TME_filtered <- TME_numeric[, !is_B_sample] #去除第16位为B的样本列
#多组学共同样本处理
#统一各组学数据样本ID格式(截取列名前15位)
colnames(CN_binary) <- substr(colnames(CN_binary), start = 1, stop = 15)
colnames(TME_numeric) <- substr(colnames(TME_numeric), start = 1, stop = 15)
colnames(MET1) <- substr(colnames(MET1), start = 1, stop = 15)
colnames(PROT) <- substr(colnames(PROT), start = 1, stop = 15)

# 提取各数据集的样本ID
samples_GEXP <- colnames(GEXP_zscore)
samples_MUTA  <- colnames(MUTA)
samples_CNV  <- colnames(CN_binary)
samples_TME <- colnames(TME_numeric)
samples_MET <- colnames(MET1)
samples_PROT <- colnames(PROT)

# 将样本列表合并为一个列表
sample_list <- list(samples_GEXP, samples_MUTA, samples_CNV, 
                    samples_TME, samples_MET, samples_PROT)
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
# 筛选各数据集并严格对齐列顺序
GEXP_common <- GEXP_zscore[, common_samples]
MUTA_common <- MUTA[, common_samples]
CNV_common  <- CN_binary[, common_samples]
TME_common  <- TME_numeric[, common_samples]
MET_common  <- MET1[, common_samples]
PROT_common <- PROT[, common_samples]

# GEXP_common <- GEXP_zscore[,colnames(GEXP_zscore) %in% common_samples ]
# MUTA_common <- MUTA[,colnames(MUTA) %in% common_samples ]
# CNV_common <- CN_binary[,colnames(CN_binary) %in% common_samples ]
# TME_common <- TME_numeric[,colnames(TME_numeric) %in% common_samples ]
# MET_common <- MET[,colnames(MET1) %in% common_samples ]
# PROT_common <- PROT[,colnames(PROT) %in% common_samples ]


# 检查行数是否一致......
ncol(GEXP_common ) == ncol(MUTA_common)  # 应返回TRUE

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

elite.CNV <- getElites(dat       = CNV_common,
                        method    = "freq", # must set as 'freq'
                        na.action = "rm", #移除NA值
                        elite.pct = 0.73) #变异频率大于0.73（2854个）

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
#   $ CNV :'data.frame':	2854 obs. of  551 variables:
#   $ TME :'data.frame':	96 obs. of  551 variables:
#   $ MET :'data.frame':	3000 obs. of  551 variables:
#   $ PROT:'data.frame':	243 obs. of  551 variables:


#使用getClustNum()确定最佳亚型数量为4
optk.nsclc <- getClustNum(data        = multiomics.data, # 6种组学数据
                         is.binary   = c(F,T,T,F,F,F), #第二、三种数据是二分类的
                         try.N.clust = 2:8, # 尝试亚型数量，从2到8
                         fig.name    = "CLUSTER NUMBER OF NSCLC")#保存的文件名
#利用10种算法进行聚类
#1.iClusterBayes算法
iClusterBayes.res <- getiClusterBayes(data        = multiomics.data,
                                      N.clust     = 4,
                                      type        =  c("gaussian", "binomial","binomial","gaussian", "gaussian","gaussian" ),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5,0.5,0.5,0.5,0.5,0.5), #先验概率，控制每个子数据集特征选择的先验概率向量
                                      sdev        = 0.05,  #提议标准差,默认0.05
                                      thin        = 3)   #MCMC抽样间隔（默认3），用于减少自相关性。
#SNF算法
# getSNF(
#   data = multiomics.data,
#   N.clust = 4,
#   type= c("gaussian", "binomial","binomial","gaussian", "gaussian","gaussian" ),
#   K = 30,
#   t = 20,
#   sigma = 0.5
# )

#同时使用s剩下9种算法进行分型
moic.res.list <- getMOIC(data        = multiomics.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"), #9种算法
                         N.clust     = 4, #聚类数选择5
                         type        = c("gaussian", "binomial","binomial","gaussian", "gaussian","gaussian" ))
#整合iClusterBayes结果
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# 保存结果
save(moic.res.list, file = "moic.res.list.rda")

#绘制一致性热图
load(file = "moic.res.list.rda")
cmoic.nsclc <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
#计算Silhouette判断分型质量。
getSilhouette(sil      = cmoic.nsclc$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)
#多组学分型热图
# β值矩阵转换为M值矩阵
indata <- multiomics.data
indata$MET<- log2(indata$MET/ (1 - indata$MET))

# 对数据进行标准化
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,NA,NA,2,2,2), # no truncation for mutation
                     centerFlag = c(T,F,F,T,T,T), # no center for mutation
                     scaleFlag  = c(T,F,F,T,T,T)) # no scale for mutation



#保存工作空间
setwd("/home/yanbin/NSCLC_subtyping/MULTI-OMICS")
save.image("NSCLC_MULTIOMICS_Clustering.RData")