install.packages("ComplexHeatmap")
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(dplyr)
library(tibble)
library(readr)
setwd("/home/yanbin/NSCLC_subtyping/CN/")
#LUAD与LUSC数据读取与合并
CN_LUAD <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/CN/TCGA-LUAD.gene-level_ascat3.tsv", header = TRUE, sep = "\t",row.names=1)
CN_LUSC <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/CN/TCGA-LUSC.gene-level_ascat3.tsv", header = TRUE, sep = "\t",row.names=1)
common_genes <- intersect(rownames(CN_LUAD), rownames(CN_LUSC))
CN <- cbind(CN_LUAD[common_genes, ], CN_LUSC[common_genes, ])
CN<- cbind(CN_LUAD,CN_LUSC)
write.csv(CN,"CN.csv")
# NSCLC聚类一 ----------------------------------------------------------------
# 数据过滤与高变基因筛选 
#过滤掉包含NA值的基因

CN_filtered <- CN %>% 
  filter_all(all_vars(!is.na(.)))
dim(CN_filtered)  # 查看最终维度
sum(is.na(CN_filtered))

# 先计算总样本数
total_samples <- ncol(CN_filtered) - 1  

# 再添加变异统计列
CN_filtered$variant_count <- apply(CN_filtered[, -1], 1, function(x) sum(x != 2))
CN_filtered$variant_freq <- CN_filtered$variant_count / total_samples
hist(CN_filtered$variant_freq,breaks=100)
# 筛选高频基因（88个）
#high_freq_genes <- CN_filtered%>% 
  #filter(variant_freq >= 0.795) %>% 
  #arrange(desc(variant_freq))

high_freq_genes <- CN_filtered[order(CN_filtered$variant_freq, decreasing = TRUE), ][1:1000, ]

# 提取原始样本数据列（排除variant_count和variant_freq）
high_freq_genes <- high_freq_genes[, colnames(CN)]

# 转换基因标识符
gene_mapping <- read_tsv("https://gdc-hub.s3.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap")
gene_mapping2 <- data.frame(id=gene_mapping$id,gene=gene_mapping$gene)
#high_freq_genes <- high_freq_genes %>% 
  #tibble::rownames_to_column(var = "ensembl_id") %>%
  #left_join(gene_mapping2, by=c("ensembl_id"="id")) %>%
  #select(gene, everything()) 
#去除重复的gene并将其设为行名
#high_freq_genes <- high_freq_genes %>%
  #distinct(gene, .keep_all = TRUE) %>%  # 保留每个 gene 的第一行
  #column_to_rownames(var = "gene")
#high_freq_genes <- high_freq_genes[,-1]



#mads <-apply(CN_filtered, 1, mad, na.rm = TRUE) # 计算每一行的MAD值
#CN_top5000 <- CN_filtered[rev(order(mads))[1:5000], ] # 提取前5000个基因

#hist(mads, breaks = 100)
#sds <-apply(CN_filtered, 1, sd, na.rm = TRUE) # 计算每一行的MAD值
#hist(sds, breaks = 100)

# 热图绘制 
setwd("/home/yanbin/NSCLC_subtyping/CN/NSCLC")
#分组
# 肺腺癌的项目代码列表
luad_codes <- c(
  "05" ,"17","35","38","44","49","4B","50","53","55","62",
  "64","67","69","71","73","75","78","80","83","86","91",
  "93","95","97","99","J2","L4","L9","ME","MN","MP","NB",
  "NJ","O1","S2","T6"
)
# 肺鳞癌的样本代码
lusc_codes <- c(
  "11", "18", "21", "22", "33", "34", "37", "39", "43", "46", 
  "51", "52", "56", "58", "60", "63", "66", "68", "6A", "70", 
  "77", "79", "82", "85", "90", "92", "94", "96", "98", 
  "J1", "L3", "LA", "MF", "ML", "NC", "NK", "O2", "OC", "UJ", "WG", "XC", "ZE"
)
#heatmap_input <- high_freq_genes[,1:(ncol(high_freq_genes)-2)]
# 提取样本条形码第5-6位代码
sample_codes <- substr(colnames(high_freq_genes), 6, 7)

# 创建分组注释向量（LUAD=蓝色，LUSC=红色，其他样本=灰色）
sample_groups <- case_when(
  sample_codes %in% luad_codes ~ "LUAD",
  sample_codes %in% lusc_codes ~ "LUSC",
)
sample_groups <- factor(sample_groups, levels = c("LUAD", "LUSC"))
#定义分组颜色
group_colors <- c("LUAD" = "#377EB8", "LUSC" = "#E41A1C")


#绘制热图
library(circlize)
col_fun <- colorRamp2(
  c(0, 2, 8), 
  c("blue", "white", "red")
)

# 添加注释并分割列
column_ha <- HeatmapAnnotation(
  group = sample_groups,
  col = list(group = group_colors),
  annotation_name_side = "left"
)

pdf("CN_top1000_heatmap.pdf",width = 12,height = 9)
Heatmap(
  as.matrix(high_freq_genes),
  name = "Copy Number",
  col = col_fun,
  top_annotation = column_ha,
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  show_column_names  = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
)
dev.off()

#第2种颜色映射
col_fun2 <- colorRamp2(c(0, 2, 4, 6, 8), c("blue", "white", "yellow", "orange", "red"))
pdf("CN_top100_heatmap2.pdf",width = 12,height = 9)
Heatmap(
  as.matrix(heatmap_input),
  name = "Copy Number",
  col = col_fun2,
  top_annotation = column_ha,
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  show_column_names  = FALSE,
  row_names_gp = gpar(fontsize = 6),
  cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
)
dev.off()
#热图绘制






# ConsensusClusterPlus聚类 
setwd("/home/yanbin/NSCLC_subtyping/CN/NSCLC")
#筛选聚类所用高变基因（突变频率前百分之10，5794个）
high_freq_genes2 <- CN_filtered[rev(order(CN_filtered$variant_freq))[1:(0.1*nrow(CN_filtered))], ] 
# high_freq_genes2 <- CN_filtered%>% 
#   filter(variant_freq >= 0.75) %>% 
#   arrange(desc(variant_freq))
input_data<-as.matrix(high_freq_genes2[,1:(ncol(high_freq_genes2)-2)])
title=tempdir()
results1 = ConsensusClusterPlus(input_data, 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_CN_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "euclidean", #改为欧式距离（否则会提示方差为0？）
                                innerLinkage="ward.D2",
                                finalLinkage="ward.D2",
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results1, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]


#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:10
x1 = 0.1; x2 = 0.9       
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  
for(i in Kvec){
  M = results1[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK1 = Kvec[which.min(PAC)]  # 理想的K值为10





# NSCLC聚类二 ----------------------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/CN/NSCLC2")
#对原始拷贝数数据进行PCA
pca_result <- prcomp(t(CN_filtered), scale = TRUE)  # scale=TRUE标准化数据

# 选择主成分（累计解释>80%方差的主成分数）
summary(pca_result)  # 查看方差贡献率
cum_var <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_pcs <- which(cum_var >= 0.8)[1]  # 选择累计方差≥80%的主成分数

# 提取主成分作为新特征
pca_data <- pca_result$x[, 1:n_pcs]

# 基于主成分进行共识聚类
input_data2<-t(pca_data)
title=tempdir()
results2 = ConsensusClusterPlus(input_data2, 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_CN_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "euclidean", 
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results1, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]
#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:10
x1 = 0.1; x2 = 0.9       
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  
for(i in Kvec){
  M = results2[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK2 = Kvec[which.min(PAC)]  # 理想的K值为4

# LUAD聚类 ------------------------------------------------------------------
library(dplyr)
LUAD_filtered <- CN_LUAD %>% 
  filter_all(all_vars(!is.na(.)))
dim(LUAD_filtered)  # 查看最终维度
sum(is.na(LUAD_filtered))

# 先计算总样本数
total_samples2 <- ncol(LUAD_filtered) - 1  

# 再添加变异统计列
LUAD_filtered$variant_count <- apply(LUAD_filtered[, -1], 1, function(x) sum(x != 2))
LUAD_filtered$variant_freq <- LUAD_filtered$variant_count / total_samples

hist(LUAD_filtered$variant_freq,breaks=100)
# 筛选高频基因
high_freq_genesLUAD <- LUAD_filtered%>% 
  filter(variant_freq >= 0.35) %>% 
  arrange(desc(variant_freq))

setwd("/home/yanbin/NSCLC_subtyping/CN/LUAD")
title=tempdir()
results1 = ConsensusClusterPlus(as.matrix(high_freq_geneLUAD), 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_CNtop5000_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "euclidean", 
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")


# LUSC聚类 ------------------------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/CN/LUSC")
LUSC_filtered <- CN_LUAD %>% 
  filter_all(all_vars(!is.na(.)))
dim(LUSC_filtered)  # 查看最终维度
sum(is.na(LUSC_filtered))

# 先计算总样本数
total_samplesLUSC <- ncol(LUSC_filtered) - 1  

# 再添加变异统计列
LUSC_filtered$variant_count <- apply(LUSC_filtered[, -1], 1, function(x) sum(x != 2))
LUSC_filtered$variant_freq <- LUSC_filtered$variant_count / total_samplesLUSC
hist(LUSC_filtered$variant_freq,breaks=100)

# 筛选高频基因
high_freq_genesLUSC <- LUSC_filtered%>% 
  filter(variant_freq >= 0.70) %>% 
  arrange(desc(variant_freq))


title=tempdir()
results1 = ConsensusClusterPlus(as.matrix(high_freq_genesLUSC), 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "LUSC_CNtop5000_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "euclidean", 
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")




#保存工作空间
setwd("/home/yanbin/NSCLC_subtyping/CN")
save.image("NSCLC_CN_Clustering.RData")

