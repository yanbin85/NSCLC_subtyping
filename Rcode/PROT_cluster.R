library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(tibble)
#LUAD与LUSC数据读取与合并
PROT_LUAD <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/PROT/TCGA-LUAD.protein.tsv", header = TRUE, sep = "\t",row.names=1)
PROT_LUSC <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/PROT/TCGA-LUSC.protein.tsv", header = TRUE, sep = "\t",row.names=1)
PROT <- cbind(PROT_LUAD,PROT_LUSC)
save(PROT,file="PROT.rda")
#移除含有na值的蛋白
PROT_filtered <- na.omit(PROT)

#原数据已进行了normalized
# #z-score对每一行进行归一化
# row_means <- apply(PROT_filtered, 1, mean, na.rm = TRUE)
# row_sds <- apply(PROT_filtered, 1, sd, na.rm = TRUE)
# #处理标准差为0的行（避免除以0）
# row_sds[row_sds == 0] <- 1  # 若标准差为0，设为1（此时Z-Score为0）
# #计算Z-Score
# PROT_zscore <- sweep(PROT_filtered, 1, row_means, "-")  # 减去均值
# PROT_zscore <- sweep(PROT_zscore, 1, row_sds, "/")    # 除以标准差

# 绘制热图 --------------------------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/PROT")


# mads <-apply(PROT_filtered, 1, mad) 
# hist(mads,breaks=100)
# heatmap_input<- PROT_zscore[rev(order(mads)), ] #按sd值排序


#分组
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
# 提取样本条形码第5-6位代码
sample_codes <- substr(colnames(heatmap_input), 6, 7)

# 创建分组注释向量（LUAD=蓝色，LUSC=红色，其他样本=灰色）
sample_groups <- case_when(
  sample_codes %in% luad_codes ~ "LUAD",
  sample_codes %in% lusc_codes ~ "LUSC",
)
sample_groups <- factor(sample_groups, levels = c("LUAD", "LUSC"))
#定义分组颜色
group_colors <- c("LUAD" = "#377EB8", "LUSC" = "#E41A1C")
library(circlize)
col_fun <- colorRamp2(
  breaks = c(-5, -1, 0, 1, 5, 10),  # 断点调整（关键阈值对齐）
  colors = c( "#0000FF",  "#6495ED","#FFFFFF",  "#FF6666","#FF0000","#8B0000" )
)
#添加注释并分割列
column_ha <- HeatmapAnnotation(
  group = sample_groups,
  col = list(group = group_colors),
  annotation_name_side = "left"
)
#绘图
pdf("PROT_heatmap.pdf",width = 12,height = 9)
Heatmap(
  as.matrix(PROT_filtered),
  name = "protein expression",
  col = col_fun,
  top_annotation = column_ha,
  row_dend_gp = gpar(lwd = 0.3), #调整树状图线条粗细
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
  show_row_names = FALSE,
  show_column_names = FALSE,
)
dev.off()


# ConsensusCluserPlus聚类 ---------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/PROT/Cluster_results")
input_data <- as.matrix(PROT_filtered)
title=tempdir()
result <- ConsensusClusterPlus(
 input_data,
  maxK = 10,
  reps = 1000,
  pItem = 0.8,
  pFeature = 0.8,  # 可调整基因抽样比例
  clusterAlg = "hc",
 innerLinkage="ward.D2",
  distance = "spearman",
  seed = 1262118388.71279,
  writeTable=TRUE,
  title ="NSCLC_PROT_ConsensusCluster",
  plot = "pdf"
)

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(result, title = "consensus_cluster", plot = "pdf")
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
  M = result[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri( M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK = Kvec[which.min(PAC)]  # 理想的K值为10

setwd("/home/yanbin/NSCLC_subtyping/PROT")
save.image("NSCLC_PROT_Clustering.RData")

# 
# > setwd("/home/yanbin/NSCLC_subtyping/PROT/test")
# > input_data <- t(as.matrix(PROT_filtered))
# > title=tempdir()
# > result <- ConsensusClusterPlus(
#   +     input_data,
#   +     maxK = 10,
#   +     reps = 1000,
#   +     pItem = 0.8,
#   +     pFeature = 0.8,  # 可调整基因抽样比例
#   +     clusterAlg = "pam",
#   +     distance = "spearman",
#   +     seed = 1262118388.71279,
#   +     writeTable=TRUE,
#   +     title ="NSCLC_PROT_ConsensusCluster",
#   +     plot = "pdf"
#   + )
