# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ConsensusClusterPlus", "genefilter"), force = TRUE)

library(ConsensusClusterPlus)
library(genefilter)
library(ComplexHeatmap)
library(tibble)
#LUAD与LUSC数据读取与合并
GEXP <- read.table("/home/data/yanbin/NSCLC_Dataset/LUNG_GEXP/TCGA.LUNG.sampleMap_HiSeqV2_PANCAN", header = TRUE, sep = "\t",row.names=1)
#save(GEXP,file="GEXP.rda")
#原数据已经过pan-cancer normalized log2(norm_count+1)转换

#z-score再次对每一行进行归一化
row_means <- apply(GEXP, 1, mean, na.rm = TRUE)
row_sds <- apply(GEXP, 1, sd, na.rm = TRUE)
#处理标准差为0的行（避免除以0）
row_sds[row_sds == 0] <- 1  # 若标准差为0，设为1（此时Z-Score为0）
#计算Z-Score
GEXP_zscore <- sweep(GEXP, 1, row_means, "-")  # 减去均值
GEXP_zscore <- sweep(GEXP_zscore, 1, row_sds, "/")    # 除以标准差
# 绘制热图 --------------------------------------------------------------------
mads <-apply(GEXP_zscore, 1, mad) # 计算每一行的MAD值
heatmap_input<- GEXP_zscore[rev(order(mads))[1:1000], ] # 提取前1000个基因

# 转换基因标识符
gene_mapping <- read_tsv("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fhugo_gencode_good_hg19_V24lift37_probemap")
# gene_mapping2 <- data.frame(id=gene_mapping$id,gene=gene_mapping$gene)
# heatmap_input<- as.data.frame(heatmap_input) %>% 
#   tibble::rownames_to_column(var = "ensembl_id") %>%
#   left_join(gene_mapping2, by=c("ensembl_id"="id")) %>%
#   select(gene, everything()) 
#去除重复的gene并将其设为行名
# heatmap_input <- heatmap_input %>%
#   distinct(gene, .keep_all = TRUE) %>%  # 保留每个 gene 的第一行
#   column_to_rownames(var = "gene")
# heatmap_input <- heatmap_input[,-1]

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

# 创建分组注释向量
sample_groups <- case_when(
  sample_codes %in% luad_codes ~ "LUAD",
  sample_codes %in% lusc_codes ~ "LUSC",
)
sample_groups <- factor(sample_groups, levels = c("LUAD", "LUSC"))

# summary(sample_groups)
# LUAD LUSC 
# 576  553 
#定义分组颜色
group_colors <- c("LUAD" = "#377EB8", "LUSC" = "#E41A1C")
# 定义更鲜艳的颜色映射（纯红/蓝 + 减少中间过渡色）
col_fun <- colorRamp2(
  breaks = c(-5, -1, 0, 1, 5, 15),  # 断点调整（关键阈值对齐）
  colors = c( "#0000FF",  "#6495ED","#FFFFFF",  "#FF6666","#FF0000","#8B0000" )
)
# 添加注释并分割列
column_ha <- HeatmapAnnotation(
  group = sample_groups,
  col = list(group = group_colors),
  annotation_name_side = "left"
)
#绘图
setwd("/home/yanbin/NSCLC_subtyping/GEXP/COUNTS_new")
pdf("GEXP2_top1000_heatmap.pdf",width = 12,height = 9)
Heatmap(
  as.matrix(heatmap_input),
  name = "gene expression",
  col = col_fun,
  top_annotation = column_ha,
  row_dend_gp = gpar(lwd = 0.3), #调整树状图线条粗细
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  show_column_names  = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
)
dev.off()



#########NSCLC前20%可变性最大的基因聚类#########
setwd("/home/yanbin/NSCLC_subtyping/GEXP/COUNTS_new/cluster_results")
#为了选择信息最丰富的基因进行类检测，将数据集减少到前20%可变性最大的基因（通过中值绝对偏差 - MAD来衡量）
cluster_input <- GEXP_zscore[rev(order(mads))[1:(0.2*nrow(GEXP_zscore))], ] # 提取前20%基因（4106个）
#聚类
title=tempdir()
results1 = ConsensusClusterPlus(as.matrix(cluster_input), 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_GEXPtop0.2_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "pearson",
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

# #根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
# Kvec = 2:20
# x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
# PAC = rep(NA,length(Kvec)) 
# names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
# for(i in Kvec){
#   M = results1[[i]]$consensusMatrix
#   Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
#   PAC[i-1] = Fn(x2) - Fn(x1)
# } 
# 
# optK1 = Kvec[which.min(PAC)]  # 理想的K值为10




#保存工作空间####
setwd("/home/yanbin/NSCLC_subtyping/GEXP/COUNTS_new")
save.image("NSCLC_GEXP2_Clustering.RData")