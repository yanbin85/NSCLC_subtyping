# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ConsensusClusterPlus", "genefilter"), force = TRUE)

library(ConsensusClusterPlus)
library(genefilter)
library(ComplexHeatmap)
library(tibble)
#LUAD与LUSC数据读取与合并
GEXP_LUAD <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/GEXP/tpm/TCGA-LUAD.star_tpm.tsv", header = TRUE, sep = "\t",row.names=1)
GEXP_LUSC <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/GEXP/tpm/TCGA-LUSC.star_tpm.tsv", header = TRUE, sep = "\t",row.names=1)
GEXP<- cbind(GEXP_LUAD,GEXP_LUSC)

#转换为矩阵
GEXP_matrix <- as.matrix(GEXP)

#z-score归一化(基因水平)
row_means <- apply(GEXP_matrix, 1, mean, na.rm = TRUE)
row_sds <- apply(GEXP_matrix, 1, sd, na.rm = TRUE)
#处理标准差为0的行（避免除以0）
row_sds[row_sds == 0] <- 1  # 若标准差为0，设为1（此时Z-Score为0）
#计算Z-Score
GEXP_zscore <- sweep(GEXP_matrix, 1, row_means, "-")  # 减去均值
GEXP_zscore <- sweep(GEXP_zscore, 1, row_sds, "/")    # 除以标准差

# 绘制热图 --------------------------------------------------------------------
mads <-apply(GEXP_zscore, 1, mad) # 计算每一行的MAD值
heatmap_input<- GEXP_zscore[rev(order(mads))[1:1000], ] # 提取前1000个基因

# 转换基因标识符
gene_mapping <- read_tsv("https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap")
gene_mapping2 <- data.frame(id=gene_mapping$id,gene=gene_mapping$gene)
heatmap_input<- as.data.frame(heatmap_input) %>% 
     tibble::rownames_to_column(var = "ensembl_id") %>%
     left_join(gene_mapping2, by=c("ensembl_id"="id")) %>%
     select(gene, everything()) 
#去除重复的gene并将其设为行名
heatmap_input <- heatmap_input %>%
     distinct(gene, .keep_all = TRUE) %>%  # 保留每个 gene 的第一行
     column_to_rownames(var = "gene")
heatmap_input <- heatmap_input[,-1]

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
  breaks = c(-4, -2, 0, 2, 5, 10),  # 多段断点
  colors = c("darkblue", "blue", "white", "#FF9999", "#FF6666", "#FF0000")
)
# 添加注释并分割列
column_ha <- HeatmapAnnotation(
     group = sample_groups,
     col = list(group = group_colors),
     annotation_name_side = "left"
     )
#绘图
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM")
pdf("GEXP_top1000_heatmap.pdf",width = 12,height = 9)
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

heatmap_input2<- GEXP_zscore[rev(order(mads))[1:100], ] # 提取前500个基因
heatmap_input2<- as.data.frame(heatmap_input2) %>% 
  tibble::rownames_to_column(var = "ensembl_id") %>%
  left_join(gene_mapping2, by=c("ensembl_id"="id")) %>%
  select(gene, everything()) 
#去除重复的gene并将其设为行名
heatmap_input2 <- heatmap_input2 %>%
  distinct(gene, .keep_all = TRUE) %>%  # 保留每个 gene 的第一行
  column_to_rownames(var = "gene")
heatmap_input2 <- heatmap_input2[,-1]

#绘图
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM")
pdf("GEXP_top1000_heatmap2.pdf",width = 12,height = 9)
Heatmap(
  as.matrix(heatmap_input2),
  name = "gene expression",
  col = col_fun,
  top_annotation = column_ha,
  row_dend_gp = gpar(lwd = 0.3), #调整树状图线条粗细
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  show_column_names  = FALSE,
  row_names_gp = gpar(fontsize = 5),
  cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
)
dev.off()
######NSCLC总基因聚类######
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM/results")
#ConsensusClusterPlus聚类
title=tempdir()
results = ConsensusClusterPlus(GEXP, 
                               maxK = 20, 
                               reps = 1000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "NSCLC_GEXP_ConsensusCluster", 
                               clusterAlg = "hc", 
                               distance = "pearson", 
                               seed = 1262118388.71279, 
                               tmyPal=NULL, 
                               writeTable=TRUE, 
                               plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]


#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:20
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK = Kvec[which.min(PAC)]  # 理想的K值为15


#########NSCLC前5%可变性最大的基因聚类#########
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM/results1")
#为了选择信息最丰富的基因进行类检测，将数据集减少到前5%可变性最大的基因（通过中值绝对偏差 - MAD来衡量）
mads <-apply(GEXP_zscore, 1, mad) # 计算每一行的MAD值
GEXP1 <- GEXP_zscore[rev(order(mads))[1:(0.05*nrow(GEXP_zscore))], ] # 提取前5%基因（3033个）
#聚类
title=tempdir()
results1 = ConsensusClusterPlus(as.matrix(GEXP1), 
                                maxK = 20, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_GEXPtop0.05_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "pearson", 
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
Kvec = 2:20
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results1[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK1 = Kvec[which.min(PAC)]  # 理想的K值为20

#########LUAD聚类######
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM/LUAD_results")
#归一化
GEXP_LUAD1<-sweep(GEXP_LUAD, 1, apply(GEXP_LUAD, 1, median, na.rm = T))
GEXP_LUAD1<-as.matrix(GEXP_LUAD1)
title=tempdir()
results_LUAD = ConsensusClusterPlus(GEXP_LUAD1, 
                                    maxK = 20, 
                                    reps = 1000, 
                                    pItem = 0.8, 
                                    pFeature = 1, 
                                    title = "LUAD_GEXP_ConsensusCluster", 
                                    clusterAlg = "hc", 
                                    distance = "pearson", 
                                    seed = 1262118388.71279, 
                                    tmyPal=NULL, 
                                    writeTable=TRUE, 
                                    plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results_LUAD, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]


#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:20
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results_LUAD[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK_LUAD = Kvec[which.min(PAC)]  # 理想的K值为20

#########PCA降维后的50%NSCLC高变基因聚类######
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM/results2")
#根据mad值选取前50%高变异基因
mads <-apply(GEXP, 1, mad) # 计算每一行的MAD值
GEXP_filtered <- GEXP[rev(order(mads))[1:30330], ] # 提取前50%的基因
rownames(GEXP_filtered)<-substring(rownames(GEXP_filtered), 1, 15)
#PCA降维
library(ggplot2)
#数据标准化（Z-score标准化，确保每个基因均值为0，标准差为1）
GEXP_scaled <- t(scale(t(GEXP_filtered)))

# run PCA
pca_result <- prcomp(t(GEXP_scaled), scale. = FALSE)
# 确定主成分数（累计解释>80%方差）
cum_var <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_pcs <- which(cum_var >= 0.8)[1]
#计算基因重要性评分（方差加权）
loadings <- pca_result$rotation[, 1:n_pcs]
variance_weights <- pca_result$sdev[1:n_pcs]^2 / sum(pca_result$sdev[1:n_pcs]^2)
gene_scores <- apply(loadings, 1, function(x) sum(x^2 * variance_weights))
#筛选Top驱动基因（取前2000）
top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:2000]
GEXP_key <- GEXP_scaled[top_genes, ]


#执行聚类（使用Pearson距离层次聚类)
title=tempdir()
results2 <- ConsensusClusterPlus(
  d = t(GEXP_key),
  maxK = 20,
  reps = 1000,
  pItem = 0.8,
  pFeature = 0.8,  # 可调整基因抽样比例
  clusterAlg = "hc",
  distance = "pearson",
  seed = 1262118388.71279,
  writeTable=TRUE,
  title ="NSCLC_GEXPhalf_ConsensusCluster",
  plot = "pdf"
)
#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results2, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]

#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:20
x1 = 0.1; x2 = 0.9        
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  
for(i in Kvec){
  M = results2[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri( M)])          # M 为计算出的共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK2 = Kvec[which.min(PAC)]  # 理想的K值为20
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM")


#保存工作空间####
setwd("/home/yanbin/NSCLC_subtyping/GEXP/TPM")
save.image("NSCLC_GEXP_Clustering.RData")