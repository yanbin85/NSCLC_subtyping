library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(tibble)
library(data.table)
library(readr)
setwd("/home/yanbin/NSCLC_subtyping/DNA_MET")
#LUAD与LUSC数据读取与合并
#MET_LUAD <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/DNA_MET/TCGA-LUAD.methylation450.tsv", header = TRUE, sep = "\t",row.names=1)
#MET_LUSC <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/DNA_MET/TCGA-LUSC.methylation450.tsv", header = TRUE, sep = "\t",row.names=1)


# 读取LUAD数据
MET_LUAD <- fread("/home/data/yanbin/NSCLC_Dataset/LUAD/DNA_MET/TCGA-LUAD.methylation450.tsv", sep = "\t", header = TRUE)
MET_LUAD <- as.data.frame(MET_LUAD)   # 转换为data.frame（可选）
rownames(MET_LUAD) <- MET_LUAD$`Composite Element REF`# 设置行名
MET_LUAD <- MET_LUAD[, -1]          
# 读取LUSC数据
MET_LUSC <- fread("/home/data/yanbin/NSCLC_Dataset/LUSC/DNA_MET/TCGA-LUSC.methylation450.tsv", sep = "\t", header = TRUE)
MET_LUSC <- as.data.frame(MET_LUSC)
rownames(MET_LUSC) <- MET_LUSC[[1]]
MET_LUSC <- MET_LUSC[, -1]
MET <- cbind(MET_LUAD,MET_LUSC)
write.csv(MET,"MET.csv")
#数据过滤
# 允许每个探针最多10%样本含NA
#max_na_samples <- ceiling(0.1 * ncol(MET))
#MET_filtered <- MET[rowSums(is.na(MET)) <= max_na_samples, ]
#去除含有NA值的位点
MET_filtered <- na.omit(MET)
#计算sd值
sds <- apply(MET_filtered, 1, sd, na.rm=TRUE)
hist(sds,breaks = 100)

# 绘制热图 --------------------------------------------------------------------
MET_highvar <- MET_filtered[rev(order(sds))[1:1000], ] # 提取前100个基因
heatmap_input <- MET_highvar%>% 
  tibble::rownames_to_column(var = "id") %>%
  left_join(gene_mapping2, by=c("id"="id")) %>%
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
sample_codes <- substr(colnames(MET_highvar), 6, 7)

# 创建分组注释向量（LUAD=蓝色，LUSC=红色，其他样本=灰色）
sample_groups <- case_when(
  sample_codes %in% luad_codes ~ "LUAD",
  sample_codes %in% lusc_codes ~ "LUSC",
)
sample_groups <- factor(sample_groups, levels = c("LUAD", "LUSC"))
#定义分组颜色
group_colors <- c("LUAD" = "#377EB8", "LUSC" = "#E41A1C")
# 添加注释并分割列
column_ha <- HeatmapAnnotation(
  group = sample_groups,
  col = list(group = group_colors),
  annotation_name_side = "left"
)
#绘图并保存
pdf("MET_top1000_heatmap.pdf",width = 12,height = 9)
Heatmap(as.matrix(MET_highvar),
        name = "Methylation",
        clustering_distance_rows = "pearson",
        cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
        top_annotation = column_ha,
        row_dend_gp = gpar(lwd = 0.3), #调整树状图线条粗细
        column_split = sample_groups,  # 按组分割列
        column_title = c("LUAD", "LUSC"),  # 分组标题
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        col = colorRamp2(c(0, 0.3, 1), c("blue", "white", "red")))  # 颜色映射
dev.off()




# ConsensuClusterPlus聚类 ---------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/DNA_MET/Cluster_results")
#选择前1%可变性最大的基因聚类(8064个)
cluster_input <- MET_filtered[rev(order(sds))[1:(0.02*nrow(MET_filtered))], ]
#删除含有NA值的基因（剩下4161个）
cluster_input <- na.omit(cluster_input)
title=tempdir()
results = ConsensusClusterPlus(as.matrix(cluster_input), 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_METtop0.01_ConsensusCluster", 
                                clusterAlg = "hc", 
                                distance = "pearson", 
                               innerLinkage="ward.D2",
                               finalLinkage="ward.D2",
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]

#根据PAC = Proportion of ambiguous clustering 模糊聚类比例确定最佳k值（有时候会失灵）
Kvec = 2:10
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK1 = Kvec[which.min(PAC)]  # 理想的K值为10



gene_mapping <- read_tsv("https://gdc-hub.s3.us-east-1.amazonaws.com/download/HM450.hg38.manifest.gencode.v36.probeMap")
gene_mapping2 <- data.frame(id=gene_mapping$`#id`,gene=gene_mapping$gene)



#保存工作空间
setwd("/home/yanbin/NSCLC_subtyping/DNA_MET")
save.image("NSCLC_MET_Clustering.RData")
