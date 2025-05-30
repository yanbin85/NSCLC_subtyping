library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(circlize)
install.packages("circlize")
#读入数据与重命名
load("/home/data/yanbin/NSCLC_Dataset/LUAD/TME/IOBR.expr.tme_combine.TCGA_LUAD.Rdata")
luad.expr.cibersort <- expr.cibersort
luad.expr.cibersort_abs <- expr.cibersort_abs
luad.expr.epic <- expr.epic
luad.expr.estimate <-expr.estimate
luad.expr.ips <- expr.ips
luad.expr.mcp <- expr.mcp
luad.expr.quantiseq <- expr.quantiseq
luad.expr.timer <- expr.timer
luad.expr.tme_combine <- expr.tme_combine
load("/home/data/yanbin/NSCLC_Dataset/LUSC/TME/IOBR.expr.tme_combine.TCGA_LUSC.Rdata")
lusc.expr.cibersort <- expr.cibersort
lusc.expr.cibersort_abs <- expr.cibersort_abs
lusc.expr.epic <- expr.epic
lusc.expr.estimate <-expr.estimate
lusc.expr.ips <- expr.ips
lusc.expr.mcp <- expr.mcp
lusc.expr.quantiseq <- expr.quantiseq
lusc.expr.timer <- expr.timer
lusc.expr.tme_combine <- expr.tme_combine
setwd("~/NSCLC_subtyping/Unsupervised Clustering/TME")
#数据合并
luad_data <- t(luad.expr.tme_combine)
colnames(luad_data) <- luad_data[1,]
luad_data <- luad_data[-1,]
lusc_data <- t(lusc.expr.tme_combine)
colnames(lusc_data) <- lusc_data[1,]
lusc_data <- lusc_data[-1,]
NSCLC <- as.data.frame(cbind(luad_data,lusc_data))
#save(NSCLC,file="TME.rda")

luad_cibersort_data <- t(luad.expr.cibersort)
colnames(luad_cibersort_data) <- luad_cibersort_data[1,]
luad_cibersort_data <- luad_cibersort_data[-1,]
lusc_cibersort_data <- t(lusc.expr.cibersort)
colnames(lusc_cibersort_data) <- lusc_cibersort_data[1,]
lusc_cibersort_data <- lusc_cibersort_data[-1,]
NSCLC_cibersort <- as.data.frame(cbind(luad_cibersort_data,lusc_cibersort_data))

# 绘制热图 --------------------------------------------------------------------
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
sample_codes <- substr(colnames(NSCLC_cibersort), 6, 7)

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
  annotation_name_side = "left",
  show_legend = FALSE
)

#删除最后三行
NSCLC_cibersort <- head(NSCLC_cibersort, -3)
# 将非数值列转换为数值型
NSCLC_cibersort <- NSCLC_cibersort %>%
  mutate(across(
    where(~ !is.numeric(.)), as.numeric))  # 尝试转换所有列为数值型
#绘图并保存
#pdf("TME_cibersort_heatmap.pdf",width = 12,height = 9)
png("TME_cibersort_heatmap.png",
    width = 12,          # 宽度（英寸）保持原PDF尺寸
    height = 9,          # 高度（英寸）
    units = "in",        # 单位使用英寸
    res = 600,           # 分辨率(600dpi满足打印级高清需求)
    type = "cairo")      # 使用抗锯齿渲染
Heatmap(as.matrix(NSCLC_cibersort),
        name = "TME cibersort",
        clustering_distance_rows = "pearson",
        cluster_columns = FALSE,  # 禁用列聚类以保持分组顺序
        top_annotation = column_ha,
        row_dend_gp = gpar(lwd = 0.3), #调整树状图线条粗细
        column_split = sample_groups,  # 按组分割列
        column_title = c("LUAD", "LUSC"),  # 分组标题
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 16),
          labels_gp = gpar(fontsize = 14),
          direction = "vertical"
        ),
        col = colorRamp2(c(0,0.05, 1), c("blue","white" , "red")))  # 颜色映射
dev.off()


# ConsensusClusterPlus聚类 --------------------------------------------------
setwd("~/NSCLC_subtyping/Unsupervised Clustering/TME/Cluster_results")
#聚类
title=tempdir()
inputdata <- as.matrix(NSCLC)
#转换为数值型
inputdata_numeric <- matrix(as.numeric(as.matrix(inputdata)), 
                         nrow = nrow(inputdata),
                         dimnames = list(rownames(inputdata), colnames(inputdata)))
#检查是否为数值型
print(sapply(inputdata_numeric, class))  
#执行聚类
title=tempdir()
results = ConsensusClusterPlus(inputdata_numeric, 
                               maxK = 10, 
                               reps = 1000, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "NSCLC_TME_ConsensusCluster", 
                               clusterAlg = "hc", 
                               innerLinkage="ward.D2",
                               finalLinkage="ward.D2",
                               distance = "pearson", #适用于连续型评分数据
                               seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]

#write.csv(results[[4]][["consensusClass"]],"TME.cluster.results.csv")
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

optK = Kvec[which.min(PAC)]  # 理想的K值为2

# 保存工作空间 ------------------------------------------------------------------
setwd("~/NSCLC_subtyping/Unsupervised Clustering/TME")
save.image("NSCLC_TME_Clustering.RData")
# 
setwd("~/NSCLC_subtyping/Unsupervised Clustering/TME/test")
聚类
 title=tempdir()
 inputdata2 <- as.matrix(NSCLC)
 results = ConsensusClusterPlus(inputdata2, 
                                maxK = 10, 
                                reps = 1000, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                title = "NSCLC_TME_ConsensusCluster", 
                                innerLinkage="ward.D2",
                                clusterAlg = "hc", 
                                distance = "euclidean", #适用于连续型评分数据
                                seed = 1262118388.71279, 
                                tmyPal=NULL, writeTable=TRUE,
                                plot = "pdf")
