library(data.table)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
setwd("/home/yanbin/NSCLC_subtyping/MUTA/")
#读取数据
MUTA_Pan <- fread(
  "/home/data/yanbin/NSCLC_Dataset/MUTA/mc3.v0.2.8.PUBLIC.nonsilentGene.xena",
  header = TRUE,          # 第一行是样本 ID
  sep = "\t",            
  check.names = FALSE,    # 禁止自动修改列名
)
setDF(MUTA_Pan, rownames = MUTA_Pan[[1]])  # 设置行名为第一列
MUTA_Pan <- MUTA_Pan[, -1, drop = FALSE]   # 移除第一列


# 提取NSCLC样本数据 -------------------------------------------------------------
# 定义非小细胞肺癌的项目代码列表
nsclc_codes <- c(
  "05", "17", "35", "38", "44", "49", "4B", "50", "53", "55", 
  "62", "64", "67", "69", "71", "73", "75", "78", "80", "83", 
  "86", "91", "93", "95", "97", "99", "J2", "L4", "L9", "ME", 
  "MN", "MP", "NB", "NJ", "O1", "S2", "T6", "11", "18", "21", 
  "22", "33", "34", "37", "39", "43", "46", "51", "52", "56", 
  "58", "60", "63", "66", "68", "6A", "70", "77", "79", "82", 
  "85", "90", "92", "94", "96", "98", "J1", "L3", "LA", "MF", 
  "ML", "NC", "NK", "O2", "OC", "UJ", "WG", "XC", "ZE"
)
# 提取所有样本 ID
all_samples <- colnames(MUTA_Pan)
# 构建正则表达式
regex_pattern <- paste0("^TCGA-(", paste(nsclc_codes, collapse = "|"), ")-.*")
# 筛选 NSCLC 样本
nsclc_samples <- grep(regex_pattern, colnames(MUTA_Pan), value = TRUE, ignore.case = TRUE)
# 提取子集
MUTA_NSCLC <- MUTA_Pan[, nsclc_samples]
#save(MUTA_NSCLC,file="MUTA_NSCLC.rda")


# 数据过滤 --------------------------------------------------------------------
#移除在所有样本中完全一致（全0或全1）的基因
MUTA_filtered <- MUTA_NSCLC[apply(MUTA_NSCLC, 1, function(x) sd(x) != 0), ]
#保留至少在5%以上样本中突变的基因（651个）
MUTA_filtered <- MUTA_NSCLC[rowSums(MUTA_NSCLC != 0) >=  0.05 * ncol(MUTA_NSCLC), ]


# 绘制热图 --------------------------------------------------------------------

library(ComplexHeatmap)
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
sample_codes <- substr(colnames(MUTA_filtered), 6, 7)

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
# 定义颜色映射（0=白色，1=红色）
col_fun <- colorRamp2(breaks = c(0, 1), colors = c("white", "red"))



# 绘制热图
pdf("MUTA_heatmap.pdf", width = 12, height = 9)  # 调整宽高
Heatmap(
  as.matrix(MUTA_filtered),
  name = "non-silent mutation",  # 图例标题
  col = col_fun,      # 颜色映射
  top_annotation = column_ha,
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  row_dend_gp = gpar(lwd = 0.3),
  cluster_columns = FALSE,  # 禁用列聚类（保持样本原始顺序）
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_side = "right",  
  heatmap_legend_param = list(
    labels = c("No Mutation", "Mutation"),  # 图例标签
    at = c(0, 1)  # 对应颜色断点
  )
)
dev.off()

#热图2
# 计算突变频率
total_samples <- ncol(MUTA_filtered) 
MUTA_sort <- MUTA_filtered
MUTA_sort$variant_count <- apply(MUTA_sort[, ], 1, function(x) sum(x == 1))
MUTA_sort$variant_freq <- MUTA_sort$variant_count / total_samples
hist(MUTA_sort$variant_freq,breaks=100)
hist(MUTA_sort$variant_freq)

# 按variant_freq降序排序，并提取前100行
heatmap_input <- MUTA_sort[order(MUTA_sort$variant_freq, decreasing = TRUE), ][1:100, ]

# 提取原始样本数据列（排除variant_count和variant_freq）
heatmap_input <- heatmap_input[, colnames(MUTA_filtered)]


pdf("MUTA_top100_heatmap.pdf", width = 12, height = 9)  
Heatmap(
  as.matrix(heatmap_input),
  name = "non-silent mutation",  # 图例标题
  col = col_fun,      # 颜色映射
  top_annotation = column_ha,
  column_split = sample_groups,  # 按组分割列
  column_title = c("LUAD", "LUSC"),  # 分组标题
  row_dend_gp = gpar(lwd = 0.3),
  cluster_columns = FALSE,  # 禁用列聚类（保持样本原始顺序）
  show_column_names  = FALSE, #不显示列名
  row_names_side = "right",  
  row_names_gp = gpar(fontsize = 6),  
  heatmap_legend_param = list(
    labels = c("No Mutation", "Mutation"),  # 图例标签
    at = c(0, 1)  # 对应颜色断点
  )
)
dev.off()

# # 绘制瀑布图 -------------------------------------------------------------------
# library(maftools)  
# library(dplyr)     
# gene_mapping <- read_tsv("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fhugo_gencode_good_hg19_V24lift37_probemap")
# # 转换为长格式（maftools要求）
# mutation_long <- MUTA_NSCLC %>%
#   tibble::rownames_to_column("Gene") %>%
#   tidyr::gather(key = "Sample", value = "Mutation", -Gene) %>%
#   filter(Mutation == 1)  # 只保留突变记录（Mutation=1）
# # 修正后的长格式数据框
# mutation_long <- mutation_long%>%
#   left_join(
#     gene_mapping %>% 
#       select(gene, chrom, chromStart, chromEnd),
#     by = c("Gene"="gene")
#   )
# maf_input <- mutation_long %>%
#   rename(
#     Hugo_Symbol = Gene,                   # 必须字段：基因名
#     Tumor_Sample_Barcode = Sample,         # 必须字段：样本ID
#     Variant_Classification = Mutation ,     # 必须字段：突变类型
#     Chromosome = chrom,                   # 染色体
#     Start_Position = chromStart,          # 突变起始位置
#     End_Position = chromEnd,               # 突变终止位置
#   ) %>%
#   mutate(
#     Variant_Classification = ifelse(Variant_Classification == 1, "Gene level non-silent mutation", "no mutation"),  # 分类标签
#   )
# maf_input$Reference_Allele = NA
# maf_input$Tumor_Seq_Allele2 = NA
# maf_input$Variant_Type = NA
# # 仅保留有突变的记录
# maf_input <- maf_input %>% filter(Variant_Classification == "Gene level non-silent mutation")
# # 创建 MAF 对象
# maf_data <- read.maf(
#   maf = maf_input,
#   clinicalData = data.frame(
#     Tumor_Sample_Barcode = colnames(MUTA_NSCLC),  # 样本ID需与MAF中的列名一致
#     Group = sample_groups  # 假设已定义分组变量（如LUAD/LUSC）
#   )
# )
# 
# 
# 
# 

# ConsensusClusterPlus聚类 --------------------------------------------------
setwd("/home/yanbin/NSCLC_subtyping/MUTA/Cluster_results/")
title=tempdir()
results <- ConsensusClusterPlus(
  d = as.matrix(MUTA_filtered),
  maxK = 10,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = "NSCLC_Mutation_ConsensusCluster",
  clusterAlg = "hc",
  distance = "binary",  #使用Jaccard距离适合二进制数据
  innerLinkage="ward.D2",
  finalLinkage="ward.D2",
  seed = 1262118388,
  writeTable = TRUE,
  plot = "pdf"
)

#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl = calcICL(results, title = "consensus_cluster", plot = "pdf")
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]][1:5,]

dim(icl[["itemConsensus"]])
icl[["itemConsensus"]][1:5,]

#保存工作空间####
setwd("/home/yanbin/NSCLC_subtyping/MUTA")
save.image("NSCLC_MUTA_Clustering.RData")
