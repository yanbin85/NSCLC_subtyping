library(readr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(ggplot2)
library(tidyverse)
setwd("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/")

# 数据处理 --------------------------------------------------------------------
#GEXP4
GEXP.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/GEXP.cluster.results.csv", header = TRUE, sep = ",",row.names=1)
#MET4
MET.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/MET.cluster.results.csv", header = TRUE, sep = ",",row.names=1)
#CN5
CN.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/CN.cluster.results.csv", header = TRUE, sep = ",",row.names=1)
#MUTA4
MUTA.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/MUTA.cluster.results.csv", header = TRUE, sep = ",",row.names=1)
#PROT4
PROT.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/PROT.cluster.results.csv", header = TRUE, sep = ",",row.names=1)
#TME4
TME.cluster.res <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/TME.cluster.results.csv",header = TRUE, sep = ",",row.names=1)
#Multiomics 5
load("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/cmoic.cluster.results.rda")
cmoic.cluster.res <- cmoic.nsclc$clust.res
cmoic.cluster.res <- cmoic.cluster.res[,-1,drop=FALSE]
colnames(cmoic.cluster.res) <- c("MOVICS clusters")

#6个组学共同样本
load("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/common_samples.rda")
rownames(GEXP.cluster.res) <- substr(rownames(GEXP.cluster.res), start = 1, stop = 15)
rownames(GEXP.cluster.res) <- gsub("\\.", "-", rownames(GEXP.cluster.res))
colnames(GEXP.cluster.res) <- c("GEXP clusters")
GEXP.cluster.res <- GEXP.cluster.res[common_samples, , drop = FALSE]


MET.cluster.res <- as.data.frame(t(MET.cluster.res))
colnames(MET.cluster.res) <- substr(colnames(MET.cluster.res), start = 1, stop = 15)
colnames(MET.cluster.res) <- gsub("\\.", "-", colnames(MET.cluster.res))
rownames(MET.cluster.res) <- c("MET clusters")
MET.cluster.res <- MET.cluster.res[, common_samples, drop = FALSE]
MET.cluster.res <- as.data.frame(t(MET.cluster.res))


rownames(CN.cluster.res) <- substr(rownames(CN.cluster.res), start = 1, stop = 15)
rownames(CN.cluster.res) <- gsub("\\.", "-", rownames(CN.cluster.res))
colnames(CN.cluster.res) <- c("CNV clusters")
CN.cluster.res <- CN.cluster.res[common_samples, , drop = FALSE]

colnames(MUTA.cluster.res) <- c("MUT clusters")
MUTA.cluster.res <- MUTA.cluster.res[common_samples, , drop = FALSE]

rownames(PROT.cluster.res) <- substr(rownames(PROT.cluster.res), start = 1, stop = 15)
rownames(PROT.cluster.res) <- gsub("\\.", "-", rownames(PROT.cluster.res))
colnames(PROT.cluster.res) <- c("PROT clusters")
PROT.cluster.res <- PROT.cluster.res[common_samples, , drop = FALSE]

TME.cluster.res <- as.data.frame(t(TME.cluster.res))
colnames(TME.cluster.res) <- substr(colnames(TME.cluster.res), start = 1, stop = 15)
colnames(TME.cluster.res) <- gsub("\\.", "-", colnames(TME.cluster.res))
rownames(TME.cluster.res) <- c("TME clusters")
TME.cluster.res <- TME.cluster.res[, common_samples, drop = FALSE]
TME.cluster.res <- as.data.frame(t(TME.cluster.res))

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

commonsample.clinical <- clinical.data[,c("disease_type","gender.demographic","project_id.project","name.project","ajcc_pathologic_stage.diagnoses","cigarettes_per_day.exposures","years_smoked.exposures","vital_status.demographic","age_at_index.demographic")]
#筛选551个样本的数据
commonsample.clinical <- as.data.frame(t(commonsample.clinical))
colnames(commonsample.clinical) <- substr(colnames(commonsample.clinical), start = 1, stop = 15)
commonsample.clinical <- commonsample.clinical[, common_samples, drop = FALSE]
commonsample.clinical <- as.data.frame(t(commonsample.clinical))
setwd("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/outputdata")
save(commonsample.clinical,file="commonsample.clinical.rda")

#生存数据
pancancer.survival <- read.table("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/inputdata/Survival_SupplementalTable_S1_20171025_xena_sp", header = TRUE,sep = "\t",row.names = 1,)
commonsample.survival <- pancancer.survival[common_samples,]
save(commonsample.survival,file="commonsample.survival.rda")
# luad_survial <- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/PHENO/survival/TCGA-LUAD.survival.tsv",
#                            header = TRUE,
#                            sep = "\t",
#                            row.names = 1,
#                            quote = "",    # 禁用引号解析
#                            fill = TRUE,   # 允许填充缺失列
#                            check.names = FALSE)  # 保留列名中的特殊字符
# lusc_survial <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/PHENO/survival/TCGA-LUSC.survival.tsv",
#                            header = TRUE,
#                            sep = "\t",
#                            row.names = 1,
#                            quote = "",    # 禁用引号解析
# #                            fill = TRUE,   # 允许填充缺失列
# #                            check.names = FALSE)  # 保留列名中的特殊字符
# #合并数据
# survial.data <- rbind(luad_survial,lusc_survial)
# #筛选551个样本的数据
# commonsample.survival <- survial.data
# commonsample.survival <- as.data.frame(t(commonsample.survival))
# colnames(commonsample.survival) <- substr(colnames(commonsample.survival), start = 1, stop = 15)
# valid_samples <- intersect(
#   colnames(commonsample.survival),  
#   common_samples              
# )
# if (length(valid_samples) > 0) {
#   commonsample.survival <- commonsample.survival[, valid_samples, drop=FALSE]
#   commonsample.survival <- as.data.frame(t(commonsample.survival))
# } else {
#   stop("错误：无共同样本。检查样本名格式是否一致。")
# }
# 
# 
# #合并临床与生存数据
# commonsample.clinical <- rownames_to_column(commonsample.clinical, var = "SampleID")
# commonsample.survival <- rownames_to_column(commonsample.survival, var = "SampleID")
# 
# commonsamples.clinical.merge <- merge(
#   x = commonsample.clinical,
#   y = commonsample.survival,
#   by="SampleID", # 按行名（样本ID）合并
#   all.x = TRUE        # 保留所有临床数据样本（左连接）
# )
# rownames(commonsamples.clinical.merge) <- commonsamples.clinical.merge$SampleID  # 恢复行名
# commonsamples.clinical.merge$SampleID <- NULL                   # 移除合并产生的冗余列
# commonsamples.clinical.merge <- commonsamples.clinical.merge[common_samples,]
# all(rownames(GEXP.cluster.res) == rownames(commonsamples.clinical.merge))  #检查行顺序
# #保存
# #save(commonsamples.clinical.merge,file="commonsamples.clinical.merge.rda")



# 7种分型与临床特征的关联矩阵热图 --------------------------------------------------------
#合并7种聚类结果与选定的临床数据列
intergration.cluster <- cbind(cmoic.cluster.res,GEXP.cluster.res)
intergration.cluster <- cbind(intergration.cluster,MUTA.cluster.res)
intergration.cluster <- cbind(intergration.cluster,CN.cluster.res)
intergration.cluster <- cbind(intergration.cluster,MET.cluster.res)
intergration.cluster <- cbind(intergration.cluster,PROT.cluster.res)
intergration.cluster <- cbind(intergration.cluster,TME.cluster.res)
intergration.cluster$gender <- commonsample.clinical$gender.demographic
intergration.cluster$project_id <- commonsample.clinical$project_id.project
intergration.cluster$pathologic_stage <- commonsample.clinical$ajcc_pathologic_stage.diagnoses

#stage_counts <- table(intergration.cluster$stage)
#  stage_counts 
# 
# Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III Stage IIIA Stage IIIB 
# 4          4        125        147          4         72         87          2         79         12 
# Stage IV 
# 15


#stage数据统一格式
intergration.cluster$pathologic_stage <- gsub("Stage IA", "Stage I", intergration.cluster$pathologic_stage)
intergration.cluster$pathologic_stage <- gsub("Stage IB", "Stage I", intergration.cluster$pathologic_stage)
intergration.cluster$pathologic_stage <- gsub("Stage IIA", "Stage II", intergration.cluster$pathologic_stage)
intergration.cluster$pathologic_stage <- gsub("Stage IIB", "Stage II", intergration.cluster$pathologic_stage)
intergration.cluster$pathologic_stage <- gsub("Stage IIIA", "Stage III", intergration.cluster$pathologic_stage)
intergration.cluster$pathologic_stage <- gsub("Stage IIIB", "Stage III", intergration.cluster$pathologic_stage)
intergration.cluster[] <- lapply(intergration.cluster, function(col) {
  col[trimws(col) == ""] <- "unknown"
  return(col)
})
#write.csv(intergration.cluster,file="intergration.cluster.csv")
#对样本排序
intergration.cluster <- intergration.cluster[order(intergration.cluster$gender), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$project_id), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$pathologic_stage), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`MOVICS clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`GEXP clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`MUT clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`CNV clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`MET clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`PROT clusters`), ]
intergration.cluster <- intergration.cluster[order(intergration.cluster$`TME clusters`), ]


#设置颜色
cluster_colors <- list(
  "MOVICS clusters" = c("1" = "#2EC4B6", "2" = "#E71D36", "3" = "#FF9F1C", "4" =  "#BDD5EA","5"="#FFA5BA"),
  "GEXP clusters" = c("1" = "#1F77B4", "2" = "#FF7F0E", "3" = "#2CA02C", "4" = "#D62728"),
  "MUT clusters" = c("1" = "#BCBD22", "2" = "#E377C2", "3" = "#7F7F7F", "4" = "#8C564B" ),
  "CNV clusters" = c("1" = "#9467BD", "2" = "#C5B0D5", "3" = "#17BECF", "4" = "#FF9896", "5" = "#AEC7E8"),
  "TME clusters" = c("1" = "#FFBB78", "2" = "#98DF8A", "3" = "#FF9896", "4" = "#C5B0D5"),
  "MET clusters" = c("1" = "#F7B6D2", "2" = "#C7C7C7", "3" = "#DBDB8D", "4" = "#9EDAE5"),
   "PROT clusters"= c("1" = "#393B79", "2" = "#DCDCDC", "3" = "#B5CF6B", "4" = "#8C6D31"),
  gender =c("male" = "#2CA02C", "female" = "#D62728"),
  project_id=c("TCGA-LUAD" =  "#D4AF37", "TCGA-LUSC" = "#6C5B7B"),
  pathologic_stage = c("Stage I"   = "#8DA0CB", 
                       "Stage II"  = "#E78AC3",
                       "Stage III" = "#FC8D62",  #
                       "Stage IV"  = "#66C2A5" ,
                       "unknown" = "white")
 )

# 转换为因子（确保离散值映射）
intergration.cluster[] <- lapply(intergration.cluster, factor)

# 创建注释对象
anno_list <- lapply(c( "MOVICS clusters","GEXP clusters", "MUT clusters", "CNV clusters", "MET clusters", "PROT clusters","TME clusters","gender", "project_id","pathologic_stage"), function(col_name) {
  HeatmapAnnotation(
    # 每个组学一列注释条
    name = col_name,
    df = intergration.cluster[col_name],
    col = setNames(list(cluster_colors[[col_name]]), col_name),
    show_legend = FALSE,  # 后续统一添加图例
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 12),
    simple_anno_size = unit(1, "cm")  # 条带高度
  )
})
# 组合所有注释条
ht <- Reduce(`%v%`, anno_list)   # 使用 %v% 操作符垂直拼接


# 添加图例
legend_order <- c( "MOVICS clusters","GEXP clusters", "MUT clusters", "CNV clusters", "MET clusters", "PROT clusters","TME clusters","gender", "project_id","pathologic_stage")
legend_list <- lapply(legend_order, function(omics) {
  Legend(
    title = omics,
    labels = names(cluster_colors[[omics]]),
    legend_gp = gpar(fill = cluster_colors[[omics]]),
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)
  )
})

# 绘制图形
setwd("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/")
pdf("intergration of 7 omics clusters.pdf", width = 16, height = 12)
draw(
  ht,
  heatmap_legend_list = legend_list,
  heatmap_legend_side = "bottom",
  padding = unit(c(0.5, 1, 0.5, 1), "cm")    # 调整边距：上、右、下、左

)
dev.off()


# 绘制百分比柱状图 ----------------------------------------------------------------
setwd("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/Percent Stacked Barplot")
#1.MET clusters
#数据准备
met_cluster_data <- data.frame(
  Group = rep(c("LUAD_male", "LUAD_female", "LUSC_male", "LUSC_female"), each = 4),
  MET_clusters = rep(1:4, times = 4),
  Value = c(0,133,0,5,     # LUAD_male
            160,4,2,1,     # LUAD_female
            1,14,0,171,    # LUSC_male
            11,0,46,3)     # LUSC_female
)
#计算百分比（按MET_cluster分组计算）
met_cluster_data <- met_cluster_data %>%
  group_by(MET_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(met_cluster_data, 
       aes(x = factor(MET_clusters), 
           y = Percentage, 
           fill = Group)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(
  #aes(label = sprintf("%.1f%%", Percentage)), 
            aes(
              label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), ""),
              y = ifelse(Percentage <5, Percentage + 0.3, Percentage)  # 小占比标签微调位置
            ),# 显示一位小数
             position = position_stack(vjust = 0.5),
            color = "black",
             size = 6) +
  scale_fill_manual(
    values = c(
      "LUAD_male" = "#ACCAEB",    
      "LUAD_female" = "#E1929D",  
      "LUSC_male" = "#8AB17C",    
      "LUSC_female" = "#E66D50"  
    ),
    breaks = c("LUAD_male", "LUAD_female", "LUSC_male", "LUSC_female") # 控制图例顺序
  )+
  labs(
    title = "MET Clusters Composition by TCGA_Project and Sex",
    x = "MET Cluster",
    y = "Percentage (%)",
    fill = "Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("met_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)

#2.MOVICS clusters
#数据准备
movics_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each = 5),
  MOVICS_clusters = rep(1:5, times = 2),
  Value = c(7,4,123,114,57,     # LUAD
            109,109,17,11,0     #LUSC
             )    
)
#计算百分比
movics_cluster_data <- movics_cluster_data %>%
  group_by(MOVICS_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(movics_cluster_data, 
       aes(x = factor(MOVICS_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(
    label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "MOVICS Clusters Composition by TCGA_Project",
    x = "MOVICS Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("movics_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)


#3.GEXP clusters
#数据准备
gexp_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each =4),
  GEXP_clusters = rep(1:4, times = 2),
  Value = c(132,13,158,2,     # LUAD
            16,221,12,0    #LUSC
  )    
)
#计算百分比
gexp_cluster_data <- gexp_cluster_data %>%
  group_by(GEXP_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(gexp_cluster_data, 
       aes(x = factor(GEXP_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(
    label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "GEXP Clusters Composition by TCGA_Project",
    x = "GEXP Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("gexp_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)


#4.MUTA clusters
#数据准备
mut_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each =4),
  MUT_clusters = rep(1:4, times = 2),
  Value = c(75,89,114,27,     # LUAD
            162,39,45,0    #LUSC
  )    
)
#计算百分比
mut_cluster_data <- mut_cluster_data %>%
  group_by(MUT_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(mut_cluster_data, 
       aes(x = factor(MUT_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "MUT Clusters Composition by TCGA_Project",
    x = "MUT Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("mut_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)


#5.CNV clusters
#数据准备
cnv_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each =5),
  CNV_clusters = rep(1:5, times = 2),
  Value = c(113,86,44,49,13,     # LUAD
            93,12,15,18,137    #LUSC
  )    
)
#计算百分比
cnv_cluster_data <- cnv_cluster_data %>%
  group_by(CNV_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(cnv_cluster_data, 
       aes(x = factor(CNV_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数,  # 显示一位小数
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "CNV Clusters Composition by TCGA_Project",
    x = "CNV Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("cnv_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)



#6.PROT clusters
#数据准备
prot_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each =4),
  PROT_clusters = rep(1:4, times = 2),
  Value = c(85,188,31,1,     # LUAD
            89,25,129,2    #LUSC
  )    
)
#计算百分比
prot_cluster_data <- prot_cluster_data %>%
  group_by(PROT_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(prot_cluster_data, 
       aes(x = factor(PROT_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数  
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "PROT Clusters Composition by TCGA_Project",
    x = "PROT Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("prot_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)


#7.TME clusters
#数据准备
tme_cluster_data <- data.frame(
  Project_id = rep(c("TCGA_LUAD", "TCGA_LUSC"), each =4),
  TME_clusters = rep(1:4, times = 2),
  Value = c(137,27,122,19,     # LUAD
            111,32,51,51    #LUSC
  )    
)
#计算百分比
tme_cluster_data <- tme_cluster_data %>%
  group_by(TME_clusters) %>%
  mutate(
    Total = sum(Value),
    Percentage = ifelse(Total > 0, round(Value/Total*100, 1), 0)
  ) %>%
  ungroup()
# 绘制百分比柱状图
ggplot(tme_cluster_data, 
       aes(x = factor(TME_clusters), 
           y = Percentage, 
           fill = Project_id)) +
  geom_col(width = 0.7, 
           color = "white", 
           linewidth = 0.3) +  # 白色分隔线
  geom_text(aes(label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "") ),# 显示一位小数
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 6) +
  scale_fill_manual(
    values = c(
      "TCGA_LUAD" = "#FCE9AC",    
      "TCGA_LUSC" = "#C7B0C7"
    ),
    breaks = c("TCGA_LUAD", "TCGA_LUSC") # 控制图例顺序
  )+
  labs(
    title = "TME Clusters Composition by TCGA_Project",
    x = "TME Cluster",
    y = "Percentage (%)",
    fill = "Project_id"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 优化Y轴间距
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题（如 "group"）
    legend.text = element_text(size = 10),  
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black"),  # Y轴刻度标签
    panel.background = element_rect(fill = "white"),  # 背景设为纯白
    plot.background = element_rect(fill = "white")     # 整体绘图区白背景（防止导出时留灰边）
  )
ggsave("tme_clusters_percent_barplot.png",width = 12,height = 12,dpi = 600)


# 生存分析 --------------------------------------------------------------------
# all(rownames(cmoic.nsclc) == rownames(commonsample.survival ))  
# surv.intergration <- data.frame(
#   OS.time  = as.numeric(commonsample.survival$OS.time),  # 生存时间（天）
#   OS.stat   = as.numeric(commonsample.survival$OS),       # 生存状态（0/1）
#   PFI.time = as.numeric(commonsample.survival$PFI.time),
#   PFI.stat =as.numeric(commonsample.survival$PFI),
#   cluster = cmoic.nsclc$clust.res$clust,                       # 分子亚型（CS1-CS5）
#   stage   = commonsample.clinical$ajcc_pathologic_stage.diagnoses,  # 临床分期（Stage I-IV）
#   row.names = rownames(commonsample.survival)
# )
# surv.info$OS.time_month <- surv.info$OS.time / 30.4368  # 使用365/12≈30.4368天/月的精确转换
# surv.info$PFI.time_month <- surv.info$PFI.time / 30.4368  # 使用365/12≈30.4368天/月的精确转换





setwd("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS")
save.image("intergartion and analysis of clusters.RData")