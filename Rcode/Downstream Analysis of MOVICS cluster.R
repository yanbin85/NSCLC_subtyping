library(dplyr)
library(MOVICS)
library(ggplot2)
library(survival)
library(survminer)
library(readr)
#读取MOVICS分型结果
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/outputdata/cmoic.cluster.results.rda")

#读取临床数据并提取生存数据
load("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/outputdata/commonsamples.clinical.merge.rda")




# 示例数据:TCGA的乳腺癌数据
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))

#生存分析 -------------------------------------------------------------------
setwd("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters")
#Cox比例风险回归模型
#数据准备
surv.info <- data.frame(
  futime  
  = as.numeric(commonsamples.clinical.merge$OS.time),  # 生存时间（天）
  fustat  
  = as.numeric(commonsamples.clinical.merge$OS),       # 生存状态（0/1）
  cluster 
  = cmoic.nsclc$clust.res$clust,                       # 分子亚型（CS1-CS5）
  stage   
  = commonsamples.clinical.merge$ajcc_pathologic_stage.diagnoses,  # 临床分期（Stage I-IV）
  row.names = rownames(commonsamples.clinical.merge)
)

all(rownames(cmoic.nsclc) == rownames(surv.info ))  #检查行顺序
# 将子分期合并为主分期
# # 删除所有 "Stage " 前缀（包括空格）
surv.info$stage <- gsub("^Stage\\s*", "", surv.info$stage, ignore.case = TRUE)
surv.info$stage <- gsub("Stage I", "I", surv.info$stage)
surv.info$stage <- gsub("IA", "I", surv.info$stage)
surv.info$stage <- gsub("IB", "I", surv.info$stage)
surv.info$stage <- gsub("II", "II", surv.info$stage)
surv.info$stage <- gsub("IIA", "II", surv.info$stage)
surv.info$stage <- gsub("IIB", "II", surv.info$stage)
surv.info$stage <- gsub("III", "III", surv.info$stage)
surv.info$stage <- gsub("IIIA", "III", surv.info$stage)
surv.info$stage <- gsub("IIIB", "III", surv.info$stage)
surv.info$stage <- gsub("IV", "IV", surv.info$stage)
surv.info[] <- lapply(surv.info, function(col) {
  col[trimws(col) == ""] <- NA
  return(col)
})   

#检查数据完整性
str(surv.info)
sum(is.na(surv.info))  # 处理缺失值（如删除或填补）
#[1] 32
surv.info <- na.omit(surv.info) #去除了没有生存数据和分期数据的样本，剩余533个

surv.info$cluster <- gsub("1", "CS1", surv.info$cluster)
surv.info$cluster <- gsub("2", "CS2", surv.info$cluster)
surv.info$cluster <- gsub("3",  "CS3", surv.info$cluster)
surv.info$cluster <- gsub("4", "CS4", surv.info$cluster)
surv.info$cluster <- gsub("5", "CS5", surv.info$cluster)

#转换为因子
surv.info$cluster <- factor(surv.info$cluster, levels = c("CS1", "CS2", "CS3", "CS4","CS5"))
surv.info$stage <- factor(surv.info$stage, levels = c("I", "II", "III", "IV")) 
table(surv.info$stage, useNA = "always")
# I   II  III   IV <NA> 
#   269  159   90   15    0 
table(surv.info$cluster, useNA = "always")
# CS1  CS2  CS3  CS4  CS5 <NA> 
#   112  112  134  119   56    0 

# 拟合多变量Cox模型
cox.ph <- coxph(
  Surv(time = futime, event = fustat) ~ cluster + stage,
  data = surv.info
)

# 查看模型摘要
# summary(cox.ph)
# Call:
#   coxph(formula = Surv(time = futime, event = fustat) ~ cluster + 
#           stage, data = surv.info)
# 
# n= 533, number of events= 213 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# clusterCS2 -0.05449   0.94697  0.21439 -0.254  0.79936    
# clusterCS3  0.07397   1.07677  0.21363  0.346  0.72916    
# clusterCS4  0.28109   1.32457  0.20605  1.364  0.17251    
# clusterCS5  0.10927   1.11546  0.24592  0.444  0.65681    
# stageII     0.50998   1.66526  0.16805  3.035  0.00241 ** 
#   stageIII    0.84905   2.33742  0.17927  4.736 2.18e-06 ***
#   stageIV     1.36919   3.93216  0.31540  4.341 1.42e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# clusterCS2     0.947     1.0560    0.6221     1.442
# clusterCS3     1.077     0.9287    0.7084     1.637
# clusterCS4     1.325     0.7550    0.8845     1.984
# clusterCS5     1.115     0.8965    0.6889     1.806
# stageII        1.665     0.6005    1.1979     2.315
# stageIII       2.337     0.4278    1.6449     3.321
# stageIV        3.932     0.2543    2.1191     7.296
# 
# Concordance= 0.627  (se = 0.023 )
# Likelihood ratio test= 35.35  on 7 df,   p=1e-05
# Wald test            = 38.19  on 7 df,   p=3e-06
# Score (logrank) test = 41.03  on 7 df,   p=8e-07



library(survminer)
# 基础森林图
pdf("forest plot.pdf", width = 12, height = 8)
ggforest(
  model = cox.ph,           # Cox模型对象
  data = surv.info,         # 使用的数据
  main = "Hazard Ratio (95% CI) for Cluster and Stage",
  fontsize = 1.0,          # 文字尺寸
  noDigits = 3             # 小数位数
)
dev.off()

#K-M生存曲线
surv.data <-  data.frame(
  futime  
  = as.numeric(commonsamples.clinical.merge$OS.time),  # 生存时间（天）
  fustat  
  = as.numeric(commonsamples.clinical.merge$OS),       # 生存状态（0/1）
  row.names = rownames(commonsamples.clinical.merge)
)

surv.nsclc <- compSurv(moic.res         = cmoic.nsclc,
                      surv.info        = surv.info,
                      convt.time       = "m", # 把天变成月
                      surv.median.line = "h", 
                      xyrs.est         = c(5,10), # 计算5年和10年生存率
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")


# 比较不同亚型间突变频率 -------------------------------------------------------------
#提取临床信息作为annotation
annCol    <- data.frame(project_id=commonsamples.clinical.merge$project_id.project,
                        gender=commonsamples.clinical.merge$gender.demographic,
                        stage=commonsamples.clinical.merge$ajcc_pathologic_stage.diagnoses,
                        row.names = rownames(commonsamples.clinical.merge)
)
#stage数据统一格式
annCol$stage<- gsub("Stage IA", "Stage I", annCol$stage)
annCol$stage <- gsub("Stage IB", "Stage I", annCol$stage)
annCol$stage <- gsub("Stage IIA", "Stage II", annCol$stage)
annCol$stage <- gsub("Stage IIB", "Stage II", annCol$stage)
annCol$stage <- gsub("Stage IIIA", "Stage III", annCol$stage)
annCol$stage <- gsub("Stage IIIB", "Stage III", annCol$stage)
annCol[] <- lapply(annCol, function(col) {
  col[trimws(col) == ""] <- "unknown"
  return(col)
})

#为annotation设置颜色
annColors <- list( gender =c("male" = "#2CA02C", "female" = "#D62728"),
                   project_id=c("TCGA-LUAD" =  "#D4AF37", "TCGA-LUSC" = "#6C5B7B"),
                   pathologic_stage = c("Stage I"   = "#8DA0CB", 
                                        "Stage II"  = "#E78AC3",
                                        "Stage III" = "#FC8D62",  #
                                        "Stage IV"  = "#66C2A5" ,
                                        "unknown" = "white"))

mut.brca <- compMut(moic.res     = cmoic.nsclc,
                    mut.matrix   = multiomics.data$MUTA, # 0/1矩阵
                    doWord       = TRUE, # 生成Word文档
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # 保留在至少5%的样本中突变的基因
                    p.adj.cutoff = 0.05, # 保留padj<0.05的基因
                    innerclust   = TRUE, # 在每个亚型中进行聚类
                    annCol       = annCol, 
                    annColors    = annColors, 
                    clust.col = c("#2CA02C",  "#FF9F1C", "#FFA5AB",  "#9D4EDD","#118AB2"),
                    width        = 8, 
                    height       = 6,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")


# 亚型间差异分析 -----------------------------------------------------------------
#读取tpm数据
LUAD_tpm<- read.table("/home/data/yanbin/NSCLC_Dataset/LUAD/GEXP/tpm/TCGA-LUAD.star_tpm.tsv", header = TRUE, sep = "\t",row.names=1)
LUSC_tpm <- read.table("/home/data/yanbin/NSCLC_Dataset/LUSC/GEXP/tpm/TCGA-LUSC.star_tpm.tsv", header = TRUE, sep = "\t",row.names=1)
#处理格式
all(rownames(LUAD_tpm)==rownames(LUSC_tpm))#检查行名是否一致
nsclc_tpm <- cbind(LUAD_tpm,LUSC_tpm)
colnames(nsclc_tpm) <- substr(colnames(nsclc_tpm), start = 1, stop = 15)
colnames(nsclc_tpm) <- gsub("\\.", "-", colnames(nsclc_tpm))
nsclc_tpm <- nsclc_tpm[,rownames(commonsamples.clinical.merge)]
nsclc_tpm$Ensembl_ID <- rownames(nsclc_tpm)
gene_id <- read.table("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/DESeq2 inputdata/gencode.v36.annotation.gtf.gene.probemap",header = TRUE, sep = "\t")
gene_id <- gene_id[ , c(1,2)]
nsclc_tpm <- merge(gene_id, nsclc_tpm, by.y  = "Ensembl_ID", by.x = "id" )#基因名转换
nsclc_tpm <- distinct(nsclc_tpm, gene, .keep_all = T)#去重（重复的只保留第一行）

# 可以用limma包中的avereps函数，或者也可以使用aggregate函数取平均
# library(limma)
# nsclc_counts <- avereps(nsclc_counts, nsclc_counts$gene)

#把基因名转换为行名
rownames(nsclc_tpm) <- nsclc_tpm$gene
nsclc_tpm <- nsclc_tpm[ , -c(1,2)]
dim(nsclc_tpm) # [1] 59427   551
head(nsclc_tpm)[1:5, 1:5]
# #把log2值转换为原始counts值
# nsclc_counts_int <- 2^(nsclc_counts) - 1
# nsclc_counts_int <- apply(nsclc_counts, 2, as.integer)
# rownames(nsclc_counts_int) <- rownames(nsclc_counts)

#差异分析limma方法
setwd( "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output")
runDEA(dea.method = "limma",
       expr       = nsclc_tpm, 
       moic.res   = cmoic.nsclc,
       prefix     = "TCGA-NSCLC")

# 基于limma结果识别上调的100个基因
marker.up <- runMarker(moic.res      = cmoic.nsclc,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-NSCLC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = nsclc_tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

#基于limma结果识别下调的50个基因
marker.dn <- runMarker(moic.res      = cmoic.nsclc,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-NSCLC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "down", # direction of dysregulation in expression
                       n.marker      = 50, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = nsclc_tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")



save.image("Downstream Analysis of MOVICS clusters.RData")
