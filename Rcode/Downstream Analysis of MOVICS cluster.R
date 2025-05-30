library(dplyr)
library(MOVICS)
library(ggplot2)
library(survival)
library(survminer)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)    
library(DOSE)
library(corrplot)
library(scRNAtoolVis)
library(GseaVis)
#读取MOVICS分型结果
load("~/NSCLC_subtyping/Unsupervised Clustering/MULTI-OMICS/outputdata/cmoic.cluster.results.rda")

#读取临床与生存数据
load("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/outputdata/commonsample.clinical.rda")
load("~/NSCLC_subtyping/INTERGRATION AND ANALYSIS/outputdata/commonsample.survival.rda")


# 示例数据:TCGA的乳腺癌数据
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))

#生存分析 -------------------------------------------------------------------
setwd("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/survival_output")
#Cox比例风险回归模型
#数据准备
#
all(rownames(cmoic.nsclc) == rownames(commonsample.survival ))  
surv.info <- data.frame(
  OS.time  = as.numeric(commonsample.survival$OS.time),  # 生存时间（天）
  OS.stat   = as.numeric(commonsample.survival$OS),       # 生存状态（0/1）
  PFI.time = as.numeric(commonsample.survival$PFI.time),
  PFI.stat =as.numeric(commonsample.survival$PFI),
  cluster = cmoic.nsclc$clust.res$clust,                       # 分子亚型（CS1-CS5）
  stage   = commonsample.clinical$ajcc_pathologic_stage.diagnoses,  # 临床分期（Stage I-IV）
  row.names = rownames(commonsample.survival)
)
surv.info$OS.time_month <- surv.info$OS.time / 30.4368  # 使用365/12≈30.4368天/月的精确转换
surv.info$PFI.time_month <- surv.info$PFI.time / 30.4368  # 使用365/12≈30.4368天/月的精确转换
all(rownames(cmoic.nsclc) == rownames(surv.info ))  #检查行顺序
#命名类别
surv.info$cluster <- gsub("1", "CS1", surv.info$cluster)
surv.info$cluster <- gsub("2", "CS2", surv.info$cluster)
surv.info$cluster <- gsub("3",  "CS3", surv.info$cluster)
surv.info$cluster <- gsub("4", "CS4", surv.info$cluster)
surv.info$cluster <- gsub("5", "CS5", surv.info$cluster)
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


# 使用survival包中的survfit函数计算kaplan-Meier生存估计
OS.fit <- survfit(Surv(OS.time_month, OS.stat) ~ cluster, data = surv.info)
PFI.fit <- survfit(Surv(PFI.time_month, PFI.stat) ~ cluster, data = surv.info)
#查看完整的生存表格
summary(OS.fit)$table
summary(PFI.fit)$table
# records n.max n.start events     rmean se(rmean)   median  0.95LCL  0.95UCL
# cluster=1     114   114     114     43  93.32458  16.98745 62.81869 37.29039       NA
# cluster=2     112   112     112     33 122.29048  15.74697 61.60306 49.38101       NA
# cluster=3     137   137     137     51  92.06291  16.26252 45.27414 28.41954       NA
# cluster=4     121   121     121     59  73.01746  15.41266 26.57967 18.66162 41.33155
# cluster=5      56    56      56     22 107.65940  21.37707 47.63970 29.01093       NA
# 绘制生存曲线
#OS
pdf(file = "Kaplan-Meier Analysis of Overall Survival by CS Subtype.pdf", 
    width = 8,        # 宽度（英寸）
    height = 8 )       # 高度（英寸）
ggsurvplot(OS.fit, 
           surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
           pval = T,# 显示统计检验p值
           conf.int = F,                      # 显示置信区间
           risk.table = T,                    # 显示风险表格，展示每个时间点的风险数量
           risk.table.col = "strata",         # 风险表格按照性别着色
           xlab = "Time(months)",               # x轴标签
           ylab = "Overall Survival probability",     # y轴标签
           legend.title = "Subtype", # 图例标题
           legend.labs = c("CS1", "CS2","CS3","CS4","CS5"), # 图例标签
           ggtheme = theme_minimal() +  # 使用极简主题
             theme(
               # 调整图例文字
               legend.text = element_text(size = 14, color = "black"),  # 图例项（CS1-CS5）字体大小
               legend.title = element_text(size = 16, face = "bold"),   # 图例标题（Subtype）字体大小
               # 坐标轴刻度标签字体（数字部分）
               axis.text.x = element_text(size = 12, color = "black"),  # X轴刻度文字
               axis.text.y = element_text(size = 12, color = "black"),  # Y轴刻度文字
               # 坐标轴标题字体（"OS Time(months)"和"Survival Probability"）
               axis.title.x = element_text(size = 14, face = "bold", margin = margin(t=10)),  # X轴标题
               axis.title.y = element_text(size = 14, face = "bold", margin = margin(r=10)),  # Y轴标题
               axis.line = element_line(color = "black", size = 0.8),  # 设置轴线为黑色粗线
               panel.grid.major = element_blank(),  # 移除主网格线
               panel.grid.minor = element_blank()   # 移除次网格线
             ),
           break.x.by = 12,                  # x轴刻度间隔
           palette = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5BA")) # 自定义颜色
dev.off()
#PFI
pdf(file = "Kaplan-Meier Analysis of Progression Free Interval by CS Subtype.pdf", 
    width = 8,        # 宽度（英寸）
    height = 8 )       # 高度（英寸）
ggsurvplot(PFI.fit, 
           surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
           pval = T,                          # 显示统计检验p值
           conf.int = F,                      # 不显示置信区间
           risk.table = T,                    # 显示风险表格，展示每个时间点的风险数量
           risk.table.col = "strata",         # 风险表格按照性别着色
           xlab = "Time(months)",               # x轴标签
           ylab = "Progression Free Survival probability",     # y轴标签
           legend.title = "Subtype", # 图例标题
           legend.labs = c("CS1", "CS2","CS3","CS4","CS5"), # 图例标签
           ggtheme = theme_minimal() +  # 使用极简主题
             theme(
               # 调整图例文字
               legend.text = element_text(size = 14, color = "black"),  # 图例项（CS1-CS5）字体大小
               legend.title = element_text(size = 16, face = "bold"),   # 图例标题（Subtype）字体大小
               # 坐标轴刻度标签字体（数字部分）
               axis.text.x = element_text(size = 12, color = "black"),  # X轴刻度文字
               axis.text.y = element_text(size = 12, color = "black"),  # Y轴刻度文字
               # 坐标轴标题字体
               axis.title.x = element_text(size = 14, face = "bold", margin = margin(t=10)),  # X轴标题
               axis.title.y = element_text(size = 14, face = "bold", margin = margin(r=10)),  # Y轴标题
               axis.line = element_line(color = "black", size = 0.8),  # 设置轴线为黑色粗线
               panel.grid.major = element_blank(),  # 移除主网格线
               panel.grid.minor = element_blank()   # 移除次网格线
             ),
           break.x.by = 12,                  # x轴刻度间隔
           palette = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5BA")) # 自定义颜色
dev.off()


#K-M生存曲线
# surv.data <-  data.frame(
#   futime
#   = as.numeric(commonsample.survival$OS.time),  # 生存时间（天）
#   fustat
#   = as.numeric(commonsample.survival$OS),       # 生存状态（0/1）
#   row.names = rownames(commonsample.survival)
# )
# 
# surv.nsclc <- compSurv(moic.res         = cmoic.nsclc,
#                        surv.info        = surv.info,
#                        convt.time       = "m", # 把天变成月
#                        surv.median.line = "h",
#                        xyrs.est         = c(5,10), # 计算5年和10年生存率
#                        fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")




#Cox比例风险回归模型
#检查数据完整性
str(surv.info)
sum(is.na(surv.info))  # 处理缺失值（如删除或填补）
#[1] 50
surv.info.cox <- surv.info
surv.info.cox <- na.omit(surv.info.cox) #去除了没有生存数据和分期数据的样本，剩余535个


#转换为因子
surv.info.cox$cluster <- factor(surv.info.cox$cluster, levels = c("CS1", "CS2", "CS3", "CS4","CS5"))
surv.info.cox$stage <- factor(surv.info.cox$stage, levels = c("I", "II", "III", "IV")) 
table(surv.info.cox$stage, useNA = "always")
# I   II  III   IV <NA> 
#   271  159   90   15    0 
table(surv.info.cox$cluster, useNA = "always")
# CS1  CS2  CS3  CS4  CS5 <NA> 
#   112  112  136  119   56    0
# 拟合多变量Cox模型
OS.cox.ph <- coxph(
  Surv(time = OS.time, event = OS.stat) ~ cluster + stage,
  data = surv.info.cox
)
PFI.cox.ph <- coxph(
  Surv(time = PFI.time, event = PFI.stat) ~ cluster + stage,
  data = surv.info.cox
)
# 查看模型摘要
summary(OS.cox.ph)
# Call:
#   coxph(formula = Surv(time = OS.time, event = OS.stat) ~ cluster + 
#           stage, data = surv.info.cox)
# 
# n= 535, number of events= 213 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# clusterCS2 -0.04662   0.95445  0.21432 -0.218  0.82780    
# clusterCS3  0.08210   1.08556  0.21358  0.384  0.70069    
# clusterCS4  0.28931   1.33551  0.20604  1.404  0.16028    
# clusterCS5  0.11703   1.12416  0.24569  0.476  0.63383    
# stageII     0.51386   1.67173  0.16808  3.057  0.00223 ** 
#   stageIII    0.85265   2.34586  0.17925  4.757 1.97e-06 ***
#   stageIV     1.37790   3.96656  0.31538  4.369 1.25e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# clusterCS2    0.9545     1.0477    0.6271     1.453
# clusterCS3    1.0856     0.9212    0.7143     1.650
# clusterCS4    1.3355     0.7488    0.8918     2.000
# clusterCS5    1.1242     0.8896    0.6945     1.820
# stageII       1.6717     0.5982    1.2025     2.324
# stageIII      2.3459     0.4263    1.6509     3.333
# stageIV       3.9666     0.2521    2.1378     7.360
# 
# Concordance= 0.628  (se = 0.024 )
# Likelihood ratio test= 35.72  on 7 df,   p=8e-06
# Wald test            = 38.6  on 7 df,   p=2e-06
# Score (logrank) test = 41.5  on 7 df,   p=6e-07
summary(PFI.cox.ph)
# Call:
#   coxph(formula = Surv(time = PFI.time, event = PFI.stat) ~ cluster + 
#           stage, data = surv.info.cox)
# 
# n= 535, number of events= 205 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# clusterCS2 -0.29098   0.74753  0.23457 -1.240 0.214792    
# clusterCS3  0.28209   1.32590  0.21195  1.331 0.183215    
# clusterCS4  0.47317   1.60508  0.20556  2.302 0.021345 *  
#   clusterCS5  0.02334   1.02361  0.26649  0.088 0.930210    
# stageII     0.64186   1.90000  0.16515  3.886 0.000102 ***
#   stageIII    0.72350   2.06165  0.18931  3.822 0.000132 ***
#   stageIV     1.43231   4.18836  0.34077  4.203 2.63e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# clusterCS2    0.7475     1.3377    0.4720     1.184
# clusterCS3    1.3259     0.7542    0.8752     2.009
# clusterCS4    1.6051     0.6230    1.0728     2.401
# clusterCS5    1.0236     0.9769    0.6072     1.726
# stageII       1.9000     0.5263    1.3746     2.626
# stageIII      2.0616     0.4850    1.4226     2.988
# stageIV       4.1884     0.2388    2.1477     8.168
# 
# Concordance= 0.654  (se = 0.021 )
# Likelihood ratio test= 42.2  on 7 df,   p=5e-07
# Wald test            = 44.21  on 7 df,   p=2e-07
# Score (logrank) test = 47.24  on 7 df,   p=5e-08

library(survminer)
# 基础森林图
pdf("OS.forest plot.pdf", width = 10, height = 8)
ggforest(
  model = OS.cox.ph,           # Cox模型对象
  data = surv.info,         # 使用的数据
  main = "Hazard Ratio (95% CI) for Cluster and Stage",
  fontsize = 1.0,          # 文字尺寸
  noDigits = 3             # 小数位数
)
dev.off()

pdf("PFI.forest plot.pdf", width = 10, height = 8)
ggforest(
  model = PFI.cox.ph,           # Cox模型对象
  data = surv.info,         # 使用的数据
  main = "Hazard Ratio (95% CI) for Cluster and Stage",
  fontsize = 1.0,          # 文字尺寸
  noDigits = 3             # 小数位数
)
dev.off()


# 比较不同亚型间突变频率 -------------------------------------------------------------
setwd("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/Mutation_output")
#提取临床信息作为annotation
annCol    <- data.frame(project_id=commonsample.clinical$project_id.project,
                        gender=commonsample.clinical$gender.demographic,
                        stage=commonsample.clinical$ajcc_pathologic_stage.diagnoses,
                        row.names = rownames(commonsample.clinical)
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
                   stage = c("Stage I"   = "#8DA0CB", 
                                        "Stage II"  = "#E78AC3",
                                        "Stage III" = "#FC8D62",  #
                                        "Stage IV"  = "#66C2A5" ,
                                        "unknown" = "white"))

mut.brca <- compMut(moic.res     = cmoic.nsclc,
                    mut.matrix   = multiomics.data$MUTA, # 0/1矩阵
                    test.method = "fisher",
                    p.adj.method = "BH",
                    doWord       = TRUE, # 生成Word文档
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # 保留在至少5%的样本中突变的基因
                    p.adj.cutoff = 0.05, # 保留padj<0.05的基因
                    innerclust   = TRUE, # 在每个亚型中进行聚类
                    annCol       = annCol, 
                    annColors    = annColors, 
                    clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5BA"),
                    width        = 8, 
                    height       = 6,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

# 肿瘤突变负荷分析 ----------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data",lib = "~/R/x86_64-pc-linux-gnu-library/4.4")
# BiocManager::install("TCGAbiolinks",lib = "~/R/x86_64-pc-linux-gnu-library/4.4")
# install.packages("/home/yanbin/R包/TCGAbiolinks_2.34.1.tar.gz",
#                  repos = NULL,
#                  type = "source",
#                  lib = "/home/yanbin/R/x86_64-pc-linux-gnu-library/4.4")
                  
library(TCGAbiolinks)
#下载maf文件并整理
LUAD.maf <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(LUAD.maf)
GDCprepare(LUAD.maf, save = T,save.filename = "TCGA-LUAD_SNP.Rdata")


LUSC.maf <- GDCquery(
  project = "TCGA-LUSC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(LUSC.maf)
GDCprepare(LUSC.maf, save = T,save.filename = "TCGA-LUSC_SNP.Rdata")

load("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/TCGA-LUAD_SNP.Rdata")
maf.luad<- data
load("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/TCGA-LUSC_SNP.Rdata")
maf.lusc<- data
maf_combine <- rbind(maf.luad,maf.lusc)
maf <- data.frame(Tumor_Sample_Barcode=maf_combine$Tumor_Sample_Barcode,
                  Hugo_Symbol=maf_combine$Hugo_Symbol,
                  Chromosome=maf_combine$Chromosome,
                  Start_Position=maf_combine$Start_Position,
                  End_Position=maf_combine$End_Position,
                  Variant_Classification=maf_combine$Variant_Classification,
                  Variant_Type=maf_combine$Variant_Type,
                  Reference_Allele=maf_combine$Reference_Allele,
                  Tumor_Seq_Allele1=maf_combine$Tumor_Seq_Allele1,
                  Tumor_Seq_Allele2=maf_combine$Tumor_Seq_Allele2)
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, start = 1, stop = 15)

tmb.nsclc <- compTMB(moic.res     = cmoic.nsclc,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5BA"),
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV",
                    width=6,
                    height =6)






#library(maftools)


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

#差异分析limma方法
setwd( "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output")
runDEA(dea.method = "limma",
       expr       = nsclc_tpm, 
       moic.res   = cmoic.nsclc,
       prefix     = "TCGA-NSCLC")
# data('pbmc.markers')
library(corrplot)
library(scRNAtoolVis)
#读取差异分析结果
CS1_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS1_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
CS2_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS2_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
CS3_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS3_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
CS4_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS4_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
CS5_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS5_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
#数据准备
CS1_DEG <- data.frame(p_val = CS1_vs_others$pvalue,
                      avg_log2FC = CS1_vs_others$log2fc,
                      p_val_adj = CS1_vs_others$padj,
                      cluster = "CS1",
                      gene = rownames(CS1_vs_others),
                      row.names = rownames(CS1_vs_others)
)
CS1_DEG <- 
CS2_DEG <- data.frame(p_val = CS2_vs_others$pvalue,
                      avg_log2FC = CS2_vs_others$log2fc,
                      p_val_adj = CS2_vs_others$padj,
                      cluster = "CS2",
                      gene = rownames(CS2_vs_others),
                      row.names = rownames(CS2_vs_others)
)
CS3_DEG <- data.frame(p_val = CS3_vs_others$pvalue,
                      avg_log2FC = CS3_vs_others$log2fc,
                      p_val_adj = CS3_vs_others$padj,
                      cluster = "CS3",
                      gene = rownames(CS3_vs_others),
                      row.names = rownames(CS3_vs_others)
)
CS4_DEG <- data.frame(p_val = CS4_vs_others$pvalue,
                      avg_log2FC = CS4_vs_others$log2fc,
                      p_val_adj = CS4_vs_others$padj,
                      cluster = "CS4",
                      gene = rownames(CS4_vs_others),
                      row.names = rownames(CS4_vs_others)
)
CS5_DEG <- data.frame(p_val = CS5_vs_others$pvalue,
                      avg_log2FC = CS5_vs_others$log2fc,
                      p_val_adj = CS5_vs_others$padj,
                      cluster = "CS5",
                      gene = rownames(CS5_vs_others),
                      row.names = rownames(CS5_vs_others)
)
DEG_valcano <- rbind(CS5_DEG,CS4_DEG,CS3_DEG,CS2_DEG,CS1_DEG)
DEG_valcano_filtered <- subset(DEG_valcano, 
       p_val < 0.05 & abs(avg_log2FC) > 1.2)
#多组火山图
CSs_valcano <- jjVolcano(
  diffData = DEG_valcano,
  topGeneN = 10,
  log2FC.cutoff = 1.2,
  pvalue.cutoff= 0.05,
  adjustP.cutoff=0.01,
  col.type = "updown",
  aesCol = c('#0099CC','#CC3333'),
  tile.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB"),
  cluster.order = rev(unique(DEG_valcano$cluster)),
  size  = 3.5,
  fontface = 'italic'
)
CSs_valcano
ggsave("volcano.png",width=16,height = 12,dpi=600)
ggsave("volcano.pdf",width=16,height = 12,dpi=600)


# 识别亚型特定生物标志物 -------------------------------------------------------------
# 基于limma结果识别上调的100个基因
# 选择由 log2FoldChange 排序的差异表达最多的基因作为每个亚型的生物标志物（默认情况下，每个亚型 200 个生物标志物）。这些生物标志物应通过显著性阈值（例如，标称P值< 0.05 并调整P值<0.05），并且不得与为其他亚型鉴定的任何生物标志物重叠。
setwd( "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/biomaker output")
marker.up <- runMarker(moic.res      = cmoic.nsclc,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-NSCLC", # MUST be the same of argument in runDEA()
                       dat.path      = "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output", # path of DEA files
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

#基于limma结果识别下调的100个基因
marker.dn <- runMarker(moic.res      = cmoic.nsclc,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-NSCLC", # MUST be the same of argument in runDEA()
                       dat.path      = "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output", # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "down", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = nsclc_tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")


# 基因集富集分析确定亚型特异性功能通路 GSEA gene set enrichment analysis-----------------------------------------------------------------
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释库
library(enrichplot)    # 可视化
library(DOSE)          # 支持富集分析
library(msigdbr)     #分子特征数据库的基因集
setwd( "~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/GSEA output")
## 物种设置
organism = 'hsa'    #  人类'hsa'   
OrgDb = 'org.Hs.eg.db'#人类"org.Hs.eg.db"

#1、CS1亚型富集分析
#提取亚型1 vs 其他亚型的log2FC排序基因列表
CS1_DEG_trans <- CS1_vs_others[,2,drop=FALSE]
CS1_DEG_trans$SYMBOL <- rownames(CS1_DEG_trans)

##创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
#转化id  
SYMBOL_ENTERZID <- bitr(rownames(CS1_DEG_trans), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb) 
CS1_DEG_trans <- merge(CS1_DEG_trans, SYMBOL_ENTERZID , by='SYMBOL') #按照SYMBOL合并注释信息


CS1_geneList <- CS1_DEG_trans$log2fc
names(CS1_geneList) <- CS1_DEG_trans$ENTREZID
CS1_geneList <- sort(CS1_geneList, decreasing = T)  #从大到小排序

##gsea富集 
CS1_KEGG_kk_entrez <- gseKEGG(geneList     = CS1_geneList,
                          organism     = organism, #人hsa 鼠mmu
                          pvalueCutoff = 0.25)  #实际为padj阈值可调整 
CS1_KEGG_kk <- DOSE::setReadable(CS1_KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')#转化id             

CS1_GO_kk_entrez <- gseGO(geneList     = CS1_geneList,
                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = OrgDb,
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整
CS1_GO_kk <- DOSE::setReadable(CS1_GO_kk_entrez, 
                           OrgDb=OrgDb,
                           keyType='ENTREZID')#转化id 

save(CS1_KEGG_kk_entrez, CS1_GO_kk_entrez, file = "CS1_GSEA_result.RData")

#KEGG结果可视化
##选取富集结果
CS1_kk_gse <- CS1_KEGG_kk
#CS1_kk_gse_entrez <- CS1_KEGG_kk_entrez
###条件筛选
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
CS1_kk_gse_cut <- CS1_kk_gse[CS1_kk_gse$pvalue<0.05 & CS1_kk_gse$p.adjust<0.25 & abs(CS1_kk_gse$NES)>1]
# CS1_kk_gse_cut_down <- CS1_kk_gse_cut[CS1_kk_gse_cut$NES < 0,]
# CS1_kk_gse_cut_up <- CS1_kk_gse_cut[CS1_kk_gse_cut$NES > 0,]
#选择展现NES前几个通路 
# CS1_down_gsea <- CS1_kk_gse_cut_down[tail(order(CS1_kk_gse_cut_down$NES,decreasing = T),10),]
# CS1_up_gsea <- CS1_kk_gse_cut_up[head(order(CS1_kk_gse_cut_up$NES,decreasing = T),10),]
CS1_diff_gsea <- CS1_kk_gse_cut[head(order(abs(CS1_kk_gse_cut$NES),decreasing = T),10),]
#### 经典的GSEA图 
# CS1_up_gsea$Description
# i=2
# CS1_gseap1 <- gseaplot2(CS1_kk_gse,
#                         CS1_up_gsea$ID[i],#富集的ID编号
#                     title = CS1_up_gsea$Description[i],#标题
#                     color = "red", #GSEA线条颜色
#                     base_size = 20,#基础字体大小
#                     rel_heights = c(1.5, 0.5, 1),#副图的相对高度
#                     subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
#                     ES_geom = "line", #enrichment score用线还是用点"dot"
#                     pvalue_table = T) #显示pvalue等信息
# ggsave(CS1_gseap1, filename = 'CS1_GSEA_up_1.pdf', width =10, height =8)
#### 合并 GSEA通路 
CS1_gseap <- gseaplot2(CS1_kk_gse,
                        CS1_diff_gsea$ID,#富集的ID编号
                    title = "CS1 KEGG Pathway",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = F) #显示pvalue等信息
ggsave(CS1_gseap, filename = "CS1_GSEA_KEGG.png",width =16,height =16,dpi = 600 )

#GO结果可视化
# 绘制点图（按p.adjust排序，分面显示GO类别）
dotplot(CS1_GO_kk,
        x = "NES",           # X轴为NES值（反映富集方向）
        color = "p.adjust",  # 颜色表示p.adjust
        showCategory = 10,   # 每类显示Top10通路
        split = "ONTOLOGY",  # 按GO类别分面
        font.size = 8) +  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  ggtitle("CS1 GO Enrichment Analysis") +
  scale_color_gradient(low = "red", high = "blue")+
  theme(plot.title = element_text(hjust = 0.5)) # 自定义颜色
ggsave("CS1_GO_Dotplot.png", width = 8, height = 12,dpi = 600) 


#2、CS2亚型富集分析
#提取亚型2 vs 其他亚型的log2FC排序基因列表
CS2_vs_others <- read.table("/home/yanbin/NSCLC_subtyping/Downstream Analysis of MOVICS clusters/limma output/consensusMOIC_TCGA-NSCLC_limma_test_result.CS2_vs_Others.txt",header=TRUE,sep="\t",row.names = 1)  
CS2_DEG_trans <- CS2_vs_others[,2,drop=FALSE]
CS2_DEG_trans$SYMBOL <- rownames(CS2_DEG_trans)

##创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
CS2_DEG_trans <- merge(CS2_DEG_trans, SYMBOL_ENTERZID , by='SYMBOL') #按照SYMBOL合并注释信息
CS2_geneList <- CS2_DEG_trans$log2fc
names(CS2_geneList) <- CS2_DEG_trans$ENTREZID
CS2_geneList <- sort(CS2_geneList, decreasing = T)  #从大到小排序

##gsea富集 
CS2_KEGG_kk_entrez <- gseKEGG(geneList     = CS2_geneList,
                              organism     = organism, #人hsa 鼠mmu
                              pvalueCutoff = 0.25)  #实际为padj阈值可调整 
CS2_KEGG_kk <- DOSE::setReadable(CS2_KEGG_kk_entrez, 
                                 OrgDb=OrgDb,
                                 keyType='ENTREZID')#转化id             

CS2_GO_kk_entrez <- gseGO(geneList     = CS2_geneList,
                          ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                          OrgDb        = OrgDb,
                          keyType      = "ENTREZID",
                          pvalueCutoff = 0.25)   #实际为padj阈值可调整
CS2_GO_kk <- DOSE::setReadable(CS2_GO_kk_entrez, 
                               OrgDb=OrgDb,
                               keyType='ENTREZID')#转化id 

save(CS2_KEGG_kk_entrez, CS2_GO_kk_entrez, file = "CS2_GSEA_result.RData")

#KEGG结果可视化
##选取富集结果
CS2_kk_gse <- CS2_KEGG_kk
#CS2_kk_gse_entrez <- CS2_KEGG_kk_entrez
###条件筛选
CS2_kk_gse_cut <- CS2_kk_gse[CS2_kk_gse$pvalue<0.05 & CS2_kk_gse$p.adjust<0.25 & abs(CS2_kk_gse$NES)>1]


#选择展现NES前几个通路 
CS2_diff_gsea <- CS2_kk_gse_cut[head(order(abs(CS2_kk_gse_cut$NES),decreasing = T),10),]

#### 合并 GSEA通路 
CS2_gseap <- gseaplot2(CS2_kk_gse,
                       CS2_diff_gsea$ID,#富集的ID编号
                       title = "CS2 KEGG Pathway",#标题
                       color = "red",#GSEA线条颜色
                       base_size = 20,#基础字体大小
                       rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                       subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                       ES_geom = "line",#enrichment score用线还是用点"dot"
                       pvalue_table = F) #显示pvalue等信息
ggsave(CS2_gseap, filename = "CS2_GSEA_KEGG.png",width =16,height =16,dpi = 600)


#GO结果可视化
# 绘制点图（按p.adjust排序，分面显示GO类别）
dotplot(CS2_GO_kk,
        x = "NES",           # X轴为NES值（反映富集方向）
        color = "p.adjust",  # 颜色表示p.adjust
        showCategory = 10,   # 每类显示Top10通路
        split = "ONTOLOGY",  # 按GO类别分面
        font.size = 8) +  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  ggtitle("CS2 GO Enrichment Analysis") +
  scale_color_gradient(low = "red", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5))  # 自定义颜色
ggsave("CS2_GO_Dotplot.png", width = 8, height = 12,,dpi = 600) 




#3、CS3亚型富集分析
#提取亚3 vs 其他亚型的log2FC排序基因列表
CS3_DEG_trans <- CS3_vs_others[,2,drop=FALSE]
CS3_DEG_trans$SYMBOL <- rownames(CS3_DEG_trans)

##创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
CS3_DEG_trans <- merge(CS3_DEG_trans, SYMBOL_ENTERZID , by='SYMBOL') #按照SYMBOL合并注释信息
CS3_geneList <- CS3_DEG_trans$log2fc
names(CS3_geneList) <- CS3_DEG_trans$ENTREZID
CS3_geneList <- sort(CS3_geneList, decreasing = T)  #从大到小排序

##gsea富集 
CS3_KEGG_kk_entrez <- gseKEGG(geneList     = CS3_geneList,
                              organism     = organism, #人hsa 鼠mmu
                              pvalueCutoff = 0.25)  #实际为padj阈值可调整 
CS3_KEGG_kk <- DOSE::setReadable(CS3_KEGG_kk_entrez, 
                                 OrgDb=OrgDb,
                                 keyType='ENTREZID')#转化id             

CS3_GO_kk_entrez <- gseGO(geneList     = CS3_geneList,
                          ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                          OrgDb        = OrgDb,
                          keyType      = "ENTREZID",
                          pvalueCutoff = 0.25)   #实际为padj阈值可调整
CS3_GO_kk <- DOSE::setReadable(CS3_GO_kk_entrez, 
                               OrgDb=OrgDb,
                               keyType='ENTREZID')#转化id 

save(CS3_KEGG_kk_entrez, CS3_GO_kk_entrez, file = "CS3_GSEA_result.RData")

#KEGG结果可视化
##选取富集结果
CS3_kk_gse <- CS3_KEGG_kk
###条件筛选
CS3_kk_gse_cut <- CS3_kk_gse[CS3_kk_gse$pvalue<0.05 & CS3_kk_gse$p.adjust<0.25 & abs(CS3_kk_gse$NES)>1]
#选择展现NES前几个通路 
CS3_diff_gsea <- CS3_kk_gse_cut[head(order(abs(CS3_kk_gse_cut$NES),decreasing = T),10),]
#### 合并 GSEA通路 
CS3_gseap <- gseaplot2(CS3_kk_gse,
                       CS3_diff_gsea$ID,#富集的ID编号
                       title = "CS3 KEGG Pathway",#标题
                       color = "red",#GSEA线条颜色
                       base_size = 20,#基础字体大小
                       rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                       subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                       ES_geom = "line",#enrichment score用线还是用点"dot"
                       pvalue_table = F) #显示pvalue等信息
ggsave(CS3_gseap, filename = "CS3_GSEA_KEGG.png",width =16,height =16,dpi = 600)

#GO结果可视化
# 绘制点图（按p.adjust排序，分面显示GO类别）
dotplot(CS3_GO_kk,
        x = "NES",           # X轴为NES值（反映富集方向）
        color = "p.adjust",  # 颜色表示p.adjust
        showCategory = 10,   # 每类显示Top10通路
        split = "ONTOLOGY",  # 按GO类别分面
        font.size = 8) +  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  ggtitle("CS3 GO Enrichment Analysis") +
  scale_color_gradient(
    low = "red", high = "blue",
    limits = c(0, 1))+
  
  theme(plot.title = element_text(hjust = 0.5))  
ggsave("CS3_GO_Dotplot.png", width = 8, height = 12,dpi = 600) 




#4、CS4亚型富集分析
#提取亚4 vs 其他亚型的log2FC排序基因列表

CS4_DEG_trans <- CS4_vs_others[,2,drop=FALSE]
CS4_DEG_trans$SYMBOL <- rownames(CS4_DEG_trans)

##创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
CS4_DEG_trans <- merge(CS4_DEG_trans, SYMBOL_ENTERZID , by='SYMBOL') #按照SYMBOL合并注释信息
CS4_geneList <- CS4_DEG_trans$log2fc
names(CS4_geneList) <- CS4_DEG_trans$ENTREZID
CS4_geneList <- sort(CS4_geneList, decreasing = T)  #从大到小排序

##gsea富集 
CS4_KEGG_kk_entrez <- gseKEGG(geneList     = CS4_geneList,
                              organism     = organism, #人hsa 鼠mmu
                              pvalueCutoff = 0.25)  #实际为padj阈值可调整 
CS4_KEGG_kk <- DOSE::setReadable(CS4_KEGG_kk_entrez, 
                                 OrgDb=OrgDb,
                                 keyType='ENTREZID')#转化id             

CS4_GO_kk_entrez <- gseGO(geneList     = CS4_geneList,
                          ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                          OrgDb        = OrgDb,
                          keyType      = "ENTREZID",
                          pvalueCutoff = 0.25)   #实际为padj阈值可调整
CS4_GO_kk <- DOSE::setReadable(CS4_GO_kk_entrez, 
                               OrgDb=OrgDb,
                               keyType='ENTREZID')#转化id 

save(CS4_KEGG_kk_entrez, CS4_GO_kk_entrez, file = "CS4_GSEA_result.RData")

#KEGG结果可视化
##选取富集结果
CS4_kk_gse <- CS4_KEGG_kk
###条件筛选
CS4_kk_gse_cut <- CS4_kk_gse[CS4_kk_gse$pvalue<0.05 & CS4_kk_gse$p.adjust<0.25 & abs(CS4_kk_gse$NES)>1]
#选择展现NES前几个通路 
CS4_diff_gsea <- CS4_kk_gse_cut[head(order(abs(CS4_kk_gse_cut$NES),decreasing = T),10),]
#### 合并 GSEA通路 
CS4_gseap <- gseaplot2(CS4_kk_gse,
                       CS4_diff_gsea$ID,#富集的ID编号
                       title = "CS4 KEGG Pathway",#标题
                       color = "red",#GSEA线条颜色
                       base_size = 20,#基础字体大小
                       rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                       subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                       ES_geom = "line",#enrichment score用线还是用点"dot"
                       pvalue_table = F) #显示pvalue等信息
ggsave(CS4_gseap, filename = "CS4_GSEA_KEGG.png",width =16,height =16,dpi = 600)

#GO结果可视化
# 绘制点图（按p.adjust排序，分面显示GO类别）
dotplot(CS4_GO_kk,
        x = "NES",           # X轴为NES值（反映富集方向）
        color = "p.adjust",  # 颜色表示p.adjust
        showCategory = 10,   # 每类显示Top10通路
        split = "ONTOLOGY",  # 按GO类别分面
        font.size = 8) +  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  ggtitle("CS4 GO Enrichment Analysis") +
  scale_color_gradient(low = "red", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5))  # 自定义颜色
ggsave("CS4_GO_Dotplot.png", width = 8, height = 12,dpi = 600) 


#5、CS5亚型富集分析
#提取亚型5vs 其他亚型的log2FC排序基因列表
CS5_DEG_trans <- CS5_vs_others[,2,drop=FALSE]
CS5_DEG_trans$SYMBOL <- rownames(CS5_DEG_trans)

##创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
CS5_DEG_trans <- merge(CS5_DEG_trans, SYMBOL_ENTERZID , by='SYMBOL') #按照SYMBOL合并注释信息
CS5_geneList <- CS5_DEG_trans$log2fc
names(CS5_geneList) <- CS5_DEG_trans$ENTREZID
CS5_geneList <- sort(CS5_geneList, decreasing = T)  #从大到小排序

##gsea富集 
CS5_KEGG_kk_entrez <- gseKEGG(geneList     = CS5_geneList,
                              organism     = organism, #人hsa 鼠mmu
                              pvalueCutoff = 0.25)  #实际为padj阈值可调整 
CS5_KEGG_kk <- DOSE::setReadable(CS5_KEGG_kk_entrez, 
                                 OrgDb=OrgDb,
                                 keyType='ENTREZID')#转化id             

CS5_GO_kk_entrez <- gseGO(geneList     = CS5_geneList,
                          ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                          OrgDb        = OrgDb,
                          keyType      = "ENTREZID",
                          pvalueCutoff = 0.25)   #实际为padj阈值可调整
CS5_GO_kk <- DOSE::setReadable(CS5_GO_kk_entrez, 
                               OrgDb=OrgDb,
                               keyType='ENTREZID')#转化id 

save(CS5_KEGG_kk_entrez, CS5_GO_kk_entrez, file = "CS5_GSEA_result.RData")

#KEGG结果可视化
##选取富集结果
CS5_kk_gse <- CS5_KEGG_kk
###条件筛选
CS5_kk_gse_cut <- CS5_kk_gse[CS5_kk_gse$pvalue<0.05 & CS5_kk_gse$p.adjust<0.25 & abs(CS5_kk_gse$NES)>1]
#选择展现NES前几个通路 
CS5_diff_gsea <- CS5_kk_gse_cut[head(order(abs(CS5_kk_gse_cut$NES),decreasing = T),10),]
#### 合并 GSEA通路 
CS5_gseap <- gseaplot2(CS5_kk_gse,
                       CS5_diff_gsea$ID,#富集的ID编号
                       title = "CS5 KEGG Pathway",#标题
                       color = "red",#GSEA线条颜色
                       base_size = 20,#基础字体大小
                       rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                       subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                       ES_geom = "line",#enrichment score用线还是用点"dot"
                       pvalue_table = F) #显示pvalue等信息
ggsave(CS5_gseap, filename = "CS5_GSEA_KEGG.png",width =16,height =16,dpi = 600)

#GO结果可视化
# 绘制点图（按p.adjust排序，分面显示GO类别）
dotplot(CS5_GO_kk,
        x = "NES",           # X轴为NES值（反映富集方向）
        color = "p.adjust",  # 颜色表示p.adjust
        showCategory = 10,   # 每类显示Top10通路
        split = "ONTOLOGY",  # 按GO类别分面
        font.size = 8) +  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  ggtitle("CS5 GO Enrichment Analysis") +
  scale_color_gradient(low = "red", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5))  # 自定义颜色
ggsave("CS5_GO_Dotplot.png", width = 8, height = 12,dpi = 600) 





# 合并所有亚型的GO富集结果
combined_GO <- rbind(
  cbind(CS1_GO_kk@result, Subtype = "CS1"),
  cbind(CS2_GO_kk@result, Subtype = "CS2"),
  cbind(CS3_GO_kk@result, Subtype = "CS3"),
  cbind(CS4_GO_kk@result, Subtype = "CS4"),
  cbind(CS5_GO_kk@result, Subtype = "CS5")
)

# 按亚型和GO类别筛选Top10通路（按p.adjust排序）
library(dplyr)
combined_GO_top10 <- combined_GO %>%
  group_by(Subtype, ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()
# 🛠️ 关键优化步骤：添加文本换行处理
library(stringr)
combined_GO_top10 <- combined_GO_top10 %>% 
  mutate(
    Description = str_wrap(Description, 
                           width = 70,  # 每行最多40字符
                           exdent = 2)  # 第二行缩进2字符
  )
# 确保亚型顺序正确
combined_GO_top10$Subtype <- factor(combined_GO_top10$Subtype, 
                                 levels = c("CS1", "CS2", "CS3", "CS4", "CS5"))
# 创建分面点图
library(ggplot2)
GO_plot1 <- ggplot(combined_GO_top10, 
            aes(x = NES, 
                y = reorder(Description, NES),  # 按NES排序y轴
                color = p.adjust,
                size = setSize)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red",  # 颜色映射p.adjust
                       limits = c(0, 0.003),         # 统一颜色标尺
                       name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count",        # 点大小映射基因数量
                        range = c(2, 6)) +
  facet_grid(ONTOLOGY ~ Subtype,                    # 矩阵分面
             scales = "free_y",                     # 自由y轴范围
             space = "free") +                      # 自动调整间距
  labs(x = "Normalized Enrichment Score (NES)",
       y = "GO Pathway",
       title = "Multi-Subtype GO Enrichment Analysis") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )
GO_plot1
# 保存高清图（可根据需要调整尺寸）
ggsave("Combined_GO_Dotplot.pdf", GO_plot1,   width = 20, height = 15, dpi = 600)





#多组结果对比——KEGG
allsubtype_genelist <- list(CS1_geneList,CS2_geneList,CS3_geneList,CS4_geneList,CS5_geneList)

lapply(1:5, function(x){
  CS_kegg <- gseKEGG(geneList     = allsubtype_genelist[[x]],
                organism = "hsa",
                keyType = "kegg",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 1,#设置 pvalueCutoff=1,保证在每个富集结果里都有共同的通路
                eps          = 0,  #精确计算极小p值
                verbose      = FALSE)
}) -> allsubtypes_gsea_kegg_list
#查看其中一个结果
View(as.data.frame(allsubtypes_gsea_kegg_list[[1]]))
View(as.data.frame(allsubtypes_gsea_kegg_list[[4]]))
GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04110",#cell type
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))

ggsave("allsubtypes_gsea_kegg1.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04060", #Cytokine-cytokine receptor interaction
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))

ggsave("allsubtypes_gsea_kegg2.png",width = 8,height = 8,dpi=900)


GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04659", #Th17 cell differentiation
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))

ggsave("allsubtypes_gsea_kegg3.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04145", #Phagosome
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))

ggsave("allsubtypes_gsea_kegg4.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04062", #Chemokine signaling pathway
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))

ggsave("allsubtypes_gsea_kegg5.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_kegg_list,
            geneSetID = "hsa04670", #Leukocyte transendothelial migration
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_kegg6.png",width = 8,height = 8,dpi=900)


#多组结果对比——GO
lapply(1:5, function(x){
  CS_go <- gseGO(geneList     = allsubtype_genelist[[x]],
                 OrgDb = org.Hs.eg.db,
                 keyType      = "ENTREZID",
                 ont = "BP",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                pvalueCutoff = 1,#设置 pvalueCutoff=1,保证在每个富集结果里都有共同的通路
                eps          = 0,  #精确计算极小p值
                verbose      = FALSE)
}) -> allsubtypes_gsea_go_list
#查看其中一个结果
View(as.data.frame(allsubtypes_gsea_go_list[[4]]))
#选择感兴趣的GO条目绘图
GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0002443",#leukocyte mediated immunity
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go1.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0002449",#lymphocyte mediated immunity
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go2.png",width = 8,height = 8,dpi=900)

GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0002274",#myeloid leukocyte activation
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go3.png",width = 8,height = 8,dpi=900)


GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0002377",#immunoglobulin production
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go4.png",width = 8,height = 8,dpi=900)


GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0019724",#B cell mediated immunity
            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go5.png",width = 8,height = 8,dpi=900)


GSEAmultiGP(gsea_list = allsubtypes_gsea_go_list,
            geneSetID = "GO:0002460",#adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains

            exp_name = c("CS1","CS2","CS3","CS4","CS5"),
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = "right",
            curve.col = ggsci::pal_lancet()(5))
ggsave("allsubtypes_gsea_go6.png",width = 8,height = 8,dpi=900)


setwd("~/NSCLC_subtyping/Downstream Analysis of MOVICS clusters")
save.image("Downstream Analysis of MOVICS clusters.RData")
