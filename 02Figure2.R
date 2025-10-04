rm(list = ls())
library(MOVICS)
library(tidyverse)
getwd()
Figure.path <- "./Output/Figure/"
OutData.path <- "./Output/Data/"
OutTable.path <- "./Output/Table/"
InData.path <- "Input/"
#Input data
mo.data <- readRDS("./Output/Data/mo.data.Rds")
cmoic <- readRDS("./Output/Data/cmoic.RDS")
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

#RunGsea
## run GSEA to identify down-regulated GO pathways using results from DESeq2
load("Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
load("Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_counts.rdata")
colnames(mrna_expr_fpkm) <- str_sub(colnames(mrna_expr_fpkm),1,16) 
colnames(mrna_expr_counts) <- str_sub(colnames(mrna_expr_counts),1,16)
runDEA(dea.method = "edger",
       expr       = mrna_expr_counts, # normalized expression data
       moic.res   = cmoic,
       prefix     = "TCGA-LUSC",
       res.path   = OutTable.path)

gsea.dn <- runGSEA(moic.res     = cmoic,
                   dea.method   = "edger",
                   prefix       = "TCGA-LUSC",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = mrna_expr_fpkm,
                   dat.path = OutTable.path,
                   res.path = OutTable.path,
                   fig.path = Figure.path,
                   dirct        = "up",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   n.path = 15,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "02 UP PATHWAY HEATMAP_using_upregulated_pathways") 
# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
gsva.res <- 
  runGSVA(moic.res      = cmoic,
          norm.expr     = mrna_expr_fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          #annCol        = annCol,
          #annColors     = annColors,
          fig.path      = Figure.path,
          fig.name      = "02 GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8,
          color         = c("#6699CC", "white"  , "#FF3C38"))
#TIME
library(GSVA)
library(ComplexHeatmap) # 用到里面的pheatmap画图然后拼图，需安装2.8以上版本的ComplexHeatmap
source("Script/pheatmap_translate.R") # 如果不想放弃较老版本的R及其对应的老ComplexHeatmap，可以只加载新版ComplexHeatmap里的这个函数，该脚本出自2.9.4版ComplexHeatmap
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(estimate)
source("Script/annTrackScale.R") # 加载函数，对数值进行标准化并对极端值做截断
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#Input
### 利用PCA计算MeTIL得分
annCol.tcga <- cmoic$clust.res
###Methy.expr 
load("./Input/TCGA_Data/output_methy/TCGA-LUSC_probe_info.rdata")
load("./Input/TCGA_Data/output_methy/TCGA-LUSC_beta_expr_and_pd.rdata")
methy_expr <- beta_expr %>%
  as.data.frame() %>%
  filter(apply(., 1, function(x) sum(is.na(x)) <= 0)) %>%
  rownames_to_column(var = "id") 
methy_expr <- aggregate(x = methy_expr,by = list(methy_expr$id), FUN = mean) 
methy_expr <- methy_expr %>%
  setNames(str_sub(names(.), 1, 16)) %>%
  column_to_rownames(var = "Group.1") %>%
  select(.,-id)
#匹配亚型
meth <- methy_expr[,rownames(annCol.tcga)] 
MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
meth.metil <- meth[MeTIL.marker,]
MeTIL <- meth[MeTIL.marker,]
MeTIL <- t(scale(t(MeTIL)))
#计算MeTIL得分
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL <- pca.MeTIL$rotation[,1]
annCol.tcga$MeTIL <- MeTIL
annCol.tcga <- annCol.tcga[,-1]

#immune
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
immune.signature <- read.table("Input/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9") 

# 免疫细胞的排序
immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# 计算immune/stromal得分
indata <- log2(mrna_expr_tpm + 1)
# 保存到文件
write.table(indata,file = "Output/Data/TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f = "Output/Data/TCGA_log2TPM_hugo.txt" , output.f="Output/Data/TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")
estimateScore("Output/Data/TCGA_log2TPM_hugo_ESTIMATE.txt","Output/Data/TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")
est.tcga <- read.table(file = "Output/Data/TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)

# 对数值进行标准化并对极端值做截断
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(mrna_expr_tpm)

tcga.immune.gsva <- gsva(as.matrix(log2(mrna_expr_tpm + 1)),
                         immune.sig.ccr,
                         method = "gsva")
colnames(annCol.tcga)[1] <- "CMOIC"

#figure
# 设置颜色
clust.col <- c("#2EC4B6", "#E71D36")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$CMOIC),]
annColors.tcga <- list() # 构建热图的图例颜色
annColors.tcga[["CMOIC"]] <- c("1" = clust.col[1],
                               "2" = clust.col[2]
)
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)

## 热图1：免疫检查点基因表达
colnames(mrna_expr_tpm) <- substr(colnames(mrna_expr_tpm),1,16) 
indata <- log2(mrna_expr_tpm[intersect(rownames(mrna_expr_tpm),imm.targets),] + 1)
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = annCol.tcga[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("CMOIC","StromalScore","ImmuneScore")],
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 0.6, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类

#pdf("CheckPoint_heatmap.pdf",width = 10,height = 10)
hm1
#dev.off()

## 热图2：肿瘤免疫微环境富集得分
colnames(tcga.immune.gsva) <- substr(colnames(tcga.immune.gsva),1,16)
hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1), # 富集得分标准化并排序
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22), # 根据不同类别的免疫细胞分割
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)

#pdf("TIMEheatmap.pdf",width = 10,height = 10)
hm2
#dev.off()

## 热图3：MeTIL得分
hm3 <- pheatmap(standarize.fun(t(annCol.tcga[,"MeTIL",drop = F]),halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(c(cyan,"black","#F12AFE"),64),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "MeTIL",
                cluster_rows = F,
                cluster_cols = F)
hm3
draw(hm1 %v% hm2 %v% hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
pdf("Output/Figure/02 TIME.pdf",width = 10,height = 10)
draw(hm1 %v% hm2 %v% hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
invisible(dev.off())


#Submap
library(pheatmap)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
mo.data <- readRDS("./Output/Data/mo.data.Rds")
#表达矩阵
dat <- log(mrna_expr_fpkm + 1)
colnames(dat) <- substr(colnames(dat),1,16)
#分组信息
ann <- read.table("Output/Data/group.txt",sep = "\t",header = T,row.names = 1)

TIDE <- dat  %>%
  dplyr::select(.,rownames(ann))
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,file.path(OutData.path,"TIDE_input.self_subtract"),sep = "\t",row.names = T,col.names = NA,quote = F)

#------------------------------------#
# 请参照文件夹中的TIDE使用教程.docx完成该部分#
#------------------------------------#

# 参照文件夹中TIDE使用教程得到输出文件TIDE_output.csv
TIDE.res <- read.csv("Input/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ann$TIDE <- TIDE.res[rownames(ann),"Responder"]
print(table(ann$TIDE,ann$clust))
# 检验免疫治疗响应性和亚型是否相关，p<0.05表示相关
print(fisher.test(table(ann$TIDE,ann$clust)))

# 自定义函数用来产生submap需要的数据格式
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("Input/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table("Input/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# 创建submap需要的数据格式 (TCGA)
load("Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
low_expr_counts <- apply(mrna_expr_fpkm, 1, function(x) sum(x < 1))
threshold <- ncol(mrna_expr_fpkm) * 0.9
tmp <- mrna_expr_fpkm[low_expr_counts < threshold, ]
colnames(tmp) <- substr(colnames(tmp),1,16)
tmp <- tmp %>%
  dplyr::select(.,rownames(ann))
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# 提出亚型的样本，顺序排列
samples.C1 <- rownames(ann[which(ann$clust == "1"),])
samples.C2 <- rownames(ann[which(ann$clust == "2"),])

sam_info <- data.frame("Clust"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# 产生输出数据的文件名
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"

# 把submap结果/130452/SubMap_SubMapResult.txt文件中的值填入相应的位置
# 输入文件中的名义p值和校正p值绘制热图
tmp <- matrix(c(0.25,0.68,0.66,0.94,0.84,0.05,0.82,0.02,
                1,1,1,1,1,0.44,1,0.19), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("C1 Pvalue","C2 Pvalue","C1 Bonferroni","C2 Bonferroni"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

A <- pheatmap(tmp, cellwidth = 30, cellheight = 30,
              cluster_rows = F,cluster_cols = F,
              color = heatmap.YlGnPe[5:1],
              gaps_row = 2,
              annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
              annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
              filename = "heatmap_submap.pdf")


#### TIDE结果 ####
ann <- rownames_to_column(ann,var = "sample")
TIDE.res <- TIDE.res  %>% 
  rownames_to_column(var = "sample") %>%
  inner_join(ann,by="sample")
pdf("Output/Figure/02TIDE.score.pdf",width = 4,height = 6)
p <- wilcox.test(TIDE.res[which(TIDE.res$clust == "1"),"TIDE"],TIDE.res[which(TIDE.res$clust == "2"),"TIDE"])$p.value
# 设置颜色
jco <- c("#2EC4B6", "#E71D36")
TIDE.res$clust <- as.factor(TIDE.res$clust)
ggplot(data = TIDE.res,aes(x = clust, y = TIDE, fill = clust))+
  scale_fill_manual(values = jco) + 
  #geom_violin(alpha=0.4, position = position_dodge(width = .75),
  #            size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = F, #outlier.size = -1, 
               color=c("#2EC4B6", "#E71D36"), lwd=0.8, alpha = 0.3)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("TIDE"[50])) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=4, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3)))), # 添加P值
           color="black") + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "up",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 18)
  )

dev.off()

fishres <- fisher.test(TIDE.res$Responder,TIDE.res$clust)
fishres <- round(fishres$p.value,3)
TIDE.res$value <- 1
pdf("Output/Figure/02TIDE.responder.pdf",width = 4,height = 6)
ggplot(TIDE.res)+
  geom_bar(aes(clust, value, fill = Responder),
           stat = "identity", position = "fill")+
  scale_fill_manual(values = c("#2EC4B6", "#E71D36"),name="Responder")+
  # annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  ggtitle(str_c("pvalue = ",fishres))+
  xlab("")+
  ylab("Percent")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 18))
dev.off()

#Immue_check
#inmmue_check
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
ann <- read.table("Output/Data/group.txt",sep = "\t",header = T)
mrna_expr_fpkm <- log(mrna_expr_fpkm + 1) %>%
  dplyr::filter(.,rownames(.) %in% c('CTLA4',"PDCD1")) 
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16) 
inmmue.check <- mrna_expr_fpkm %>%
  dplyr::select(ann$sample) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(ann,by = "sample") %>%
  column_to_rownames(var = "sample")
p_CTLA4 <- wilcox.test(inmmue.check[which(inmmue.check$clust == "1"),"CTLA4"],inmmue.check[which(inmmue.check$clust == "2"),"CTLA4"])$p.value
p_PD_1 <- wilcox.test(inmmue.check[which(inmmue.check$clust == "1"),"PDCD1"],inmmue.check[which(inmmue.check$clust == "2"),"PDCD1"])$p.value
#CTLA4
pdf("Output/Figure/CTLA4_Expression.pdf",width = 4,height = 6)
# 设置颜色
jco <- c("#2EC4B6", "#E71D36")
inmmue.check$clust <- as.factor(inmmue.check$clust)
ggplot(data = inmmue.check,aes(x = clust, y = CTLA4, fill = clust))+
  scale_fill_manual(values = jco) + 
  #geom_violin(alpha=0.4, position = position_dodge(width = .75),
  #            size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = F, #outlier.size = -1, 
               color=c("#2EC4B6", "#E71D36"), lwd=0.8, alpha = 0.3)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("CTLA4 Expressionlog2(FPKM+1)")) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=3, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(p_CTLA4<0.001, "< 0.001", paste0("= ",round(p,3)))), # 添加P值
           color="black") + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "up",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 18)
  )
dev.off()

#PD-1
pdf("Output/Figure/PD-1_Expression.pdf",width = 4,height = 6)
# 设置颜色
jco <- c("#2EC4B6", "#E71D36")
inmmue.check$clust <- as.factor(inmmue.check$clust)
ggplot(data = inmmue.check,aes(x = clust, y = PDCD1, fill = clust))+
  scale_fill_manual(values = jco) + 
  #geom_violin(alpha=0.4, position = position_dodge(width = .75),
  #            size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = F, #outlier.size = -1, 
               color=c("#2EC4B6", "#E71D36"), lwd=0.8, alpha = 0.3)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("PD-1 Expressionlog2(FPKM+1)")) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=3, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(p_PD_1<0.001, "< 0.001", paste0("= ",round(p,3)))), # 添加P值
           color="black") + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "up",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 18)
  )
dev.off()


##immune cross talk
library(IOBR)
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
expr <- log(mrna_expr_tpm + 1)
#C1 
cmoic <- readRDS("./Output/Data/cmoic.RDS")
exprC1 <- expr %>%
  as.data.frame() %>%
  setNames(substr(colnames(.),1,16)) %>%
  dplyr::select(rownames(cmoic$clust.res[cmoic$clust.res$clust == 1,]))
# MCPcounter
im_mcpcounterc1 <- deconvo_tme(eset = exprC1,
                             method = "mcpcounter"
)  
# CIBERSORT
im_cibersortc1 <- deconvo_tme(eset = exprC1,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
tme_combineC1 <- im_cibersortc1[,1:23] %>% 
  inner_join(im_mcpcounterc1[,c(1,11)], by="ID") 
surv <- read.table("./Input/TCGA-LUSC.survival.tsv", sep = "\t", header = T)
survC1 <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% tme_combineC1$ID) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::rename("fustat" = "OS","futime" = "OS.time") %>%
  rownames_to_column(var = "ID")
coxC1 <- merge(survC1,tme_combineC1,by = "ID") %>%
  column_to_rownames(var = "ID")


#批量多个单因素cox回归
#cox
unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],6)
  CI<-data.frame(lower=round(Gsum$conf.int[,3],2),upper=round(Gsum$conf.int[,4],2)) %>% 
    dplyr::mutate(ci=str_c(lower,"-",upper)) %>% 
    .[,"ci"]
  lowCI<-round(Gsum$conf.int[,3],2)
  highCI<-round(Gsum$conf.int[,4],2)
  coef<-round(Gsum$coefficients[,1],2)
  lowCoef<-round(log(Gsum$conf.int[,3]),2)
  highCoef<-round(log(Gsum$conf.int[,4]),2)
  unicox<-data.frame("Characteristics"=xVar,
                     "group"=rownames(Gsum$coefficients),
                     "coef"=coef,
                     "lowCoef"=lowCoef,
                     "highCoef"=highCoef,
                     "hr"=HR,
                     "lowCI"=lowCI,
                     "highCI"=highCI,
                     "CI95"=CI,
                     "pvalue"=Pvalue)
  return(unicox)
}
coxC1 <- data.table::fread("./Output/Data/COXC1IMMUE.txt")
#选择需要进行单因素分析变量名称
##准备BaSurv
BaSurv<-Surv(time = coxC1$futime,event = coxC1$fustat)
##提取变量名
varnames<-colnames(coxC1)[4:26]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(coxC1,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
poscol <- "#FB9A99" #正相关用红色连线
negcol <- "#C6DBEF" #负相关用蓝色连线

mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") #cluster的颜色，如果有更多类，就给更多的颜色
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) #用于绘制网络图

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
colnames(coxC1) <- str_replace_all(colnames(coxC1),"_"," ")
univar <- univar[,c("Characteristics","hr","lowCI","highCI","pvalue")]
colnames(univar) <- c("ID","HR","CI_up","CI_LOW","log_rank_p")
bb <- univar
bb$ID <- str_replace_all(bb$ID,"_"," ")
#用pvalue控制节点圆的大小
bb$weight <- abs(log10(bb$log_rank_p))
#用HR标圆心点的颜色
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

tme_combineC1 <- read.table("./Output/Data/tme_combineC1.txt",sep = "\t",header = T,row.names = 1,check.names = F)
write.table(tme_combineC1,file = "./Output/Data/tme_combineC1.txt",sep = "\t")
corr <- cor(tme_combineC1, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, #hclust聚为4类，根据数据的具体情况调整
         mar = c(1,0,1,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)
cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(tme_combineC1) 
head(p.corr[, 1:5])
#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)

#去掉相关性较弱的连接
corpvlue <- corpvlue[corpvlue$pvalue < 0.0001,] #只保留pvalue < 0.0001的
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)

#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)

#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)

cellcluster <- as.data.frame(t(tme_combineC1))
#cellcluster[1:5,1:5]

hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus4 <- cutree(hc, 4)) #分4类

A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
C <- as.character(rownames(as.data.frame(subset(clus4,clus4==3))))
D <- as.character(rownames(as.data.frame(subset(clus4,clus4==4))))
cls <- list(A,B,C,D)

nodes <- as.data.frame(unlist(cls))
nodes$type <- c(rep("B",9),rep("A",4),rep("C",5),rep("D",5))
names(nodes) <- c("media","type.label")

nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes
# 合并生存分析的数据和细胞分类的数据
summary(nodes$media %in% bb$ID) #检查细胞名是否一致
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge

nodes$Fraction <- abs(nodes$weight_HR)
nodes$id <- paste("S", 01:23, sep = "")
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes

#建立nodes和links的连接id，把细胞名换成ID
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,c('Macrophages M1'='S8','NK cells activated'='S14','T cells CD4 memory activated'='S17','T cells CD8'='S20','T cells follicular helper'='S21','Dendritic cells activated'='S3','Mast cells activated'='S10','Mast cells resting'='S11','Fibroblasts'='S6','Macrophages M0'='S7','T cells CD4 memory resting'='S18'))
corpvlue$to <- revalue(corpvlue$to,c('Dendritic cells resting'='S4','Macrophages M1'='S8','Monocytes'='S12','NK cells activated'='S14','B cells memory'='S1','B cells naive'='S2','NK cells resting'='S15','T cells CD4 memory activated'='S17','T cells CD4 naive'='S19','T cells CD8'='S20','T cells follicular helper'='S21','T cells gamma delta'='S22','T cells regulatory'='S23','Dendritic cells activated'='S3','Eosinophils'='S5','Mast cells activated'='S10','Mast cells resting'='S11','Neutrophils'='S13','Fibroblasts'='S6','Macrophages M0'='S7','Macrophages M2'='S9','Plasma cells'='S16','T cells CD4 memory resting'='S18'))
(links <- corpvlue)

#利用nodes和links构建网络的input文件
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2],"C"=mycol[3],"D"=mycol[4]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重
pdf("./Output/Figure/C1Immune_network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
legend("topright", #图例的位置
       c("Cell cluster-A", "Cell cluster-B", "Cell cluster-C", "Cell cluster-D"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()

##immune cross talk
library(IOBR)
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
expr <- log(mrna_expr_tpm + 1)
#C1 
cmoic <- readRDS("./Output/Data/cmoic.RDS")
exprC1 <- expr %>%
  as.data.frame() %>%
  setNames(substr(colnames(.),1,16)) %>%
  dplyr::select(rownames(cmoic$clust.res[cmoic$clust.res$clust == 1,]))
# MCPcounter
im_mcpcounterc1 <- deconvo_tme(eset = exprC1,
                               method = "mcpcounter"
)  
# CIBERSORT
im_cibersortc1 <- deconvo_tme(eset = exprC1,
                              method = "cibersort",
                              arrays = F,
                              perm = 1000
)
tme_combineC1 <- im_cibersortc1[,1:23] %>% 
  inner_join(im_mcpcounterc1[,c(1,11)], by="ID") 
surv <- read.table("./Input/TCGA-LUSC.survival.tsv", sep = "\t", header = T)
survC1 <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% tme_combineC1$ID) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::rename("fustat" = "OS","futime" = "OS.time") %>%
  rownames_to_column(var = "ID")
coxC1 <- merge(survC1,tme_combineC1,by = "ID") %>%
  column_to_rownames(var = "ID")


#批量多个单因素cox回归
#cox
unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],6)
  CI<-data.frame(lower=round(Gsum$conf.int[,3],2),upper=round(Gsum$conf.int[,4],2)) %>% 
    dplyr::mutate(ci=str_c(lower,"-",upper)) %>% 
    .[,"ci"]
  lowCI<-round(Gsum$conf.int[,3],2)
  highCI<-round(Gsum$conf.int[,4],2)
  coef<-round(Gsum$coefficients[,1],2)
  lowCoef<-round(log(Gsum$conf.int[,3]),2)
  highCoef<-round(log(Gsum$conf.int[,4]),2)
  unicox<-data.frame("Characteristics"=xVar,
                     "group"=rownames(Gsum$coefficients),
                     "coef"=coef,
                     "lowCoef"=lowCoef,
                     "highCoef"=highCoef,
                     "hr"=HR,
                     "lowCI"=lowCI,
                     "highCI"=highCI,
                     "CI95"=CI,
                     "pvalue"=Pvalue)
  return(unicox)
}
coxC1 <- data.table::fread("./Output/Data/COXC1IMMUE.txt")
#选择需要进行单因素分析变量名称
##准备BaSurv
BaSurv<-Surv(time = coxC1$futime,event = coxC1$fustat)
##提取变量名
varnames<-colnames(coxC1)[4:26]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(coxC1,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
poscol <- "#FB9A99" #正相关用红色连线
negcol <- "#C6DBEF" #负相关用蓝色连线

mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") #cluster的颜色，如果有更多类，就给更多的颜色
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) #用于绘制网络图

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
colnames(coxC1) <- str_replace_all(colnames(coxC1),"_"," ")
univar <- univar[,c("Characteristics","hr","lowCI","highCI","pvalue")]
colnames(univar) <- c("ID","HR","CI_up","CI_LOW","log_rank_p")
bb <- univar
bb$ID <- str_replace_all(bb$ID,"_"," ")
#用pvalue控制节点圆的大小
bb$weight <- abs(log10(bb$log_rank_p))
#用HR标圆心点的颜色
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

tme_combineC1 <- read.table("./Output/Data/tme_combineC1.txt",sep = "\t",header = T,row.names = 1,check.names = F)
write.table(tme_combineC1,file = "./Output/Data/tme_combineC1.txt",sep = "\t")
corr <- cor(tme_combineC1, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, #hclust聚为4类，根据数据的具体情况调整
         mar = c(1,0,1,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)
cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(tme_combineC1) 
head(p.corr[, 1:5])
#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)

#去掉相关性较弱的连接
corpvlue <- corpvlue[corpvlue$pvalue < 0.0001,] #只保留pvalue < 0.0001的
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)

#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)

#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)

cellcluster <- as.data.frame(t(tme_combineC1))
#cellcluster[1:5,1:5]

hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus4 <- cutree(hc, 4)) #分4类

A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
C <- as.character(rownames(as.data.frame(subset(clus4,clus4==3))))
D <- as.character(rownames(as.data.frame(subset(clus4,clus4==4))))
cls <- list(A,B,C,D)

nodes <- as.data.frame(unlist(cls))
nodes$type <- c(rep("B",9),rep("A",4),rep("C",5),rep("D",5))
names(nodes) <- c("media","type.label")

nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes
# 合并生存分析的数据和细胞分类的数据
summary(nodes$media %in% bb$ID) #检查细胞名是否一致
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge

nodes$Fraction <- abs(nodes$weight_HR)
nodes$id <- paste("S", 01:23, sep = "")
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes

#建立nodes和links的连接id，把细胞名换成ID
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,c('Macrophages M1'='S8','NK cells activated'='S14','T cells CD4 memory activated'='S17','T cells CD8'='S20','T cells follicular helper'='S21','Dendritic cells activated'='S3','Mast cells activated'='S10','Mast cells resting'='S11','Fibroblasts'='S6','Macrophages M0'='S7','T cells CD4 memory resting'='S18'))
corpvlue$to <- revalue(corpvlue$to,c('Dendritic cells resting'='S4','Macrophages M1'='S8','Monocytes'='S12','NK cells activated'='S14','B cells memory'='S1','B cells naive'='S2','NK cells resting'='S15','T cells CD4 memory activated'='S17','T cells CD4 naive'='S19','T cells CD8'='S20','T cells follicular helper'='S21','T cells gamma delta'='S22','T cells regulatory'='S23','Dendritic cells activated'='S3','Eosinophils'='S5','Mast cells activated'='S10','Mast cells resting'='S11','Neutrophils'='S13','Fibroblasts'='S6','Macrophages M0'='S7','Macrophages M2'='S9','Plasma cells'='S16','T cells CD4 memory resting'='S18'))
(links <- corpvlue)

#利用nodes和links构建网络的input文件
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2],"C"=mycol[3],"D"=mycol[4]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重
pdf("./Output/Figure/C1Immune_network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
legend("topright", #图例的位置
       c("Cell cluster-A", "Cell cluster-B", "Cell cluster-C", "Cell cluster-D"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()

##immune cross talk
library(IOBR)
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
expr <- log(mrna_expr_tpm + 1)
#C2 
cmoic <- readRDS("./Output/Data/cmoic.RDS")
exprC2 <- expr %>%
  as.data.frame() %>%
  setNames(substr(colnames(.),1,16)) %>%
  dplyr::select(rownames(cmoic$clust.res[cmoic$clust.res$clust == 2,]))
# MCPcounter
im_mcpcounterC2 <- deconvo_tme(eset = exprC2,
                               method = "mcpcounter"
)  
# CIBERSORT
im_cibersortC2 <- deconvo_tme(eset = exprC2,
                              method = "cibersort",
                              arrays = F,
                              perm = 1000
)
tme_combineC2 <- im_cibersortC2[,1:23] %>% 
  inner_join(im_mcpcounterC2[,c(1,11)], by="ID") 
surv <- read.table("./Input/TCGA-LUSC.survival.tsv", sep = "\t", header = T)
survC2 <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% tme_combineC2$ID) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::rename("fustat" = "OS","futime" = "OS.time") %>%
  rownames_to_column(var = "ID")
coxC2 <- merge(survC2,tme_combineC2,by = "ID") %>%
  column_to_rownames(var = "ID")


#批量多个单因素cox回归
#cox
unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],6)
  CI<-data.frame(lower=round(Gsum$conf.int[,3],2),upper=round(Gsum$conf.int[,4],2)) %>% 
    dplyr::mutate(ci=str_c(lower,"-",upper)) %>% 
    .[,"ci"]
  lowCI<-round(Gsum$conf.int[,3],2)
  highCI<-round(Gsum$conf.int[,4],2)
  coef<-round(Gsum$coefficients[,1],2)
  lowCoef<-round(log(Gsum$conf.int[,3]),2)
  highCoef<-round(log(Gsum$conf.int[,4]),2)
  unicox<-data.frame("Characteristics"=xVar,
                     "group"=rownames(Gsum$coefficients),
                     "coef"=coef,
                     "lowCoef"=lowCoef,
                     "highCoef"=highCoef,
                     "hr"=HR,
                     "lowCI"=lowCI,
                     "highCI"=highCI,
                     "CI95"=CI,
                     "pvalue"=Pvalue)
  return(unicox)
}
#write.table(coxC2,file = "./Output/Data/COXC2IMMUE.txt",sep = "\t")
coxC2 <- data.table::fread("./Output/Data/COXC2IMMUE.txt")
#选择需要进行单因素分析变量名称
##准备BaSurv
BaSurv<-Surv(time = coxC2$futime,event = coxC2$fustat)
##提取变量名
varnames<-colnames(coxC2)[4:26]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(coxC2,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
poscol <- "#FB9A99" #正相关用红色连线
negcol <- "#C6DBEF" #负相关用蓝色连线

mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") #cluster的颜色，如果有更多类，就给更多的颜色
library(reshape2)
library(corrplot)
library(plyr)

library(igraph) #用于绘制网络图

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
colnames(coxC2) <- str_replace_all(colnames(coxC2),"_"," ")
univar <- univar[,c("Characteristics","hr","lowCI","highCI","pvalue")]
colnames(univar) <- c("ID","HR","CI_up","CI_LOW","log_rank_p")
bb <- univar
bb$ID <- str_replace_all(bb$ID,"_"," ")
#用pvalue控制节点圆的大小
bb$weight <- abs(log10(bb$log_rank_p))
#用HR标圆心点的颜色
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

tme_combineC2 <- read.table("./Output/Data/tme_combineC2.txt",sep = "\t",header = T,row.names = 1,check.names = F)
#write.table(tme_combineC2,file = "./Output/Data/tme_combineC2.txt",sep = "\t")
corr <- cor(tme_combineC2, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, #hclust聚为4类，根据数据的具体情况调整
         mar = c(1,0,1,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)
cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(tme_combineC2) 
head(p.corr[, 1:5])
#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)

#去掉相关性较弱的连接
corpvlue <- corpvlue[corpvlue$pvalue < 0.0001,] #只保留pvalue < 0.0001的
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)

#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)

#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)

cellcluster <- as.data.frame(t(tme_combineC2))
#cellcluster[1:5,1:5]

hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus4 <- cutree(hc, 4)) #分4类

A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
C <- as.character(rownames(as.data.frame(subset(clus4,clus4==3))))
D <- as.character(rownames(as.data.frame(subset(clus4,clus4==4))))
cls <- list(A,B,C,D)

nodes <- as.data.frame(unlist(cls))
nodes$type <- c(rep("B",9),rep("A",4),rep("C",5),rep("D",5))
names(nodes) <- c("media","type.label")

nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes
# 合并生存分析的数据和细胞分类的数据
summary(nodes$media %in% bb$ID) #检查细胞名是否一致
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge

nodes$Fraction <- abs(nodes$weight_HR)
nodes$id <- paste("S", 01:23, sep = "")
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes

#建立nodes和links的连接id，把细胞名换成ID
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,c('Dendritic cells resting'='S4','Macrophages M1'='S8','Monocytes'='S12','NK cells activated'='S14','B cells memory'='S1','B cells naive'='S2','NK cells resting'='S15','T cells CD4 memory activated'='S17','T cells CD4 naive'='S19','T cells CD8'='S20','T cells follicular helper'='S21','T cells gamma delta'='S22','T cells regulatory'='S23','Dendritic cells activated'='S3','Eosinophils'='S5','Mast cells activated'='S10','Mast cells resting'='S11','Neutrophils'='S13','Fibroblasts'='S6','Macrophages M0'='S7','Macrophages M2'='S9','Plasma cells'='S16','T cells CD4 memory resting'='S18'))
corpvlue$to <- revalue(corpvlue$to,c('Dendritic cells resting'='S4','Macrophages M1'='S8','Monocytes'='S12','NK cells activated'='S14','B cells memory'='S1','B cells naive'='S2','NK cells resting'='S15','T cells CD4 memory activated'='S17','T cells CD4 naive'='S19','T cells CD8'='S20','T cells follicular helper'='S21','T cells gamma delta'='S22','T cells regulatory'='S23','Dendritic cells activated'='S3','Eosinophils'='S5','Mast cells activated'='S10','Mast cells resting'='S11','Neutrophils'='S13','Fibroblasts'='S6','Macrophages M0'='S7','Macrophages M2'='S9','Plasma cells'='S16','T cells CD4 memory resting'='S18'))
(links <- corpvlue)

#利用nodes和links构建网络的input文件
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2],"C"=mycol[3],"D"=mycol[4]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重
pdf("./Output/Figure/C2Immune_network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
legend("topright", #图例的位置
       c("Cell cluster-A", "Cell cluster-B", "Cell cluster-C", "Cell cluster-D"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()



#immune_landscape
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
load(file = "./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
cmoic <- readRDS("./Output/Data/cmoic.RDS")
expr <- log2(mrna_expr_tpm+0.1)
expr <- expr[apply(expr,1,sd)>0.5,]
expr <- expr %>%
  setNames(substr(colnames(.),1,16)) %>%
  dplyr::select(cmoic$clust.res$samID)
# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr,
                             method = "mcpcounter"
)
## 
## >>> Running MCP-counter

# EPIC
im_epic <- deconvo_tme(eset = expr,
                       method = "epic",
                       arrays = F
)
## 
## >>> Running EPIC

# xCell
im_xcell <- deconvo_tme(eset = expr,
                        method = "xcell",
                        arrays = F
)
## 
## >>> Running xCell
## [1] "Num. of genes: 9854"
## Estimating ssGSEA scores for 489 gene sets.
## [1] "Calculating ranks..."
## [1] "Calculating absolute values from ranks..."
# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
## 
## >>> Running CIBERSORT

# IPS
im_ips <- deconvo_tme(eset = expr,
                      method = "ips",
                      plot = F
)
## 
## >>> Running Immunophenoscore

# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr,
                            method = "quantiseq",
                            scale_mrna = T
)
## 
## Running quanTIseq deconvolution module
## Gene expression normalization and re-annotation (arrays: FALSE)
## Removing 17 noisy genes
## Removing 15 genes with high expression in tumors
## Signature genes found in data set: 134/138 (97.1%)
## Mixture deconvolution (method: lsei)
## Deconvolution sucessful!
# ESTIMATE

im_estimate <- deconvo_tme(eset = expr,
                           method = "estimate"
)
## 
## >>> Running ESTIMATE
## [1] "Merged dataset includes 9232 genes (1180 mismatched)."
## [1] "1 gene set: StromalSignature  overlap= 136"
## [1] "2 gene set: ImmuneSignature  overlap= 139"

# TIMER
im_timer <- deconvo_tme(eset = expr
                        ,method = "timer"
                        ,group_list = rep("lusc",dim(expr)[2])
)
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID") %>%
  column_to_rownames(var = "ID")
group <- read.table("./Input/group.txt",sep = "\t",header = T)
group$sample <- substr(group$sample,1,16)
#write.table(tme_combine,file = "./Output/Data/tme_cobine.txt",sep = "\t")
tme_combine <- data.table::fread("./Output/Data/tme_cobine.txt") %>%
  column_to_rownames( var = 'V1')
immMethod <- sapply(strsplit(colnames(tme_combine),"-",fixed = T),"[",2)

library(ComplexHeatmap) 
# 最新版ComplexHeatmap好像有一些bug，使用里面的pheatmap函数会报错
# 因此我们从脚本直接加载pheatmap函数
source("pheatmap_translate.R") # 位于当前文件夹，出自ComplexHeatmap_2.7.9.tar.gz

ht_opt$message = FALSE

# 创建注释
annCol <- data.frame(Type = group$clust,
                     row.names = rownames(group),
                     stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(tme_combine),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

annColors <- list("Type" = c("1" = "#2EC4B6","2" = "#E71D36"))
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
# 数据标准化
indata <- t(tme_combine)
indata <- indata[rowSums(indata) > 0,] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序

samorder <- group[order(group$clust),]
samorder <- samorder$sample
# 拆分各算法的结果
plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "MCPcounter"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "EPIC"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "xCell"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "IPS"),,drop = F]),]
plotdata6 <- plotdata[rownames(annRow[which(annRow$Methods == "quantiseq"),,drop = F]),]
plotdata7 <- plotdata[rownames(annRow[which(annRow$Methods == "estimate"),,drop = F]),]
plotdata8 <- plotdata[rownames(annRow[which(annRow$Methods == "TIMER"),,drop = F]),]
# 分别画7次热图（参数基本同pheatmap里的pheatmap）
hm1 <- pheatmap(mat = as.matrix(plotdata1[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                annotation_col = annCol[samorder,,drop = F],
                annotation_colors = annColors,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "MCPcounter") # 为子热图的图例命名

hm2 <- pheatmap(mat = as.matrix(plotdata2[1:6,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "EPIC")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "xCell")

hm4 <- pheatmap(mat = as.matrix(plotdata4[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "CIBERSORT")

hm5 <- pheatmap(mat = as.matrix(plotdata5[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "IPS")

hm6 <- pheatmap(mat = as.matrix(plotdata6[c(1:3,5:9),samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "Quantiseq")

hm7 <- pheatmap(mat = as.matrix(plotdata7[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "Estimate")

hm8 <- pheatmap(mat = as.matrix(plotdata8[,samorder]),
                border_color = NA,
                color = heatmap.BlBkRd, 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Type)[2],
                name = "TIMER")

pdf("./Output/Figure/immune heatmap by ComplexHeatmap.pdf", width = 10,height = 15) # 保存前请注意RGUI里不能有任何显示的图像，否则不会pdf打不开
draw(hm1 %v% hm2 %v% hm6 %v% hm7 %v% hm8, # 垂直连接子热图
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部
invisible(dev.off())

library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
cmoic <- readRDS("")
# phenotype data
pdata <- 
head(pdata)

# signature
(load("signature.RData"))
head(signature)

# expression set
eset <- read.csv("easy_input_expr.csv", row.names = 1)
eset[1:3, 1:3]

# phenotype data的样本需与expression set的样本一致
pdata <- pdata[pdata$ID%in%colnames(eset),]
dim(pdata)
eset <- eset[,colnames(eset)%in%pdata$ID]
dim(eset)

# 筛选，每个geneset的基因需要起码有2个以上在expression set中，否则运行后面的代码会报错
mingenecounts <- 2
print(lapply(signature,function(x) summary(x%in%rownames(eset))))
signature <- signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE)) >= mingenecounts]