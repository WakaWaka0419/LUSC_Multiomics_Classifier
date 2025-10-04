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
