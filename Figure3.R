Sys.setenv(LANG = "en")
library(ComplexHeatmap)
library(RColorBrewer)
library(maftools)
library(dplyr)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(NMF)
library(RTN)
library(snow)
library(ComplexHeatmap)
library(ClassDiscovery)
library(RColorBrewer)
library(gplots)
library(MOVICS)
library(ComplexHeatmap) # 用于绘制热图
library(circlize) # 用于热图颜色设置
library(ChAMPdata) # 用于提供甲基化注释文件
library(data.table) # 用于读取大文件
library(genefu) # 用于获取乳腺癌PAM50分型
data("pam50.robust")
data("probe.features")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
getwd()
Figure.path <- "./Output/Figure/"
OutData.path <- "./Output/Data/"
OutTable.path <- "./Output/Table/"
InData.path <- "Input/"
Mut.fig <- "./Output/Figure/mut_fig/"
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
cmoic <- readRDS("./Output/Data/cmoic.RDS")
# compare
load(file = "./Input/TCGA_Data/output_snv/TCGA-LUSC_maf_clin.rdata")
load(file = "./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
comparetmb <- snv
snv$Tumor_Sample_Barcode <- paste0(snv$Tumor_Sample_Barcode,"-01A")
tmb.LUSC  <- compTMB(moic.res     = cmoic,
                     maf          = snv,
                     rmDup        = TRUE, # remove duplicated variants per sample
                     rmFLAGS      = FALSE, # keep FLAGS mutations
                     exome.size   = 38, # estimated exome size
                     test.method  = "nonparametric", # statistical testing method
                     fig.name     = "./Output/Table/DISTRIBUTION OF TMB AND TITV")

comparefga <- read.table("./Input/tumor_seg.txt",sep = "\t",header = T)
comparefga$Sample <- substr(comparefga$Sample,1,16)
comparefga <- comparefga[,-5]
colnames(comparefga) <- c("sample","chrom","start","end","value")
fga <- compFGA(moic.res     = cmoic,
                    segment      = comparefga,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "./Output/Figure/BARPLOT OF FGA")

drug  <- compDrugsen(moic.res    = cmoic,
                         norm.expr   = mrna_expr_fpkm[,cmoic$clust.res$samID], # double guarantee sample order
                         drugs       = c("Gemcitabine","Erlotinib"), # a vector of names of drug in GDSC
                         tissueType  = "lung", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 
#Mut
# 直接读取
mut <- read.table(file = "./Input/TCGA_Data/output_snv/onco_matrix.txt",header=T,sep="\t",check.names = F) 
mut <- ifelse(mut == "",0,1) 
mut.name <- paste0(colnames(mut),"-01A")
colnames(mut) <- mut.name
#TMB
LUSC_maf <- read.maf(snv,clin_snv,isTCGA = T)
tmb <- tmb(LUSC_maf, captureSize = 38, logScale = T) %>%
  dplyr::select("Tumor_Sample_Barcode","total_perMB_log") %>%
  dplyr::rename("log10TMB" = "total_perMB_log") %>%
  tibble::column_to_rownames(var = "Tumor_Sample_Barcode")
rownames(tmb) <- paste0(rownames(tmb),"-01A")
#mut.sig
#设置颜色
maf <- snv
maf <- as.data.frame(maf)
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
rmSilence = T #是否移除沉默突变（根据实际情况还可能移除其他SNP类型）;移除则改为T
if (rmSilence) {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP" & maf$Variant_Classification != "Silent"),]) #仅考虑SNP并移除沉默突变
} else {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP"),]) #仅使用SNP (突变签名研究本身只考虑SNP)
}
snp.count <- mut.to.sigs.input(mut.ref = maf, 
                               sample.id = "Tumor_Sample_Barcode", 
                               chr = "Chromosome", 
                               pos = "Start_Position", 
                               ref = "Reference_Allele", # 参考位点碱基
                               alt = "Tumor_Seq_Allele2", # 突变位点碱基
                               bsg = BSgenome.Hsapiens.UCSC.hg38) # hg19参考基因组
write.table(snp.count,"./Output/Data/snp.count.txt",sep = "\t",row.names = T,col.names = NA)
mutsig <- read.table("./Output/Data/mutsig.weightMatrix.txt", row.names = 1, header = T)
rownames(mutsig) <- paste0(rownames(mutsig),"-01A")
cna.region <- read.table("./Input/559127/all_lesions.conf_90.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna.gene <- read.table("./Input/559127/all_thresholded.by_genes.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 设置亚型颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
blue   <- "#5bc0eb"
red    <- "#f25f5c"

# 处理突变签名数据
subt <- as.data.frame(cmoic$clust.res[,2])%>%
   dplyr::rename("clust" = "cmoic$clust.res[, 2]")
rownames(subt) <- cmoic$clust.res$samID
mutsig <- mutsig[,c("Signature.1","Signature.2","Signature.4","Signature.5","Signature.13")] # 文章中用到3种类型的signature， SBS1 (age-related), SBS2 and SBS13 (APOBEC activity-related) and SBS5 (ERCC2 mutation-related)
mutsig$APOBEC <- mutsig$Signature.2 + mutsig$Signature.13 # APOBEC相关的签名由签名2和3叠加
mutsig$CMOIC <- subt[rownames(mutsig),"clust"] # 添加亚型结果
mutsig <- mutsig[order(mutsig$CMOIC,-mutsig$APOBEC,decreasing = F),] # 确定整个热图的排序，按照亚型升序以及APOBEC降序排列
mutsig <- na.omit(mutsig)
mutsig <- dplyr::filter(mutsig,rownames(mutsig) != "TCGA-77-8146-01A")

# 挑选要展示的基因
mutgene <- c("TP53","TTN","MUC16", "RYR2", "LRP1B",
             "USH2A" ,"SYNE1" ,"FAM135B" ,"KMT2D" ,"PAPPA2" ,
             "CNTNAP5" ,"PCDH11X" ,"NEB" ,"CNTNAP2" ,"DMD" ,"ABCA13" ,"HMCN1")

# 制作oncoprint的输入数据
onco.input <- mut[mutgene,rownames(mutsig)]
onco.input[onco.input == 1] <- "Mutated" # 二值矩阵中1记为突变
onco.input[onco.input != "Mutated"] <- "" # 非“突变”给予空值
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Mutated = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A60000", col = "#A60000")) 
  }
)
col = c("Mutated" ="#A60000") # 突变颜色，注意这里只给了主图像的图例

my_ann <- subt[rownames(mutsig),,drop = F]
my_annotation = HeatmapAnnotation(df = my_ann, 
                                  col = list(clust = c("1" = clust.col[1],
                                                       "2" = clust.col[2])))

# 突变主区域的上部注释（突变负荷柱状图）
top_anno <- anno_barplot(as.numeric(tmb[rownames(mutsig),"log10TMB"]),
                         border = FALSE,
                         gp = gpar(fill = "#3379B4",border =NA,lty="blank"), 
                         height = unit(2.5, "cm"))

# 突变主区域的上部注释（突变签名柱状图）
tmp <- mutsig[,c("Signature.2","Signature.13","Signature.1","Signature.5")] # 只取和APOBEC有关的签名
tmp$Others <- 1 - rowSums(tmp) # 计算其他签名的比例
top_anno2 <- anno_barplot(as.matrix(tmp),
                          border = FALSE,
                          gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90"), 
                                    border = NA, # 无边框
                                    lty = "blank"),
                          height = unit(2, "cm")) # 高度

tmp <- as.data.frame(t(mut[mutgene,rownames(mutsig)]))
mut.order <- names(sort(colSums(tmp),decreasing = T)) # 根据突变频数高低排序展示突变的顺序
tmp$CMOIC <- subt[rownames(tmp),"clust"]
pct <- NULL # 计算各个基因突变的百分比
for (i in mut.order) {
  tmp1 <- tmp[,c(i,"CMOIC")]
  tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$CMOIC))[2,]/sum(tmp1[,1])
  pct <- rbind.data.frame(pct,tmp1)
}
rownames(pct) <- mut.order

# 添加右侧百分比堆积柱状图
right_anno <- anno_barplot(as.matrix(pct),
                           which = "row",
                           border = FALSE,
                           gp = gpar(fill = clust.col,border=NA,lty="blank"), 
                           bar_width = 0.6,
                           width = unit(1.8, "cm"),
                           height = unit(1, "cm"))

op1 <- oncoPrint(onco.input[mut.order,rownames(my_ann)], # 排序的突变矩阵
                 alter_fun = alter_fun,  # 主区域的函数，包括各单元格大小、背景颜色等等
                 col = col, # 突变颜色
                 bottom_annotation = NULL, # 无底部注释
                 top_annotation = c(HeatmapAnnotation(TMB = top_anno), # 顶部第一个注释：TMB
                                    my_annotation, # 顶部第二个注释：亚型
                                    HeatmapAnnotation(MutSig = top_anno2)), # 顶部第三个注释：突变签名
                 column_order = rownames(my_ann), # 样本的排序，根据突变签名的顺序
                 right_annotation = rowAnnotation(PCT = right_anno), # 右侧堆叠柱状图注释
                 show_pct = T, # 展示左侧的百分比
                 column_title = "", # 不显示主题
                 show_heatmap_legend = T, # 展示图例
                 column_split = my_ann$CMOIC, # 根据亚型切分热图
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))

op1

# 选择要展示的拷贝数区域
lesion.sig <- c(
                "8p11.23-Amp",
                "17q12-Amp",
                "3p13-Del",
                "9p23-Del",
                "9p21.3-Del",
                "10q23.31-Del",
                "22q13.32-Del")

cna <- cna.region[1:(nrow(cna.region)/2),c(1,8,9:(ncol(cna.region)-1))] # 选取有效列
cna <- cna[-c(19,53,76,117),]
rownames(cna) <- paste0(gsub(" ","",cna$Descriptor),"-", substr(rownames(cna),1,3)) # 重命名行以确定扩增和缺失的位点
cna.modified <- cna[1:nrow(cna),3:ncol(cna)]
colnames(cna.modified) <- substr(colnames(cna.modified),1,16)
onco.input2 <- cna.modified[lesion.sig,rownames(mutsig)] # 选取要展示的拷贝数变异
tmp1 <- onco.input2[1:2,] # 前6个为扩增
tmp1[tmp1 == 1] <- "Gain" # 数值大于0即为Gain
tmp1[tmp1 == 2] <- "Gain"
tmp1[tmp1 == 0] <- ""

tmp2 <- onco.input2[3:7,] # 后9个为缺失
tmp2[tmp2 == 1] <- "Loss"
tmp2[tmp2 == 2] <- "Loss"
tmp2[tmp2 == 0] <- ""
onco.input2 <- rbind.data.frame(tmp1,tmp2)

alter_fun2 = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = red, col = red)) 
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = blue, col = blue)) 
  }
)
col2 = c("Gain" = red,
         "Loss" = blue)
# 确定展示的顺序（看自己喜好，我这里是按照臂的顺序来的）
lesion.order <- c("3p13-Del",
                  "8p11.23-Amp",
                  "9p23-Del",
                  "9p21.3-Del",
                  "10q23.31-Del",
                  "17q12-Amp",
                  "22q13.32-Del")
tmp <- as.data.frame(t(cna.modified[lesion.order,rownames(mutsig),]))
tmp[tmp > 0] <- 1 # 所有大于1的均改为1以便计算变异频数
tmp$CMOIC <- as.character(subt[rownames(tmp),"clust"])
pct <- NULL
for (i in lesion.order) {
  tmp1 <- tmp[,c(i,"CMOIC")]
  tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$CMOIC))[2,]/sum(tmp1[,1])
  pct <- rbind.data.frame(pct,tmp1)
}
rownames(pct) <- lesion.sig
order

# 右侧堆叠百分比柱状图
right_anno2 <- anno_barplot(as.matrix(pct),
                            which = "row",
                            border = FALSE,
                            gp = gpar(fill = clust.col,
                                      border = NA,
                                      lty = "blank"), 
                            bar_width = 0.6,
                            width = unit(1.8, "cm"),
                            height = unit(1, "cm"))


# 同样的方式绘制热图
op2 <- oncoPrint(onco.input2[lesion.order,rownames(my_ann)], 
                 alter_fun = alter_fun2, 
                 col = col2, 
                 bottom_annotation = NULL, 
                 top_annotation = NULL,
                 column_order = rownames(my_ann),
                 right_annotation = rowAnnotation(PCT = right_anno2),
                 row_order = lesion.order, 
                 show_pct = T,
                 column_title = "", 
                 show_heatmap_legend = T, 
                 column_split = my_ann$CMOIC,
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))
op2


cna <- cna.gene
colnames(cna)[3:526] <- substr(colnames(cna)[3:526],1,16)
cna <- cna[c("IFNA1","MTAP","CDKN2A","CDKN2B"),rownames(mutsig)] # 文章筛选了4个基因
onco.input3 <- cna
# 由于上面的分析中缺失的部分不存在High balanced loss，所以直接让2也属于Gain而没有添加High balanced gain
# 这里小伙伴如果自己的数据满足要求，上面的拷贝数也可以分为4类
onco.input3[onco.input3 == 1] <- "Gain"
onco.input3[onco.input3 == 2] <- "High_balanced_gain"
onco.input3[onco.input3 == 0] <- ""
onco.input3[onco.input3 == -1] <- "Loss"
onco.input3[onco.input3 == -2] <- "High_balanced_loss"

alter_fun3 = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[5], col = brewer.pal(6,"Paired")[5])) 
  },
  High_balanced_gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[6], col = brewer.pal(6,"Paired")[6]))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[1], col = brewer.pal(6,"Paired")[1])) 
  },
  High_balanced_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[2], col = brewer.pal(6,"Paired")[2]))
  }
)
col3 = c("Gain" = brewer.pal(6,"Paired")[5],
         "High_balanced_gain" =brewer.pal(6,"Paired")[6],
         "Loss" = brewer.pal(6,"Paired")[1],
         "High_balanced_loss" =brewer.pal(6,"Paired")[2])

op3 <- oncoPrint(onco.input3[,rownames(my_ann)], 
                 alter_fun = alter_fun3,  
                 col = col3, 
                 bottom_annotation = NULL, 
                 top_annotation = NULL,
                 column_order = rownames(my_ann), 
                 right_annotation = NULL,
                 show_pct = T, 
                 column_title = "", 
                 show_heatmap_legend=T, 
                 column_split = my_ann$CMOIC,
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))
op3

# 构建额外图例
lgd.mutsig = Legend(labels = c("SBS2","SBS13","SBS1","SBS5","Others"), 
                    title = "MutSig", 
                    legend_gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90")))

lgd.cna.region = Legend(labels = c("Gain","Loss"), 
                        title = "CNA (arm-level)", 
                        legend_gp = gpar(fill = c(red,blue)))

lgd.cna.gene = Legend(labels = c("Gain","High_balanced_gain","Loss","High_balanced_loss"), 
                      title = "CNA (gene-level)", 
                      legend_gp = gpar(fill = brewer.pal(6,"Paired")[c(5,6,1,2)]))              

lgd_list <- list(lgd.mutsig, lgd.cna.region, lgd.cna.gene)

# 合并热图
pdf("./Output/Figure/mutational landscape in TCGA.pdf", width = 10,height = 10)
draw(op1 %v% op2 %v% op3, # 垂直叠加热图
     annotation_legend_list = lgd_list) # 添加自定义的图例
invisible(dev.off())

#GRN
data("tfsData")
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
# 加载基因表达以及样本数值信息

load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_tpm.rdata")
tcgaLUSC <- mrna_expr_tpm
pheno <- read.table("./Input/group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(pheno) <- substr(rownames(pheno),1,16)
# 加载调控子
tfs <- read.table("./Input/easyinput_regulon.txt",header = T)
colnames(tcgaLUSC) <- substr(colnames(tcgaLUSC),1,16) 
tcgaLUSC <- dplyr::select(tcgaLUSC,rownames(pheno))
# 取共有的基因名
regulatoryElements <- intersect(tfs$regulon, rownames(tcgaLUSC))

# 运行TNI构建程序
# we used the R package “RTN” to reconstruct transcriptional regulatory networks (regulons)
rtni_tcgaLUSC <- tni.constructor(expData = as.matrix(log2(tcgaLUSC + 1)), # 样图计算时候没有取对数, 
                                 regulatoryElements = regulatoryElements)

# 通过置换以及bootstrap计算reference regulatory network.
# mutual information analysis and Spearman rank-order correlation deduced the possible associations between a regulator and all potential target from the transcriptome expression profile, and permutation analysis was utilized to erase associations with an FDR > 0.00001. Bootstrapping strategy removed unstable associations through one thousand times of resampling with consensus bootstrap greater than 95%. 
# 这里量力而行设置多核，或者直接单核运算
options(cluster=snow::makeCluster(spec = 8, "SOCK")) # 打开4核并行计算（不确定是不是4核，不过我windows只用4，服务器我开12）
rtni_tcgaLUSC <- tni.permutation(rtni_tcgaLUSC, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_tcgaLUSC <- tni.bootstrap(rtni_tcgaLUSC, nBootstraps = 1000)
stopCluster(getOption("cluster")) # 关闭并行计算
# 计算DPI-filtered regulatory network
# Data processing inequality filtering eliminated the weakest associations in triangles of two regulators and common targets
rtni_tcgaLUSC <- tni.dpi.filter(rtni_tcgaLUSC, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)

# 保存TNI对象以便后续分析
save(rtni_tcgaLUSC, file="./Output/Data/rtni_tcgaLUSC.RData")

# load("rtni_tcgaLUSC.RData")
# 计算每个样本的regulon活性
# Individual regulon activity was estimated by two-sided GSEA
rtnigsea_tcgaLUSC <- tni.gsea2(rtni_tcgaLUSC, regulatoryElements = regulatoryElements)
MIBC_regact <- tni.get(rtnigsea_tcgaLUSC, what = "regulonActivity")

# 保存活性对象
save(MIBC_regact,file = "MIBC_regact.RData") 
# 设置颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
blue <- "#5bc0eb"
gold <- "#ECE700"

plotdata <- standarize.fun(t(MIBC_regact$differential),halfwidth = 1.5) # 标准化regulon的活性
annCol.tcga <- pheno[order(pheno$clust),,drop = F] # 构建样本注释信息，并对亚型进行排序
annColors.tcga <- list()
annColors.tcga[["clust"]] <- c("1" = clust.col[1],
                               "2" = clust.col[2])
hcg <- hclust(distanceMatrix(as.matrix(MIBC_regact$differential[rownames(annCol.tcga),]), "euclidean"), "ward.D")
hm <- pheatmap(plotdata[hcg$order,rownames(annCol.tcga)],
               border_color = NA, # 热图单元格无边框
               color = colorpanel(64,low=blue,mid = "black",high=gold),
               cluster_rows = F, # 行不聚类
               cluster_cols = F, # 列聚类
               show_rownames = T, # 显示行名
               show_colnames = F, # 不显示列名
               gaps_col = cumsum(table(annCol.tcga$clust))[1:2], # 亚型分割
               cellwidth = 0.8, # 固定单元格宽度
               cellheight = 10, # 固定单元格高度
               name = "LUSC Regulon", # 图例名字
               annotation_col = annCol.tcga[,"clust",drop = F], # 样本注释
               annotation_colors = annColors.tcga["clust"]) # 样本注释的对应颜色

pdf("./Output/Figure/regulon heatmap.pdf", width = 8,height = 6)
draw(hm) # 输出热图
invisible(dev.off())


#PMAPscore
# 加载数据
data("maffile", "gene_Ucox_res", "sur", "gene_symbol_Entrez", "roc_data")
load("./Input/TCGA_Data/output_snv/TCGA-LUSC_maf_clin.rdata")
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
source("./Script/utils.R") # 修正的源代码
cmoic <- readRDS("./Output/Data/cmoic.RDS")
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
mrna_expr <- mrna_expr_fpkm %>%
  dplyr::select(.,cmoic$clust.res$samID)

#maffile
snv$Tumor_Sample_Barcode <- paste0(snv$Tumor_Sample_Barcode,"-01A")
clin_snv$submitter_id <- paste0(clin_snv$submitter_id,"-01A")
snv <- dplyr::filter(snv,snv$Tumor_Sample_Barcode == colnames(mrna_expr))
maffile <- read.maf(snv,clin_snv,isTCGA = T)
# 清理数据
maf_data <- maffile@data[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")]
# maf_data <- maf_data[maf_data$Hugo_Symbol %in% intersect(rownames(gene_Ucox_res), gene_symbol_Entrez$Hugo_Symbol), ]
#gene_Ucox_res
load("./Output/Data/raw_LUSC.RData")
## Elites_surv
select.tumor <- colnames(mrna_expr)
surv <- read.table("./Input/TCGA-LUSC.survival.tsv", sep = "\t", header = T)
surv <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% select.tumor) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::rename("fustat" = "OS","futime" = "OS.time")
mrna.expr <- getElites(dat       = mrna_expr,
                       method    = "cox",
                       surv.info = surv, 
                       elite.num = 1000, # 保留MAD前100的基因
                       elite.pct = 0.1) # 此时这个参数就不起作用了
gene_Ucox_res <- mrna.expr$unicox.res%>%
  column_to_rownames(var = "gene") %>%
  dplyr::select("HR","lower","upper","pvalue") %>%
  dplyr::rename("HR.95L" = "lower","HR.95H" = "upper")

# gene_Ucox_res <- gene_Ucox_res[, "HR"]
sur <- surv %>%
  dplyr::rename("STATUS" = "fustat","TIME" = "futime")

# 下载并生成KEGG通路基因集，只需运行一次
id <- unique(keggLink("pathway", "hsa"))
id <- gsub("path:", "", id)
for (x in id){
  if (!file.exists(file.path("Inputdata", "KEGG", paste0(x, ".xml")))){
    download.file(url = file.path("https://rest.kegg.jp/get", x, "kgml"),
                  destfile = file.path("Inputdata", "KEGG", paste0(x, ".xml")), quiet = T)
  }
  Sys.sleep(0.2) # 避免下载次数过多导致服务器拒绝访问
}
## 将kgml格式的通路文件转换为Rdata
makeSPIAdata(kgml.path = "InputData/KEGG", organism = "hsa", out.path = "InputData")


# 设置所用基因集为自定义基因集
data.dir = "./Input/KEGG/InputData/"

## 整理突变信息，转换为基因×样本的矩阵
mut_status <- get_mut_status(maf_data = maf_data,  # 突变信息
                             nonsynonymous = TRUE) # 矩阵中是否只保留非同义突变
## 计算各通路的突变相关评分
pfs_score <- get_pfs_score(mut_status = mut_status,                 # 突变矩阵
                           percent = 0.03,                          # 突变频率低于此数的基因不会被纳入计算
                           gene_Ucox_res = gene_Ucox_res,           # 基因预后信息
                           gene_symbol_Entrez = gene_symbol_Entrez, # Gene SYMBOL和EntrezID的对应表格
                           verbose = F,
                           data.dir = data.dir)                     # 使用自定义基因集

## 生成预后相关的签名
final_signature <- get_final_signature(pfs_score = pfs_score,  # 各通路得分
                                       sur = sur)              # 生存信息

## 计算风险得分（multiple_score列）
km_data <- neo_get_km_data(mut_sam = mut_status,                # 突变信息
                           gene_Ucox = gene_Ucox_res,           # 基因单变量cox回归系数
                           symbol_Entrez = gene_symbol_Entrez,  # Gene Symbol和EntrezID对应表
                           sur = sur,                           # 生存信息
                           path_Ucox_mul = NULL,                # 通路多变量cox回归系数，不填，由R包根据其他输入参数求出
                           sig = final_signature,               # get_final_signature生成的预后签名
                           TRAIN = T,                           # 是否由R包根据其他输入参数求出通路多变量cox回归系数
                           data.dir = data.dir)                           

## 划分高低分组
cls <- neo_get_sam_cla(km_data = km_data,     # 风险得分
                       cut_off = -0.986)      # 高低分组阈值（此阈值为函数默认值）# 设置所用基因集为自定义基因集

## 整理突变信息，转换为基因×样本的矩阵
mut_status <- get_mut_status(maf_data = maf_data,  # 突变信息
                             nonsynonymous = TRUE) # 矩阵中是否只保留非同义突变
## 计算各通路的突变相关评分
pfs_score <- get_pfs_score(mut_status = mut_status,                 # 突变矩阵
                           percent = 0.03,                          # 突变频率低于此数的基因不会被纳入计算
                           gene_Ucox_res = gene_Ucox_res,           # 基因预后信息
                           gene_symbol_Entrez = gene_symbol_Entrez, # Gene SYMBOL和EntrezID的对应表格
                           verbose = F,
                           data.dir = data.dir)                     # 使用自定义基因集

## 生成预后相关的签名
final_signature <- get_final_signature(pfs_score = pfs_score,  # 各通路得分
                                       sur = sur)              # 生存信息

## 计算风险得分（multiple_score列）
km_data <- neo_get_km_data(mut_sam = mut_status,                # 突变信息
                           gene_Ucox = gene_Ucox_res,           # 基因单变量cox回归系数
                           symbol_Entrez = gene_symbol_Entrez,  # Gene Symbol和EntrezID对应表
                           sur = sur,                           # 生存信息
                           path_Ucox_mul = NULL,                # 通路多变量cox回归系数，不填，由R包根据其他输入参数求出
                           sig = final_signature,               # get_final_signature生成的预后签名
                           TRAIN = T,                           # 是否由R包根据其他输入参数求出通路多变量cox回归系数
                           data.dir = data.dir)                           

## 划分高低分组
cls <- neo_get_sam_cla(km_data = km_data,     # 风险得分
                       cut_off = -0.986)      # 高低分组阈值（此阈值为函数默认值）

## 绘制km曲线图
get_km_survival_curve(km_data = km_data, TRAIN = TRUE, risk.table=TRUE) # 根据中位数高低分组，绘制km曲线
dev.copy2pdf(file = "km_curve.pdf", width = 5, height = 5)

## 绘制roc曲线图
roc_data$risk_score <- km_data$multiple_score
get_roc_curve(roc_data, print.auc = TRUE, main = "Objective Response") # 绘制由风险得分判断是否响应的roc曲线
dev.copy2pdf(file = "roc_curve.pdf", width = 5, height = 5)

## 绘制OncoPlot
load(file = file.path(system.file("extdata", package = "PMAPscore"), "hsaSPIA.RData"))   # PMAP包提供的基因集
load(file = file.path(data.dir, "hsaSPIA.RData"))                                        # 自定义的基因集
path_gene <- getPathGene(path.info = path.info, gene_symbol_Entrez = gene_symbol_Entrez) # 获取所有通路涉及的基因
risk_score <- km_data$multiple_score
names(risk_score) <- rownames(km_data)
cut_off <- median(risk_score)
maffile@clinical.data <- maffile@clinical.data[,-1]
get_Oncoplots(maffile = maffile,                  # 队列maf文件
              path_gene = path_gene,              # 各通路相关基因
              mut_status = mut_status,            # 突变状态矩阵
              risk_score = risk_score,            # 风险得分
              cut_off = cut_off,                  # 风险得分高低分组阈值
              final_signature = final_signature,  # 最终所使用的签名
              pathway_name = "Gap junction",
              isTCGA = T)      # 需要展示的通路名称
dev.copy2pdf(file = "Oncoplots.pdf", width = 8, height = 6)

## 绘制风险得分热图
plot.data <- km_data
plot.data <- dplyr::arrange(plot.data, plot.data$multiple_score)                  # 按照风险评分升序排列
annCol <- plot.data[, c("event", "survival", "multiple_score")]                   # 制作列注释信息
plot.data <- plot.data[, setdiff(colnames(plot.data), colnames(annCol))]          # 从绘图数据剔除列注释信息
annCol$survival = NULL
annColor <- list()                                                                # 指定列注释颜色
annColor[["event"]] = setNames(object = c("white", "black"), nm = c("0", "1"))
annColor[["multiple_score"]] = colorRampPalette(c("#F7FBFF", "#2171B5"))(100)
annColor[["clust"]] = setNames(object = c("#DD492E","#40548A"), nm = c("1", "2"))

plot.data <- t((plot.data))
range = min(apply(plot.data, 1, max))   # 查看各列最大值，选取各通路最大值的最小值，
plot.data[plot.data>range] = range      # 作为热图主矩阵的上限
group <- read.table("./Input/group.txt", sep = "\t",header = T)
group$sample <- substr(group$sample,1,16)
plot.data <- plot.data %>%
  as.data.frame(.) %>%
  dplyr::select(.,colnames(mrna_expr)) %>%
  as.matrix(.) 
annCol <- annCol %>%  
  dplyr::filter(.,rownames(.) %in% colnames(mrna_expr)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(group,by = "sample") %>%
  column_to_rownames(var = "sample")
annCol <-arrange(annCol,desc(clust))
new_order <- match(rownames(annCol), colnames(plot.data))
plot.data <- plot.data[, new_order]
pheatmap(mat = plot.data, 
         color = colorRampPalette(c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F"))(100), # 热图主矩阵的颜色
         show_rownames = T, show_colnames = F,                                                # 显示通路名，不显示行名
         cluster_cols = F, cluster_rows = T,                                                  # 对通路进行距离，不对样本进行聚类
         annotation_col = annCol, 
         annotation_colors = annColor,
         name = "PFS score")
dev.copy2pdf(file = "PFSscore.Heatmap.pdf", width = 10, height = 6)


#Immunomodulator
# 设置热图颜色
heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")

# 设置感兴趣基因集
immunomodulator <- read.table("./Input/immunomodulator.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 数据处理 #
# 该部分是为产生最终用于绘图的文件
## 表达谱
cmoic <- readRDS("./Output/Data/cmoic.RDS")
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
mrna_expr_fpkm <- log(mrna_expr_fpkm + 1)
expr <- mrna_expr_fpkm %>%
  dplyr::select(.,cmoic$clust.res$samID)
is.element(rownames(immunomodulator),rownames(expr)) # 所有基因都在表达谱内
cmoic$clust.res$clust <- ifelse(cmoic$clust.res$clust == 1,"C1","C2")
subtype <- cmoic$clust.res$clust
sinfo <- data.frame(row.names = colnames(expr),
                    subtype = subtype)
expr <- expr[rownames(immunomodulator),]

## 甲基化谱
load("./Input/TCGA_Data/output_methy/TCGA-LUSC_beta_expr_and_pd.rdata")
meth <- beta_expr
meth <- as.data.frame(na.omit(meth))
colnames(meth) <- substr(colnames(meth),1,16)
meth <- dplyr::select(meth,cmoic$clust.res$samID)
probeOfInterest <- probe.features[which(probe.features$gene %in% rownames(immunomodulator)),]
probeOfInterest <- probeOfInterest[intersect(rownames(probeOfInterest), rownames(meth)),]
is.element(rownames(immunomodulator), probeOfInterest$gene) # 有些基因没有对应的甲基化探针
meth <- meth[rownames(probeOfInterest),]
meth$gene <- probeOfInterest$gene
meth <- as.data.frame(apply(meth[,setdiff(colnames(meth), "gene")], 2, function(x) tapply(x, INDEX=factor(meth$gene), FUN=median, na.rm=TRUE)))

## 拷贝数变异 (-2,-1,0,1,2: 2 copy del, 1 copy del, no change, amplification, high-amplification)
cna <- read.table("./Input/559127/all_thresholded.by_genes.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna$gene <- sapply(strsplit(rownames(cna),"|",fixed = T),"[",1)
cna <- cna[!duplicated(cna$gene),]; cna <- cna[,setdiff(colnames(cna),"gene")]
is.element(rownames(immunomodulator),rownames(cna)) # 有些基因没有对应的拷贝数结果
cna <- cna[intersect(rownames(cna),rownames(immunomodulator)),]
cna[cna > 1] <- 1 # 统一扩增
cna[cna < -1] <- -1 # 统一缺失
colnames(cna)[2:ncol(cna)] <- substr(colnames(cna)[2:ncol(cna)],1,16)
## 提取共同样本
comsam <- intersect(colnames(expr), colnames(meth))
comsam <- intersect(comsam, colnames(cna))
sinfo <- sinfo[comsam,,drop = F]
expr <- expr[,comsam]
meth <- meth[,comsam]
cna <- cna[,comsam]

# 亚型数目（注意，亚型必须以C1，C2，C3...等命名）
(n.subt <- length(unique(sinfo$subtype))) # 获取亚型数目
subt <- unique(sinfo$subtype) # 获取亚型名

# 初始化绘图矩阵
expMat <- as.data.frame(t(expr[rownames(immunomodulator),]))
expMat$subtype <- sinfo[rownames(expMat), "subtype"]
expMat <- as.data.frame(t(apply(expMat[,setdiff(colnames(expMat), "subtype")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(expMat$subtype), 
                                         FUN = median, 
                                         na.rm = TRUE)))) # 对同一亚型内的样本取中位数
corExpMeth <- ampFreq <- delFreq <- 
  as.data.frame(matrix(NA,
                       nrow = nrow(immunomodulator),
                       ncol = n.subt, 
                       dimnames = list(rownames(immunomodulator), 
                                       unique(sinfo$subtype))))

## 表达谱与甲基化的相关性
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(expr)) | !is.element(i, rownames(meth))) { # 如果存在任意一方有缺失的基因
    corExpMeth[i,] <- NA # 则保持矩阵为NA
  } else { # 否则取出亚型样本，做表达和甲基化的相关性
    for (j in subt) {
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      expr.subset <- as.numeric(expr[i, sam])
      meth.subset <- as.numeric(meth[i, sam])
      ct <- cor.test(expr.subset, meth.subset, method = "spearman") # 这里采用speaman相关性
      corExpMeth[i, j] <- ct$estimate
    }
  }
}

## 扩增/缺失频率
for (i in rownames(immunomodulator)) {
  if(!is.element(i, rownames(cna))) { # 同理，如果存在拷贝数中缺失某基因，则保持NA
    ampFreq[i,] <- NA 
    delFreq[i,] <- NA
  } else { # 否则
    # 计算i在总样本中的频率
    ampFreqInAll <- sum(as.numeric(cna[i,]) == 1)/ncol(cna) # 总样本中扩增的数目除以总样本数
    delFreqInAll <- sum(as.numeric(cna[i,]) == -1)/ncol(cna) # 总样本中缺失的数目除以总样本数
    for (j in subt) {
      # 计算i在亚型j中的频率
      sam <- rownames(sinfo[which(sinfo$subtype == j),,drop = F])
      cna.subset <- cna[, sam]
      ampFreqInSubt <- sum(as.numeric(cna.subset[i,]) == 1)/length(sam) # 该亚型中扩增的数目除以该亚型样本数
      delFreqInSubt <- sum(as.numeric(cna.subset[i,]) == -1)/length(sam) # 该亚型中缺失的数目除以该亚型样本数
      
      ampFreqInDiff <- ampFreqInSubt - ampFreqInAll # 根据原本，用亚型特异性扩增比例减去总扩增比例
      delFreqInDiff <- delFreqInSubt - delFreqInAll # 同理
      
      ampFreq[i, j] <- ampFreqInDiff
      delFreq[i, j] <- delFreqInDiff
    }
  }
}


# 创建列注释
annCol <- data.frame(subtype = subt,
                     row.names = subt)
annCol <- annCol[order(annCol$subtype),,drop = F] # 按照亚型排序
annColors <- list()
annColors[["subtype"]] <- c("C1" = "#DD492E",
                            "C2" = "#40548A")
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = annColors,
                              gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                              simple_anno_size     = unit(3.5, "mm"), # 注释高3.5毫米
                              show_legend          = F, # 不显示亚型的图例，因为一目了然
                              show_annotation_name = F, # 不显示该注释的名称
                              border               = FALSE) # 不显示注释的外边框

# 创建行注释
annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" # 这里字符少一些，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" # 这里换行，不会挤在一起，可以后期AI
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "black","N/A" = "#888888","Stimulatory" = "#E59E02"))
left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", # 这里是行注释（默认为列）
                               gp                   = gpar(col = "grey80"), # 每个单元格边框为灰色
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), # 注释宽3.5毫米
                               show_annotation_name = F,
                               border               = F)

## 绘制表达谱热图（参数下同）
col_expr <- colorRamp2(seq(min(na.omit(expMat)), max(na.omit(expMat)), length = 5), heatmap.BlWtRd) # 创建热图颜色（将热图输入矩阵的最大最小值取5个点，分配颜色红蓝色板；注意矩阵中可能存在的NA值）
hm.expr <- Heatmap(matrix             = as.matrix(expMat),
                   col                = col_expr,
                   border             = NA, # 无热图外边框
                   rect_gp = gpar(col = "grey80"), # 热图单元格边框为灰色
                   cluster_rows       = F, # 行不聚类
                   cluster_columns    = F, # 列不聚类
                   show_row_names     = T, # 显示行名
                   row_names_side     = "left", # 行名显示在左侧
                   row_names_gp       = gpar(fontsize = 10), # 行名字号为10
                   show_column_names  = F, # 不显示列名（可后期在颜色内AI使得亚型一目了然）
                   column_names_side  = "top", # 列名显示在顶部
                   row_split          = annRow$Category, # 行按照Category进行分割（因子顺序）
                   top_annotation     = top_anno, # 热图顶部注释
                   left_annotation    = left_anno, # 热图左侧注释
                   name               = "mRNA\nExpression", # 热图颜色图例的名称
                   width              = ncol(expMat) * unit(4, "mm"), # 热图单元格宽度（稍大于高度，因为所有注释都放在底部，平衡图形纵横比）
                   height             = nrow(expMat) * unit(3.5, "mm")) # 热图单元格高度

col_corExprMeth <- colorRamp2(seq(min(na.omit(corExpMeth)), max(na.omit(corExpMeth)), length = 5), heatmap.BlWtRd)
hm.corExprMeth <- Heatmap(matrix             = as.matrix(corExpMeth),
                          col                = col_corExprMeth,
                          border             = NA,
                          rect_gp = gpar(col = "grey80"),
                          cluster_rows       = F,
                          cluster_columns    = F,
                          show_row_names     = F,
                          row_names_side     = "left",
                          row_names_gp       = gpar(fontsize = 10),
                          show_column_names  = F,
                          column_names_side  = "top",
                          row_split          = annRow$Category,
                          row_title          = NULL,
                          top_annotation     = top_anno,
                          name               = "Expression\nvs. Methylation",
                          width              = ncol(expMat) * unit(4, "mm"),
                          height             = nrow(expMat) * unit(3.5, "mm"))

col_ampFreq <- colorRamp2(seq(min(na.omit(ampFreq)), max(na.omit(ampFreq)), length = 5), heatmap.BlWtRd)
hm.ampFreq <- Heatmap(matrix             = as.matrix(ampFreq),
                      col                = col_ampFreq,
                      border             = NA,
                      rect_gp = gpar(col = "grey80"),
                      cluster_rows       = F,
                      cluster_columns    = F,
                      show_row_names     = F,
                      row_names_side     = "left",
                      row_names_gp       = gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split          = annRow$Category,
                      row_title          = NULL,
                      top_annotation     = top_anno,
                      name               = "Amplification\nFrequency",
                      width              = ncol(expMat) * unit(4, "mm"),
                      height             = nrow(expMat) * unit(3.5, "mm"))

col_delFreq <- colorRamp2(seq(min(na.omit(delFreq)), max(na.omit(delFreq)), length = 5), heatmap.BlWtRd)
hm.delFreq <- Heatmap(matrix             = as.matrix(delFreq),
                      col                = col_delFreq,
                      border             = NA,
                      rect_gp = gpar(col = "grey70"),
                      cluster_rows       = F,
                      cluster_columns    = F,
                      show_row_names     = F,
                      row_names_side     = "left",
                      row_names_gp       = gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split          = annRow$Category,
                      row_title          = NULL,
                      top_annotation     = top_anno,
                      name               = "Deletion\nFrequency",
                      width              = ncol(expMat) * unit(4, "mm"),
                      height             = nrow(expMat) * unit(3.5, "mm"))

pdf(file = "./Output/Figure/03complexheatmap of immunomodulator.pdf", width = 8,height = 12)
draw(hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq, # 水平衔接各个子热图
     heatmap_legend_side = "bottom") # 热图颜色图例显示在下方
invisible(dev.off())



pdata <- read.table("./Input/group.txt", header = T, sep = "\t")
pdata$sample   <- substr(pdata$sample,1,16)
colnames(pdata) <- c("ID","TMEcluster")
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
eset <- mrna_expr_fpkm %>%
  setNames(substr(colnames(.),1,16)) %>%
  dplyr::select(pdata$ID)
eset <- as.matrix(eset)
eset <- log(eset + 1)
eset <- as.data.frame(eset)
# 筛选，每个geneset的基因需要起码有2个以上在expression set中，否则运行后面的代码会报错
mingenecounts <- 2
print(lapply(signature,function(x) summary(x%in%rownames(eset))))
signature <- signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE)) >= mingenecounts]

goi <- names(signature)

#相关的signature-gene-set如下：
print(goi)

# 计算gene set score
for (sig in goi) {
  pdata[, sig] <- NA
  genes <- signature[[sig]]
  genes <- genes[genes %in% rownames(eset)]
  tmp <- eset[genes, , drop=FALSE]
  
  pc <- prcomp(t(tmp),retx=TRUE)
  pdata[, sig] <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(tmp)))
}
head(pdata)

# 把其中两列的名字由缩写改为全称
# 按照这个写法，你还可以修改其他列的名字，它将最后出现在图的横坐标里
colnames(pdata) <- sub("APM","Antigen_processing_machinery",
                       colnames(pdata))
colnames(pdata) <- sub("DDR","DNA_damage_response",
                       colnames(pdata))

# 用TMEscoreA和TMEscoreB计算TMEscore
pdata$TMEscore <- pdata$TMEscoreA - pdata$TMEscoreB

# 进行数据标准化
# 共用一个坐标时可以使数据更加集中
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)

# 重排列的顺序，非必须
pdata <- pdata[,c(1,2,ncol(pdata),3:c(ncol(pdata)-1))]
head(pdata)

#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")

pdata_melt <- read.csv("very_easy_input.csv", header = T)
head(pdata_melt)
pdata_melt$TMEcluster <- as.factor(pdata_melt$TMEcluster)
# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
               # color = TMEcluster
                )) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#2EC4B6", "#E71D36")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top") 


# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

 Bggsave(p, filename = "TME-relevant-sigature-boxplot.pdf", width = 13, height = 6)
