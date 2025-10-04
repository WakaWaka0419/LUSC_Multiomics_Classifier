library(tidyverse)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(timeROC)
library(dplyr)
library(tibble)
library(stringr)
library(survival)
library(survivalsvm)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#第一列为样本名，第二三列为生存时间和生存状态，后面的列均为基因。
cmoic <- readRDS("./Output/Data/cmoic.RDS")
load("./Output/Data/gse37745.Rdata")
load("./Output/Data/gse157010.Rdata")
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
surv <- read.table("Input/TCGA-LUSC.survival.tsv", sep = "\t", header = T)
select.tumor <- cmoic$clust.res$samID
tcga.surv <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% select.tumor) %>%
  column_to_rownames(var = "sample") 
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
TCGA.expr <- dplyr::select(mrna_expr_fpkm,select.tumor)
common.gene <- intersect(rownames(TCGA.expr),rownames(gse37745.expr))
common.gene <- intersect(common.gene,rownames(gse157010.expr))
TCGA.expr <- dplyr::filter(TCGA.expr,rownames(TCGA.expr) %in% common.gene)
gse37745.expr <- dplyr::filter(gse37745.expr,rownames(gse37745.expr) %in% common.gene)
gse157010.expr <- as.data.frame(gse157010.expr)
gse157010.expr <- dplyr::filter(gse157010.expr,rownames(gse157010.expr) %in% common.gene)
subtype <- read.table("./Input/group.txt",sep = "\t",header = T)
subtype$sample <- substr(subtype$sample,1,16)
subtype <- column_to_rownames(subtype,var = "sample")
TCGA.expr <- log(TCGA.expr + 1)
#### 差异分析 ####
diff.func <- function(x){
  vt <- as.vector(as.matrix(TCGA.expr[x,]))
  vt.high <- vt[subtype == "1"]
  vt.low <- vt[subtype == "2"]
  wil.res <- wilcox.test(vt.high,vt.low)
  logFC <- mean(vt.high)-mean(vt.low)
  res <- data.frame(gene=rownames(TCGA.expr)[x],logFC=logFC,pvalue=wil.res$p.value)
  return(res)
}
allDiff <- lapply(1:nrow(TCGA.expr),diff.func)
allDiff <- do.call(rbind,allDiff)
allDiff <- allDiff %>% dplyr::mutate(adj.pval=p.adjust(pvalue))
diff <- allDiff %>% dplyr::filter(abs(logFC)>0.5,adj.pval<0.01)
#cox
unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],3)
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
# 单因素分析
tcga.surv <- rownames_to_column(tcga.surv,var = "sample")
unicox.dat <- TCGA.expr %>% 
  dplyr::filter(rownames(.)%in%diff$gene) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  inner_join(tcga.surv[,c("sample","OS","OS.time")],by="sample") %>% 
  dplyr::rename_all(funs(str_replace_all(.,"-","_")))
##准备BaSurv
BaSurv<-Surv(time = unicox.dat$OS.time,event = unicox.dat$OS)
##提取变量名
names(unicox.dat)
varnames<-colnames(unicox.dat)[2:(ncol(unicox.dat)-2)]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(unicox.dat,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
write.csv(univar,file = "./Output/Data/univar.csv",row.names = F)
cox.res <- univar[univar$pvalue <= 0.005,"group"]

library(openxlsx)
library(seqinr)

library(plyr)
library(survival)

library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
Sys.setenv(LANG = "EN")
setwd("./Script/PrognosticML_4.0/")
source("./Codes/ML.R")
## 训练集机器模型数据
Train_surv <- unicox.dat[,c("sample","OS","OS.time")] %>% 
  column_to_rownames(var = "sample")
Train_expr <- unicox.dat[,c("sample",univar$Characteristics[univar$pvalue<=0.005])] %>% 
  column_to_rownames(var = "sample") %>% 
  as.matrix() %>%
  t() %>%
  as.data.frame()
train_os <- "OS"
train_os.time <- "OS.time"
FinalModel <- c("panML", "multiCox")[1]
load("D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/gse37745.Rdata")
load("D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/gse157010.Rdata")
gse37745.cli$OS.time <- as.numeric(gse37745.cli$OS.time)
gse37745.cli <- gse37745.cli %>%
  dplyr::filter(., OS.time >= 30) %>%
  dplyr::mutate(Cohort = "GSE37745")

gse157010.cli <- gse157010.cli %>%
  dplyr::filter(., OS.time >= 30) %>%
  dplyr::mutate(Cohort = "GSE157010")
Test_surv <- rbind(gse157010.cli,gse37745.cli)  
gse37745.expr <- rownames_to_column(gse37745.expr,var = "gene")
gse157010.expr <- as.data.frame(gse157010.expr)
gse157010.expr <- rownames_to_column(gse157010.expr,var = "gene")
Test_expr <- merge(gse157010.expr,gse37745.expr,by = "gene") %>%
  column_to_rownames(var = "gene")

comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]


comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 


Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) 
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_surv$Cohort), function(x) summary(apply(x, 2, var))) # ????scale????

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------

methods <- read.xlsx("./Codes/41467_2022_28421_MOESM4_ESM.xlsx", startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 4
timeVar = "OS.time"; statusVar = "OS" 


## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed = 0419) 
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_expr = Train_set, 
                                 Train_surv = Train_surv, 
                                 mode = "Variable",       
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() 
set.seed(seed = 040419) 
for (method in methods){ 
  # method <- "CoxBoost+plsRcox" 
  cat(match(method, methods), ":", method, "\n") 
  method_name = method 
  method <- strsplit(method, "\\+")[[1]] 
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2], 
                                  Train_expr = Train_expr[, selected.var], 
                                  Train_surv = Train_surv, 
                                  mode = "Model",       
                                  classVar = classVar)  
  }
  
  
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model,"D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/model.rds")



if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ 
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) 
    return(tmp)
  })
}


saveRDS(coxmodel,"D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/coxmodel.rds") 

## Evaluate the model -----------------------------------------------------
model <- readRDS(file.path(res.path, "model.rds")) 
# model <- readRDS(file.path(res.path, "coxmodel.rds")) 

methodsValid <- names(model) 
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = rbind.data.frame(Train_set,Test_set), 
                                    type = "lp") 
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat,"D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/RS_mat.txt",sep = "\t", row.names = T, col.names = NA, quote = F) 
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) 
}

fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) 
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features" 
write.table(fea_df,  "D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/fea_df.txt",sep = "\t", row.names = F, col.names = T, quote = F)

Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], 
                                  Test_expr = Test_set, 
                                  Test_surv = Test_surv, 
                                  Train_expr = Train_set, 
                                  Train_surv = Train_surv, 
                                  Train_name = "TCGA", 
                                  #Train_expr = NULL,
                                  #Train_surv = NULL, 
                                  cohortVar = "Cohort", 
                                  timeVar = timeVar,
                                  statusVar = statusVar) 
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, "D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Output/Data/cindex_mat.txt",sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) 
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] 

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
names(CohortCol) <- colnames(Cindex_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, 
                    avg_Cindex = avg_Cindex, 
                    CohortCol = CohortCol, 
                    barCol = "steelblue", 
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 
setwd("../")
pdf("./Output/Figure/heatmap of cindex.pdf", width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") 
invisible(dev.off())
fina.model <- data.frame(
  gene = model[["RSF+StepCox[forward]"]][["subFeature"]],
  coefficients = model[["StepCox[forward]"]][["coefficients"]]
)
saveRDS(fina.model, file = "./Output/Data/fina.model.Rds")
## 读取泛癌表达谱和注释文件
# 文件较大，占内存，请谨慎操作
panexpr <- fread("./Input/tcga_RSEM_gene_tpm",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
annopb <- read.delim("./Input/gencode.v23.annotation.gene.probemap",row.names = 1,check.names = F,stringsAsFactors = F,header = T,sep = "\t")

# 匹配基因名
panexpr <- as.data.frame(panexpr)
rownames(panexpr) <- panexpr$sample; panexpr <- panexpr[,-1]
comgene <- intersect(rownames(annopb), rownames(panexpr))
panexpr <- panexpr[comgene,]; annopb <- annopb[comgene,]
panexpr$genename <- annopb$gene; panexpr <- panexpr[!duplicated(panexpr$genename),]
rownames(panexpr) <- panexpr$genename; panexpr <- panexpr[,-ncol(panexpr)]
panexpr[1:3,1:3]

## 读取目标基因及其所在的分组信息，用于最后画热图
final.model <- readRDS("./Output/Data/final.model.Rds")
gene <- as.data.frame(final.model$subFeature)
# 提取目标基因的表达矩阵
panexpr <- panexpr[gene$`final.model$subFeature`,]
gc() # 释放内存

# 还原表达谱（以下步骤请根据自己表达谱的情况来）
panexpr <- 2^panexpr - 0.001 # 原始数据为log2(x+0.001)
panexpr[panexpr < 0] <- 0 # 小于0的值拉到0
# 重新对数化
panexpr <- log2(panexpr + 1) 

# 剔除表达量比较奇怪的基因，出自FigureYa35batch_bestSeparationV3_update
panexpr <- panexpr[,apply(panexpr, 1, sd) > 0] # 取方差大于1的基因

## 读取样本注释和生存信息
pansurv <- read.table("./Input/pancancerSurvivalData.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(pansurv) <- paste0(rownames(pansurv),"-01") # 创建原位癌barcode

# 取出例文癌症
tumors <- unique(pansurv$type) # 提取有生存信息的癌种
comsam <- intersect(colnames(panexpr),rownames(pansurv)) # 提取共享的样本
pansurv <- pansurv[comsam,]
panexpr <- panexpr[,comsam]
tumors <- unique(pansurv$type)
gene_group <- setNames(gene,"Symbol")

# cox分析的数据初始化
survland.coxhr <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化cox分析HR结果
survland.coxp <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化cox分析p值结果
survland.coxplot <- matrix(0,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化绘图数据

# logrank分析的数据初始化
survland.logrankhr <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化logrank分析hr值结果
survland.logrankp <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化logrank分析p值结果
survland.logrankplot <- matrix(0,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) # 初始化绘图数据

# 循环计算每一个癌症
for(t in tumors) {
  for (g in gene_group$Symbol) { # 循环计算每一个基因
    sam <- rownames(pansurv[which(pansurv$type == t),]) #提取当前癌症的sample ID
    expr <- as.numeric(panexpr[g,sam]) # 提取当前基因的表达量
    
    expr.surv <- data.frame(futime = pansurv[sam,"OS.time"], # 提取当前癌症的生存信息
                            fustat = pansurv[sam,"OS"], # 提取当前癌症的生存信息
                            expr = expr, # 基因表达量
                            stringsAsFactors = F)
    
    ## 方法一：cox
    cox <- coxph(Surv(futime,fustat) ~ expr, data = expr.surv) # cox分析
    coxSummary <- summary(cox)
    hr <- as.numeric(coxSummary$coefficients[,"exp(coef)"])[1] # 提出HR
    pvalue <- as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1] # 提出p值
    survland.coxhr[g,t] <- hr
    survland.coxp[g,t] <- pvalue
    
    # 为画图准备矩阵
    if(pvalue < 0.05) { # 如果p值显著的话存储数据
      survland.coxplot[g,t] <- ifelse(hr > 1, 1, -1) # HR>1为风险因素，记为“1”，HR<1为保护因素，记为-1，其余默认为0
    }
    
    ## 方法二：logrank
    # 用中值（median）为样本分组，如果想用最佳分组，可参考FigureYa35batch_bestSeparationV3_update
    expr.surv$group = ifelse(expr > median(expr),"high","low")
    expr.surv$group <- factor(expr.surv$group, levels = c("low", "high"))
    
    data.survdiff <- survdiff(Surv(futime,fustat) ~ group, data = expr.surv)
    pvalue <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    hr <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    survland.logrankhr[g,t] <- hr
    survland.logrankp[g,t] <- pvalue
    
    # 为画图准备矩阵
    if(pvalue < 0.05) { # 如果p值显著的话存储数据
      survland.logrankplot[g,t] <- ifelse(hr > 1, 1, -1) # HR>1为风险因素，记为“1”，HR<1为保护因素，记为-1，其余默认为0
    }
  }
}

## 保存到文件，便于DIY其他形式的图
# 或者以更大范围的基因（例如全基因组）作为输入，然后用HR和pvalue筛选出genes associated with the overall survival of patients (worse/better survival) in at least one cancer type或你关心的几种癌症或所有癌症
# cox
write.table(survland.coxplot, file = "cox_genes associated with the OS.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(survland.coxhr,file = "cox HR in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(survland.coxp,file = "cox pvalue in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# logrank
write.table(survland.logrankplot, file = "logrank_genes associated with the OS.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(survland.logrankhr,file = "logrank HR in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(survland.logrankp,file = "logrank pvalue in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


# 自定义颜色
red <- "#D02E20"
blue <- "#4D76B7"
green <- "#50B544"
yellow <- "#F8C77A"
cyan <- "#5AC8F9"

annRow <- gene_group # 行基因注释
rownames(annRow) <- annRow$Symbol
annColors <- list("Function" = c("Readers" = green,
                                 "Writers" = yellow,
                                 "Eraser" = cyan))

# cox
pheatmap1(survland.coxplot,
         border_color = "grey50",
         show_rownames = T, # 显示行名
         show_colnames = T, # 显示列明
         cluster_rows = F, # 行不聚类
         cluster_cols = F, # 列不聚类
         color = c(blue,"grey95",red),
         #annotation_row = annRow[,"Function",drop = F],
         #annotation_colors = annColors,
         legend_breaks = c(-1,0,1), # 修改图例的显示位置
         legend_labels = c("Protective","p>0.05","Risky"), # 修改图例标签
         cellwidth = 10, # 单元格宽度
         cellheight = 10, # 单元格高度
         filename = "survival landscape using cox.pdf", # 保存文件
         width = 8, # 图片宽度
         height = 6) # 图片高度
