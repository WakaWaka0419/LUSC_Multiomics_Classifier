setwd("~/wsx/projects/jlx/20230902/3-model/")
library(tidyverse)
library(data.table)
library(ggpubr)
library(pheatmap)
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
# library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(timeROC)
library(survival)
library(survminer)
load("../1-data/01hrd.tcga.RData")
load("../1-data/01geo37745.RData")
load("../1-data/01geo157010.RData")
load("../2-multi.omics/02multiOmics.RData")

tcga.fpkm <- tcga.fpkm %>% 
  dplyr::select(hrd.surv$sample) %>% 
  dplyr::filter(rownames(.)%in%rownames(expr37745),
                rownames(.)%in%rownames(expr157010))
all(colnames(tcga.fpkm)==hrd.surv$sample)

#### 差异分析 ####
diff.func <- function(x){
  vt <- as.vector(as.matrix(tcga.fpkm[x,]))
  vt.high <- vt[hrd.surv$hrd.group=="HRD high"]
  vt.low <- vt[hrd.surv$hrd.group=="HRD low"]
  wil.res <- wilcox.test(vt.high,vt.low)
  logFC <- mean(vt.high)-mean(vt.low)
  res <- data.frame(gene=rownames(tcga.fpkm)[x],logFC=logFC,pvalue=wil.res$p.value)
  return(res)
}
allDiff <- lapply(1:nrow(tcga.fpkm),diff.func)
allDiff <- do.call(rbind,allDiff)
allDiff <- allDiff %>% dplyr::mutate(adj.pval=p.adjust(pvalue))
diff <- allDiff %>% dplyr::filter(abs(logFC)>0.5,adj.pval<0.05)
## 火山图
colors <-  c("#e889bd","grey","#8ea0c9")
volca.dat <- allDiff %>% 
  dplyr::mutate(yplot=-log10(adj.pval),
                state=ifelse(logFC>0.5 & adj.pval<0.05,"up",ifelse(logFC<(-0.5) & adj.pval<0.05,"down","unchange")),
                # logFC=ifelse(logFC>1,1,ifelse(logFC<(-1),-1,logFC))
  )
pdf("../figure/03volcano.pdf",height = 6,width = 6,onefile = FALSE)
ggscatter(data=volca.dat,x="logFC",y="yplot",color = "state",alpha = 0.8,
          size = 2,palette = colors,xlab = "logFC",ylab = "-log10(Pvalue)")+
  grids(linetype = "solid")+
  geom_vline(xintercept = c(-0.5,0.5),lty=4,col="black",lwd=1)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=1)+
  theme_classic()+
  theme(text = element_text(size=14),
        legend.position = c(0.2,0.85),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 0.5))
dev.off()
## 热图
anno <- hrd.surv %>% 
  dplyr::select(sample,age,gender,stage,hrd.group) %>% 
  dplyr::arrange(hrd.group) %>% 
  column_to_rownames(var = "sample")
pheat.dat <- as.data.frame(tcga.fpkm) %>% 
  dplyr::filter(rownames(.)%in%diff$gene) %>% 
  dplyr::select(rownames(anno))
all(colnames(pheat.dat)==rownames(anno))
pdf("../figure/03heatmap.pdf",height = 6,width = 8,onefile = FALSE)
pheatmap::pheatmap(as.matrix(pheat.dat),show_rownames = F, show_colnames = F,
                   color = colorRampPalette(colors = c("#8ea0c9", "white","#e889bd"))(100),
                   # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
                   annotation_colors=list(hrd.group=c(`HRD high`="#e889bd", `HRD low`="#8ea0c9")),
                   cluster_rows = T,cluster_cols = F,scale = "row",
                   annotation_col = anno, annotation_names_col = T,
                   use_raster=T,fontsize = 20)
dev.off()
# PCA分析
pca.dat <- as.data.frame(tcga.fpkm) %>% 
  dplyr::filter(rownames(.)%in%diff$gene)
dat.pca <- prcomp(t(pca.dat),center = T,scale. = T)
plot.pca <- predict(dat.pca)[,1:2]
plot.pca <- data.frame(plot.pca) %>%
  rownames_to_column(var = "sample") %>%
  inner_join(hrd.surv[,c("sample","hrd.group")],by="sample")
pdf("../figure/03pca.pdf",height = 6,width = 6,onefile = FALSE)
ggplot(plot.pca,aes(x = PC1, y = PC2, color = hrd.group)) +
  scale_color_manual(values = c("#e889bd","#8ea0c9"))+
  geom_point(size=3) +
  stat_ellipse()+
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = c(0.2,0.85),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 0.5))
dev.off()

#### 单因素cox回归选择变量 ####
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
unicox.dat <- tcga.fpkm %>% 
  dplyr::filter(rownames(.)%in%diff$gene) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  inner_join(hrd.surv[,c("sample","os","os.time")],by="sample") %>% 
  dplyr::rename_all(funs(str_replace_all(.,"-","_")))
##准备BaSurv
BaSurv<-Surv(time = unicox.dat$os.time,event = unicox.dat$os)
##提取变量名
names(unicox.dat)
varnames<-colnames(unicox.dat)[2:(ncol(unicox.dat)-2)]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(unicox.dat,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
write.csv(univar,file = "03univar.csv",row.names = F)

#### 机器学习模型 ####
source("D:/Bioinfo/BRCA- 2/3.model/functions.R")
## 训练集机器模型数据
Train_surv <- unicox.dat[,c("sample","OS","OS.time")] %>% 
  column_to_rownames(var = "sample")
Train_expr <- unicox.dat[,c("sample",univar$Characteristics[univar$pvalue<0.005])] %>% 
  column_to_rownames(var = "sample") %>% 
  as.matrix()
train_os <- "OS"
train_os.time <- "OS.time"
## 数据标准化
Train_set <- scale(Train_expr, center = T, scale = T) 
## 构建筛选特征的算法和参数
sv.methods <- data.frame(tar1.method=c(rep("Enet",11),rep("StepCox",3),"CoxBoost","SuperPC","plsRcox","RSF","GBM"),
                         tar1.args=c(seq(0,1,0.1),"backward","forward","both",rep(NA,5))) %>% 
  dplyr::mutate(tar1.method.name=ifelse(is.na(tar1.args),tar1.method,
                                        ifelse(tar1.method=="Enet",str_c(tar1.method,"[alpha=",tar1.args,"]"),
                                               ifelse(tar1.method=="StepCox",str_c(tar1.method,"[direction=",tar1.args,"]"),tar1.method))),
                tar1.method.name=ifelse(tar1.method.name=="Enet[alpha=1]","Lasso",
                                        ifelse(tar1.method.name=="Enet[alpha=0]","Ridge",tar1.method.name)))
## 筛选特征
preTrain.var <- list() # 用于保存各算法筛选的变量
set.seed(seed = 123456) # 设置建模种子，使得结果可重复
for (x in 1:nrow(sv.methods)){
  method.name <- sv.methods$tar1.method.name[x]
  try({
    preTrain.var[[method.name]] = RunML(method = sv.methods$tar1.method[x], # 变量筛选所需要的机器学习方法
                                        method_param = sv.methods$tar1.args[x], # 变量筛选所需要的参数
                                        Train_expr = Train_set,         # 训练集有潜在预测价值的变量
                                        Train_surv = Train_surv,     # 训练集分类标签
                                        mode = "Variable",             # 运行模式，Variable(筛选变量)和Model(获取模型)
                                        timeVar = train_os.time,
                                        statusVar = train_os)      # 用于训练的分类变量，必须出现在Train_class中
  })
}
preTrain.var <- preTrain.var[sapply(preTrain.var,length)<ncol(Train_set)] # 去除没有筛掉变量的方法
preTrain.var[["simple"]] <- colnames(Train_set) # 新增未经筛选的变量集（以便后续代码撰写）

#### 训练模型 ####
model.train <- list() # 用于保存各模型的所有信息
set.seed(seed = 123456) # 设置建模种子，使得结果可重复
# Train_set_bk = Train_set # RunML有一个函数(plsRglm)无法正常传参，需要对训练集数据进行存档备份
## 算法包括：Lasso、Ridge、Enet、StepCox、survivalSVM、CoxBoost、SuperPC、plsRcox、RSF、GBM
## 构建筛选特征的算法和参数
mt.methods <- data.frame(tar2.method=c(rep("Enet",11),rep("StepCox",3),"survivalSVM","CoxBoost","SuperPC","RSF","GBM"),
                         tar2.args=c(seq(0,1,0.1),"backward","forward","both",rep(NA,5)))
two.stage.method <- data.frame(tar1=names(preTrain.var)) %>% 
  purrr::pmap_df(~mt.methods %>% dplyr::mutate(tar1.method.name=(...))) %>% 
  dplyr::mutate(tar1.method=str_split(tar1.method.name,"\\[",simplify = T)[,1],
                tar2.method.name=ifelse(is.na(tar2.args),tar2.method,
                                        ifelse(tar2.method=="Enet",str_c(tar2.method,"[alpha=",tar2.args,"]"),
                                               ifelse(tar2.method=="StepCox",str_c(tar2.method,"[direction=",tar2.args,"]"),tar2.method))),
                tar2.method.name=ifelse(tar2.method.name=="Enet[alpha=1]","Lasso",
                                        ifelse(tar2.method.name=="Enet[alpha=0]","Ridge",tar2.method.name)),
                method.name=ifelse(tar1.method.name=="simple",tar2.method.name,str_c(tar1.method.name," + ",tar2.method.name))) %>% 
  dplyr::filter(tar1.method!=tar2.method)
## 构建各种模型组合
model.train <- list()
for (x in 12:nrow(two.stage.method)){
  cat(x, " Run method :", two.stage.method$method.name[x], "\n")
  method.name = two.stage.method$method.name[x] # 本轮算法名称
  
  Variable = preTrain.var[[two.stage.method$tar1.method.name[x]]] # 根据方法名称的第一个值，调用先前变量筛选的结果
  Train_set_2 = Train_set[, Variable]   # 对训练集取子集，因为有一个算法原作者写的有点问题，无法正常传参
  try({
    model.train[[method.name]] <- RunML(method = two.stage.method$tar2.method[x],     # 根据方法名称第二个值，调用构建的函数分类模型
                                        method_param = two.stage.method$tar2.args[x], # 变量筛选所需要的参数
                                        Train_expr = Train_set_2,         # 训练集有潜在预测价值的变量
                                        Train_surv = Train_surv,     # 训练集分类标签
                                        mode = "Model",             # 运行模式，Variable(筛选变量)和Model(获取模型)
                                        timeVar = train_os.time,
                                        statusVar = train_os)      
  })
  # 如果最终模型纳入的变量数小于预先设定的下限，则认为该算法输出的结果是无意义的
  # if(length(ExtractVar(model.train[[method.name]])) <= min.selected.var) {
  #   model.train[[method.name]] <- NULL
  # }
}
save(model.train, file = "./Output/Data/TIMEmodel.train.RData") 
## 计算预测和cindex
pred.cindex.train.res <- lapply(model.train,function(x){CalRiskScore(fit = x,new_data = Train_set,new_surv = Train_surv,
                                                                     timeVar = train_os.time, statusVar = train_os, time.unit = "day")})
pred.train.res <- lapply(pred.cindex.train.res,function(x){x$predict.res})
cindex.train.res <- lapply(pred.cindex.train.res,function(x){x$cindex})
auc.train.res <- lapply(pred.cindex.train.res,function(x){x$auc})
feature.select.res <- lapply(model.train, function(fit){ExtractVar(fit)})
## 整理数据
auc.train.df <- lapply(auc.train.res, function(x)(x$AUC[c(2,3,5)]))
model.train.pd <- docall(rbind,auc.train.df)
model.train.pd <- as.data.frame(model.train.pd) %>% 
  dplyr::mutate(cindex=unlist(cindex.train.res))
colnames(model.train.pd) <- c("tyear3","tyear5","tyear10","tcindex")

#### 验证集GEO37745 ####
Test37745_surv <- gse37745.cli[,c("OS","OS.time")] %>% 
  dplyr::filter(!is.na(OS),!is.na(OS.time))
rownames(gse37745.expr) <- str_replace_all(rownames(gse37745.expr),"-","_")
Test37745_expr <- gse37745.expr %>%
  dplyr::filter(rownames(.)%in%univar$Characteristics[univar$pvalue<0.005]) %>% 
  dplyr::select(rownames(Test37745_surv)) %>% 
  t() 
Test37745_os <- "OS"
Test37745_os.time <- "OS.time"
## 数据标准化
Test37745_set <- scale(Test37745_expr, center = T, scale = T) 
## 预测结果
Test37745_surv$OS.time <- as.numeric(Test37745_surv$OS.time)
pred.cindex.Test37745.res <- lapply(model.train,function(x){CalRiskScore(fit = x,new_data = Test37745_set,new_surv = Test37745_surv,
                                                                         timeVar = Test37745_os.time, statusVar = Test37745_os, time.unit = "day")})
pred.Test37745.res <- lapply(pred.cindex.Test37745.res,function(x){x$predict.res})
cindex.Test37745.res <- lapply(pred.cindex.Test37745.res,function(x){x$cindex})
auc.Test37745.res <- lapply(pred.cindex.Test37745.res,function(x){x$auc})
## 整理数据
auc.v37745.df <- lapply(auc.Test37745.res, function(x)(x$AUC[c(2,3,5)]))
model.v37745.pd <- docall(rbind,auc.v37745.df)
model.v37745.pd <- as.data.frame(model.v37745.pd) %>% 
  dplyr::mutate(cindex=unlist(cindex.Test37745.res))
colnames(model.v37745.pd) <- c("v37745year3","v37745year5","v37745year10","v37745cindex")

#### 验证集GEO157010 ####
Test157010_surv <- gse157010.cli[,c("OS","OS.time")] %>% 
  dplyr::filter(!is.na(OS),!is.na(OS.time))
rownames(gse157010.expr) <- str_replace_all(rownames(gse157010.expr),"-","_")
Test157010_expr <- gse157010.expr %>%
  dplyr::filter(rownames(.)%in%univar$Characteristics[univar$pvalue<0.005]) %>% 
  dplyr::select(rownames(Test157010_surv)) %>% 
  t() 
Test157010_os <- "OS"
Test157010_os.time <- "OS.time"
## 数据标准化
Test157010_set <- scale(Test157010_expr, center = T, scale = T) 
## 预测结果
pred.cindex.Test157010.res <- lapply(model.train,function(x){CalRiskScore(fit = x,new_data = Test157010_set,new_surv = Test157010_surv,
                                                                          timeVar = Test157010_os.time, statusVar = Test157010_os, time.unit = "month")})
pred.Test157010.res <- lapply(pred.cindex.Test157010.res,function(x){x$predict.res})
cindex.Test157010.res <- lapply(pred.cindex.Test157010.res,function(x){x$cindex})
auc.Test157010.res <- lapply(pred.cindex.Test157010.res,function(x){x$auc})
## 整理数据
auc.v157010.df <- lapply(auc.Test157010.res, function(x)(x$AUC[c(2,3,5)]))
model.v157010.pd <- docall(rbind,auc.v157010.df)
model.v157010.pd <- as.data.frame(model.v157010.pd) %>% 
  dplyr::mutate(cindex=unlist(cindex.Test157010.res))
colnames(model.v157010.pd) <- c("v157010year3","v157010year5","v157010year10","v157010cindex")

## 整理数据，绘制热图
model.pheat.data <- Reduce(cbind,list(model.train.pd,model.v37745.pd,model.v157010.pd))
model.pheat.data <- model.pheat.data %>% 
  dplyr::mutate(year1=purrr::pmap_dbl(across(ends_with("year3")),~mean(c(...))),
                year3=purrr::pmap_dbl(across(ends_with("year5")),~mean(c(...))),
                year5=purrr::pmap_dbl(across(ends_with("year10")),~mean(c(...))),
                cindex=purrr::pmap_dbl(across(ends_with("cindex")),~mean(c(...))),
                nVars=sapply(feature.select.res,function(x){length(x)})) %>% 
  dplyr::arrange(desc(cindex))
model.anno <- data.frame(dataSet=c(rep("TCGA",4),rep("GSE37745",4),rep("GSE157010",4),rep("combine",4),"nVars"),
                         type=c(rep(c("3 year","5 year","10 year","c-index"),4),"nVars"))
rownames(model.anno) <- colnames(model.pheat.data)
all(colnames(model.pheat.data)==rownames(model.anno))
model.pheat.data1 <- model.pheat.data[1:50,1:16]
annoRow <- model.pheat.data %>% 
  dplyr::select(17) %>% 
  dplyr::slice(1:50)
pdf("../figure/03heatmap.model.1.pdf",height = 12,width = 8,onefile = FALSE)
pheatmap::pheatmap(as.matrix(model.pheat.data1),show_rownames = T, show_colnames = F,
                   color = colorRampPalette(colors = c("#8ea0c9", "white","#e889bd"))(100),
                   # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
                   # annotation_colors=list(hrd.group=c(`HRD high`="#e889bd", `HRD low`="#8ea0c9")),
                   cluster_rows = F,cluster_cols = F,scale = "none",
                   annotation_col = model.anno, annotation_row = annoRow,
                   annotation_names_col = T,
                   gaps_col =  c(12),
                   use_raster=T,fontsize = 10)
dev.off()
model.pheat.data1 <- model.pheat.data[51:nrow(model.pheat.data),1:16]
annoRow <- model.pheat.data %>% 
  dplyr::select(17) %>% 
  dplyr::slice(51:nrow(model.pheat.data))
pdf("../figure/03heatmap.model.2.pdf",height = 26,width = 8,onefile = FALSE)
pheatmap::pheatmap(as.matrix(model.pheat.data1),show_rownames = T, show_colnames = F,
                   color = colorRampPalette(colors = c("#8ea0c9", "white","#e889bd"))(100),
                   # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
                   # annotation_colors=list(hrd.group=c(`HRD high`="#e889bd", `HRD low`="#8ea0c9")),
                   cluster_rows = F,cluster_cols = F,scale = "none",
                   annotation_col = model.anno,  annotation_row = annoRow,
                   annotation_names_col = T,
                   gaps_col =  c(12),
                   use_raster=T,fontsize = 10)
dev.off()

#### 绘制生存图和ROC曲线 ####
train.surv <- Train_surv %>% dplyr::mutate(pred=pred.train.res[[rownames(model.pheat.data)[1]]])
test37745.surv <- Test37745_surv %>% dplyr::mutate(pred=pred.Test37745.res[[rownames(model.pheat.data)[1]]])
test157010.surv <- Test157010_surv %>% 
  dplyr::mutate(pred=pred.Test157010.res[[rownames(model.pheat.data)[1]]])
## 生存曲线
km.title <- c("TCGA","GSE37745","GSE157010")
km.list <- list(train.surv,test37745.surv,test157010.surv)
km.res <- list()
for (x in 1:length(km.list)) {
  cutPoint <- surv_cutpoint(km.list[[x]],time = "OS.time",event = "OS",variables = c("pred"))
  km.list[[x]] <- km.list[[x]] %>% dplyr::mutate(group=ifelse(pred>summary(cutPoint)[1,1],"high","low"))
  fit <- survfit(Surv(OS.time,OS)~group,data=km.list[[x]])
  p <- ggsurvplot(fit, data = km.list[[x]],
                  conf.int = T,
                  pval = TRUE,
                  pval.coord = c(1,1),
                  fun = "pct",
                  risk.table = T,
                  # risk.table.height = 0.2,
                  palette = c("#e889bd","#8ea0c9"),
                  size = 1,
                  pval.size = 6,
                  font.legend = 12,
                  linetype = "strata",
                  risk.table.title = "",
                  title=km.title[x],
                  ##调色板
                  legend = c(0.14,0.22),
                  legend.labs = c("high score","low score"),
                  legend.title = "group",
                  ggtheme = theme_classic()+
                    theme(legend.background = element_rect(fill = "transparent",colour = NA),
                          plot.title = element_text(size=25,hjust=0.5),
                          legend.text=element_text(size=15),
                          panel.border = element_rect(color = "black",fill = NA,size = 0.5))
  ) 
  km.res <- rlist::list.append(km.res,p)
}
pdf("./Output/Figure/surv.plot.pdf",height = 6,width = 15,onefile = FALSE)
arrange_ggsurvplots(km.res, print = TRUE, ncol = 3, nrow = 1)
dev.off()
## auc曲线
km.list[[3]][["OS.time"]] <- as.integer(km.list[[3]][["OS.time"]])
for (x in 1:length(km.list)) {
  roc.res <- timeROC(T=km.list[[x]][["OS.time"]],delta=km.list[[x]][["OS"]],marker=km.list[[x]][["pred"]],
                     cause=1,weighting="marginal",
                     times=seq(0,3650,365))
  pdf(str_c("./Output/Figure/03AUC.",km.title[x],".pdf"),height = 6,width = 6)
  plotAUCcurve(roc.res, FP = 2, add = F, conf.int = F, conf.band = F, col = "#8ea0c9")
  dev.off()
}
