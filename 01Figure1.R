Sys.setenv(LANG = "en")
library(MOVICS)
library(tibble)
library(dplyr)
#路径设置
rm(list = ls())
getwd()
Figure.path <- "./Figure/"
OutData.path <- "./Data/"
OutTable.path <- "./Table/"
InData.path <- "Input/"
#合并多组学数据
# 使用 load() 函数加载数据文件
mo.data <- readRDS("Data/mo.data.Rds")
#确定最佳亚型
optk.GBM <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = file.path(Figure.path,"01GBM_Cluster_Number.pdf"))
#单一分型结果
iClusterBayes.res <- getMOIC(data        = mo.data,
                             N.clust     = 2,
                             methodslist = "iClusterBayes", # 指定算法
                             type        = c("gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
                             n.burnin    = 1800,
                             n.draw      = 1200,
                             prior.gamma = c(0.5, 0.5, 0.5,0.5),
                             sdev        = 0.05,
                             thin        = 3)
#9种分型
# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering","IntNMF", "CIMLR", "MoCluster"), # 9种算法
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian",  "gaussian","binomial"))
#合并分型结果
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
# 保存下结果
saveRDS(moic.res.list, file = "./Data/moic.res.list.Rds")
#分型结果可视化
cmoic <- getConsensusMOIC(moic.res.list = moic.res.list,
                          fig.name      = "Figure/01CONSENSUS HEATMAP",
                          distance      = "euclidean",
                          linkage       = "average")
saveRDS(cmoic,file = "./Data/cmoic.RDS")
#检查分型质量
getSilhouette(sil      = cmoic$sil, # a sil object returned by getConsensusMOIC()
              fig.path = "./Figure/",
              fig.name = "01SILHOUETTE",
              height   = 5.5,
              width    = 5)

# β值矩阵转换为M值矩阵、
indata <- mo.data
indata[[4]] <- log2(indata[[4]] / (1 - indata[[4]]))

# 对数据进行标准化
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation
#我们这里就用贝叶斯分型的结果进行展示，首先是提取每个组学的结果，然后每个组学中选择前10个分子进行标注
feat   <- moic.res.list$iClusterBayes$feat.res
feat1  <- feat[which(feat$dataset == "dat1"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "dat2"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "dat3"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "dat4"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "dat5"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)

## 为每个组学的热图自定义颜色，不定义也可
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mirna.col <-  c("yellow","blue")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, mirna.col,meth.col, mut.col)
## extract PAM50, pathologic stage and age for sample annotation
cli_name <- colnames(mo.data[[1]])
cli.data <- data.table::fread("Input/TCGA-GBM.GDC_phenotype.tsv") %>%
  dplyr::rename("sample" = "submitter_id.samples") %>%
  dplyr::filter(sample %in% cli_name) %>%
  column_to_rownames(var = "sample") 

annCol    <- cli.data[,c("tumor_stage.diagnoses", "gender.demographic", "person_neoplasm_cancer_status",
                         "age_at_initial_pathologic_diagnosis","race.demographic","vital_status.demographic",
                         "site_of_resection_or_biopsy.diagnoses","synchronous_malignancy.diagnoses"), drop = FALSE]
annCol    <- annCol %>%
  mutate(person_neoplasm_cancer_status  = ifelse(person_neoplasm_cancer_status == "","not reported",person_neoplasm_cancer_status)) %>%
  rename("Stage" = "tumor_stage.diagnoses",
         "Gender" = "gender.demographic",
         "New Tumor Event" = "person_neoplasm_cancer_status",
         "Age" = "age_at_initial_pathologic_diagnosis",
         "Race" = "race.demographic",
         "OS" = "vital_status.demographic",
         "Tumor Location" = "site_of_resection_or_biopsy.diagnoses",
         "Malignancy Diagnoses" = "synchronous_malignancy.diagnoses") %>% 
  mutate(Stage = case_when(
    Stage %in% "stage i" ~ "Stage I",
    Stage %in% "stage ia" ~ "Stage I",
    Stage %in% "stage ib" ~ "Stage I",
    Stage %in% "stage ii" ~ "Stage II",
    Stage %in% "stage iia" ~ "Stage II",
    Stage %in% "stage iib" ~ "Stage II",
    Stage %in% "stage iii" ~ "Stage III",
    Stage %in% "stage iiia" ~ "Stage III",
    Stage %in% "stage iiib" ~ "Stage III",
    Stage %in% "stage iv" ~ "Stage IV"
  )) %>%
  mutate(Age = ifelse(is.na(Age),50,Age)) %>%
  select(.,Age,Stage,Gender,`Tumor Location`)

# generate corresponding colors for sample annotation
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age),
                                                           median(annCol$Age), 
                                                           max(annCol$Age)), 
                                                colors = c("#BCC6DD", "#B5CBE2","#455D99")),
                  Stage  = c("Stage I" = "#8ab1d2",
                             "Stage II"   = "#E58579",
                             "Stage III"   = "#D9BDD8",
                             "Stage IV"   = "#9180AC",
                             "not reported" = "#999999"),
                  Gender  = c("female" = "#42BCB2",
                              "male"   = "#235689"),
                  "New Tumor Event" = c("WITH TUMOR"    = "#A8817A",
                                        "TUMOR FREE"    = "#E8BE74",
                                        "not reported"    = "#999999"),
                  "Tumor Location" = c(
                    "Lower lobe, lung" = "#BC3C2999",
                    "Lung, NOS" = "#0072B599",
                    "Main bronchus" = "#E1872799" ,
                    "Middle lobe, lung" = "#20854E99" ,
                    "Overlapping lesion of lung" = "#7876b199" ,
                    "Upper lobe, lung" = "#6F99AD99" 
                  )
)
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","miRNA value","M value","Mutated"),
             clust.res     = cmoic$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = annCol, # no annotation for samples
             annColors     = annColors, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "Output/Figure/01COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

# survival comparison
select.tumor <- colnames(mo.data[[1]])
surv <- read.table("Input/TCGA-GBM.survival.tsv", sep = "\t", header = T)
surv <- surv[,c(1,2,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% select.tumor) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::rename("fustat" = "OS","futime" = "OS.time")

surv.compare <- compSurv(moic.res         = cmoic,
                      surv.info        = surv,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(0,1), # estimate 5 and 10-year survival
                      fig.name         = "Output/Figure/01KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

# mutational frequency comparison
mut_tcga <- compMut(moic.res     = cmoic,
                    mut.matrix   = mo.data[[5]], # 0/1矩阵
                    doWord       = TRUE, # 生成Word文档
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.2, # 保留在至少5%的样本中突变的基因
                    p.cutoff = 0.05, # 保留padj<0.05的基因
                    p.adj.cutoff = 0.05,
                    innerclust   = TRUE, # 在每个亚型中进行聚类
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 8, 
                    height       = 4,
                    fig.path = Figure.path,
                    res.path = OutTable.path,
                    fig.name     = "01ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
#compare cli
## extract PAM50, pathologic stage and age for sample annotation
cli_name <- colnames(mo.data[[1]])
cli.data <- data.table::fread("Input/TCGA-GBM.GDC_phenotype.tsv") %>%
  dplyr::rename("sample" = "submitter_id.samples") %>%
  dplyr::filter(sample %in% cli_name) %>%
  column_to_rownames(var = "sample") 
select.tumor <- colnames(mo.data[[1]])
surv <- read.table("Input/TCGA-GBM.survival.tsv", sep = "\t", header = T)
surv <- surv[,c(1,4)] %>%
  #mutate(sample = sample.tumor) %>%
  filter(sample %in% select.tumor) 

annCol    <- cli.data[,c("tumor_stage.diagnoses", "gender.demographic", "person_neoplasm_cancer_status",
                         "age_at_initial_pathologic_diagnosis","race.demographic","vital_status.demographic",
                         "site_of_resection_or_biopsy.diagnoses","synchronous_malignancy.diagnoses"), drop = FALSE]
Compare_cli <- annCol %>%
  mutate(person_neoplasm_cancer_status  = ifelse(person_neoplasm_cancer_status == "","not reported",person_neoplasm_cancer_status)) %>%
  rename("Stage" = "tumor_stage.diagnoses",
         "Gender" = "gender.demographic",
         "New Tumor Event" = "person_neoplasm_cancer_status",
         "Age" = "age_at_initial_pathologic_diagnosis",
         "Race" = "race.demographic",
         "OS" = "vital_status.demographic",
         "Tumor Location" = "site_of_resection_or_biopsy.diagnoses",
         "Malignancy Diagnoses" = "synchronous_malignancy.diagnoses") %>% 
  mutate(Stage = case_when(
    Stage %in% "stage i" ~ "Stage I",
    Stage %in% "stage ia" ~ "Stage I",
    Stage %in% "stage ib" ~ "Stage I",
    Stage %in% "stage ii" ~ "Stage II",
    Stage %in% "stage iia" ~ "Stage II",
    Stage %in% "stage iib" ~ "Stage II",
    Stage %in% "stage iii" ~ "Stage III",
    Stage %in% "stage iiia" ~ "Stage III",
    Stage %in% "stage iiib" ~ "Stage III",
    Stage %in% "stage iv" ~ "Stage IV"
  )) %>%
  rownames_to_column(var = "sample") %>%
  left_join(cmoic$clust.res,join_by("sample" == "samID")) %>%
  left_join(surv,by = "sample") %>%
  column_to_rownames(var = "sample")

compClinvar1 <- edit(compClinvar)
clin.GBM <- compClinvar1(moic.res      = cmoic,
                         var2comp      = Compare_cli[,-9], # data.frame needs to summarize (must has row names of samples)
                         #strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("Gender","OS",
                                           "New Tumor Event","Race",
                                           "Tumor Location","Malignancy Diagnoses"), # features that are considered categorical variables
                         nonnormalVars = c("OS.time","Age"), # feature(s) that are considered using nonparametric test
                         exactVars     = "Stage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")

##cluster_compare
moic.res.list <- readRDS("./Output/Data/moic.res.list.Rds")
SNF <- moic.res.list$SNF$clust.res
PinsPlus <- moic.res.list$PINSPlus$clust.res
NEMO <- moic.res.list$NEMO$clust.res
COCA <- moic.res.list$COCA$clust.res
LRAcluster <- moic.res.list$LRAcluster$clust.res
ConsensusCluster <- moic.res.list$ConsensusClustering$clust.res
IntNMF <- moic.res.list$IntNMF$clust.res
CIMLR <- moic.res.list$CIMLR$clust.res
MoCluster <- moic.res.list$MoCluster$clust.res
iClusterBayes <- moic.res.list$iClusterBayes$clust.res
cmoic <- readRDS("./Output/Data/cmoic.RDS")
cluster.compare <- data.frame(
  row.names = rownames(iClusterBayes),
  Cmoic = cmoic$clust.res$clust,
  SNF = SNF$clust,
  PinsPlus = PinsPlus$clust,
  NEMO = NEMO$clust,
  COCA = COCA$clust,
  LRAcluster = LRAcluster$clust,
  ConsensusCluster = ConsensusCluster$clust,
  IntNMF = IntNMF$clust,
  CIMLR = CIMLR$clust,
  MoCluster = MoCluster$clust,
  iClusterBayes = iClusterBayes$clust
)
annColors <- list()
annColors[["Cmoic"]] <- c("1" = "#2EC4B6", "2" = "#E71D36")
annColors[["SNF"]] <- c("1" = "#DFEBAF","2" = "#2CA8E1")
annColors[["PinsPlus"]] <- c("1" = "#E80035","2" = "#DBDEDD")
annColors[["NEMO"]] <- c("1" = "#E8536B","2" = "#F6B879")
annColors[["COCA"]] <- c("1" = "#A2D9F1","2" = "#BAC8E5")
annColors[["LRAcluster"]] <- c("1" = "#F2A1A2","2" = "#68BE8B")
annColors[["ConsensusCluster"]] <- c("1" = "#A8BFDE","2" = "#DA8A88")
annColors[["IntNMF"]] <- c("1" = "#F9C8C4","2" = "#FCD179")
annColors[["CIMLR"]] <- c("1" = "#B3BC95", "2" = "#ADBED6")
annColors[["MoCluster"]] <- c("1" = "#AFB6D2", "2" = "#E4DBE8")
annColors[["iClusterBayes"]] <- c("1" = "#EF767A", "2" = "#48C0AA")

load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-GBM_mrna_expr_tpm.rdata")
colnames(mrna_expr_tpm) <- substr(colnames(mrna_expr_tpm),1,16) 
indata <- log2(mrna_expr_tpm[1:20,rownames(cmoic$clust.res)] + 1)
annCol <-arrange(cluster.compare,desc(Cmoic))
new_order <- match(rownames(annCol), colnames(indata))
plot.data <- indata[, new_order]
hm1 <- pheatmap(standarize.fun(plot.data,halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = annCol,
                annotation_colors = annColors,
                #color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 0.6, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类
