library(tinyarray)
gse <- geo_download("GSE157010")
class(gse)
length(gse)
exp <- gse$exp
exp[1:4,1:4] # 标准的表达矩阵，行是基因，列是样本
ids <- AnnoProbe::idmap('GPL570') # 配合AnnoProbe
## Setting options('download.file.method.GEOquery'='auto')
## Setting options('GEOquery.inmemory.gpl'=FALSE)
gse157010 <- trans_array(exp, ids)
exp1[1:4,1:4]
gse157010_cli <- gse$pd %>%
  dplyr::select(.,"os_mo:ch1","os_event:ch1") %>%
  rename(c("os_mo:ch1" ="futime","os_event:ch1" = "fustat"))
gse157010_cli$futime <- as.numeric(gse157010_cli$futime) * 30
gse157010_cli$fustat <- as.numeric(gse157010_cli$fustat)
marker.up <- runMarker(moic.res      = cmoic,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-LUSC", # MUST be the same of argument in runDEA()
                       dat.path      = OutTable.path, # path of DEA files
                       res.path      = Figure.path, # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 2000, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = mrna_expr_fpkm, # use normalized expression as heatmap input
                       #annCol        = annCol, # sample annotation in heatmap
                       #annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
# run NTP in Yau cohort by using up-regulated biomarkers
gse157010.ntp.pred <- runNTP(expr       = gse157010,
                             templates  = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag  = T, # scale input data (by default)
                             centerFlag = T, # center input data (by default)
                             doPlot     = TRUE, # to generate heatmap
                             fig.name   = "NTP HEATMAP FOR GSE135222") 
surv.gse157010 <- compSurv(moic.res = gse157010.ntp.pred,
                           surv.info = gse157010_cli,
                           convt.time = "m",
                           #surv.cut = c(0,3),
                           surv.median.line = "hv")
load("./Output/Data/gse37745.Rdata")
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
colnames(gse37745.cli) <- c("futime","fustat")
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
gse37745.cli$futime <- as.numeric(gse37745.cli$futime)
gse37745.cli <- gse37745.cli[gse37745.cli$futime > 30,]
marker.up <- runMarker(moic.res      = cmoic,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-LUSC", # MUST be the same of argument in runDEA()
                       dat.path      = OutTable.path, # path of DEA files
                       res.path      = Figure.path, # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 1000, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = mrna_expr_fpkm, # use normalized expression as heatmap input
                       #annCol        = annCol, # sample annotation in heatmap
                       #annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
# run NTP in Yau cohort by using up-regulated biomarkers
gse37745.ntp.pred <- runNTP(expr       = gse37745.expr,
                             templates  = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag  = T, # scale input data (by default)
                             centerFlag = T, # center input data (by default)
                             doPlot     = TRUE, # to generate heatmap
                             fig.name   = "NTP HEATMAP FOR GSE37745") 
surv.gse157010 <- compSurv(moic.res = gse37745.ntp.pred,
                           surv.info = gse37745.cli,
                           convt.time = "m",
                           #surv.cut = c(0,3),
                           surv.median.line = "hv")


