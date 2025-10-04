geo_check <- data.table::fread("./Input/GEO/GSE135222_GEO_RNA-seq_omicslab_exp.tsv")
gse135222_cli <- readxl::read_xls("./Input/GEO/GSE135222.CLI.xls")
geo_check$gene_id <- substr(geo_check$gene_id,1,15)
ids <- bitr(geo_check$gene_id,fromType ='ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Hs.eg.db')
geo_check <- geo_check %>%
  dplyr::filter(geo_check$gene_id %in% ids$ENSEMBL) %>%
  left_join(ids, join_by(gene_id == ENSEMBL)) 
geo_check <- aggregate(x = geo_check,by = list(geo_check$SYMBOL), FUN = mean) 
geo_check <- geo_check %>%
  column_to_rownames(var = "Group.1")
geo_check <- geo_check[,2:28]
gse135222 <- log(geo_check + 1)
gse135222_cli$`Sample ID` <- paste0("NSCLC",gse135222_cli$`Sample ID`)
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
cmoic <- readRDS("./Output/Data/cmoic.RDS")
colnames(mrna_expr_fpkm) <- substr(colnames(mrna_expr_fpkm),1,16)
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
gse135222.ntp.pred <- runNTP(expr       = gse135222,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = T, # scale input data (by default)
                       centerFlag = T, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR GSE135222") 
pred.clust <- gse135222.ntp.pred$clust.res %>%
  left_join(gse135222_cli,join_by("samID" == "Sample ID"))
pred.clust$value <- 1
pdf("Output/Figure/02immune_check.pdf",width = 4,height = 6)
ggplot(pred.clust)+
  geom_bar(aes(clust, value, fill = `Clinical benefit`),
           stat = "identity", position = "fill")+
  scale_fill_manual(values = c("#2EC4B6", "#E71D36"),name="Responder")+
  # annotate("text", x = 2, y = 0.85, label="*", size = 5)+
 # ggtitle(str_c("pvalue = ",fishres))+
  xlab("")+
  ylab("Percent")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 18))
dev.off()



