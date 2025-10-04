library(dplyr)
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(tibble)
library(ggplot2)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load("./Input/TCGA_Data/output_mRNA_lncRNA_expr/TCGA-LUSC_mrna_expr_fpkm.rdata")
colnames(mrna_expr_fpkm)  <- substr(colnames(mrna_expr_fpkm),1,16)
cmoic <- readRDS("./Output/Data/cmoic.RDS")
mrna_expr_fpkm <- dplyr::select(mrna_expr_fpkm,cmoic$clust.res$samID)
mrna_expr_fpkm <- log(mrna_expr_fpkm + 1)
fnSig <- "D:/Bioinfo/Bioinfo_project/MOVICS_LUSC/Input/pcbc-stemsig.tsv"
w <- read.delim(fnSig, header = FALSE, row.names = 1 ) %>% as.matrix() %>% drop()
w[1:10]
#
mrna_expr_fpkm <- expr
X <- mrna_expr_fpkm %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "gene_id") %>%
  filter( gene_id %in% names(w) ) %>%
  column_to_rownames( "gene_id" ) %>% as.matrix()


# Reduce the signature to the common set of genes.
stopifnot( all( rownames(X) %in% names(w)))
w <- w[ rownames(X) ]
w[1:5]

# Score the Matrix `X` using Spearman correlation. 
s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)
s[1:5]

# Then output scores to a file mRNA_StemScore.tsv.
write.table( cbind(s), file="./MB_mRNA_StemScore.tsv", sep="\t", quote=FALSE, col.names=FALSE )
ss <- data.frame(Sample.ID = names(s), mRNAsi = s)
input <- Group %>% left_join(ss, join_by("ID" == "Sample.ID"))
input <- input[order(input$mRNAsi),]
save(input,file = "./Output/Data/RNAinput.Rdata")


input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$mRNAsi <- as.numeric(input$mRNAsi)
input$index <- 1:nrow(input)

#mRNAsi
p <- t.test(input[which(input$group == "C"),"mRNAsi"],
                 input[which(input$group == "S"),"mRNAsi"])$p.value
pdf("Output/Figure/mRNAsi.pdf",width = 4,height = 6)
# 设置颜色
jco <- c("#DD492E","#40548A")
input$group <- as.factor(input$group)
ggplot(data = input,aes(x = group, y = mRNAsi, fill = group))+
  scale_fill_manual(values = jco) + 
  #geom_violin(alpha=0.4, position = position_dodge(width = .75),
  #            size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = F, #outlier.size = -1, 
               color=c("#DD492E","#40548A"), lwd=0.8, alpha = 0.3)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("mRNAsi Score")) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=1.2, # 根据自己的数据调节p value的位置
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

load("./Input/TCGA_Data/output_methy/TCGA-LUSC_beta_expr_and_pd.rdata")
load("./Input/TCGA_Data/output_methy/TCGA-LUSC_probe_info.rdata")
methy_expr <- beta_expr %>%
  as.data.frame() %>%
  filter(apply(., 1, function(x) sum(is.na(x)) <= 0)) 
colnames(methy_expr) <- substr(colnames(methy_expr),1,16)
methy_expr <- dplyr::select(methy_expr,cmoic$clust.res$samID)
(load("./Input/pcbc-stemsig.p219.Rda"))
w <- mm$w
w[1:5]
w <- w[names(w) %in% rownames(methy_expr)]
X <- methy_expr[as.character(names(w)),]
X <- as.matrix(X)
X[1:5, 1:5]
ss <- t(w) %*% X
ss[1,1:3]

## Scale the scores to be between 0 and 1
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(t(ss))
colnames(ss) <- "mDNAsi"  
head(ss)
save(ss, file = "./Output/Data/MB_mDNAsi.Rda")

ss <- data.frame(Sample.ID = rownames(ss), mDNAsi = ss$mDNAsi)
input <- cmoic$clust.res %>% left_join(ss, join_by("samID" == "Sample.ID"))
input <- input[order(input$mDNAsi),]
p <- wilcox.test(input[which(input$clust == "1"),"mDNAsi"],
                 input[which(input$clust == "2"),"mDNAsi"])$p.value
pdf("Output/Figure/mDNAsi.pdf",width = 4,height = 6)
# 设置颜色
jco <- c("#2EC4B6", "#E71D36")
input$clust <- as.factor(input$clust)
ggplot(data = input,aes(x = clust, y = mDNAsi, fill = clust))+
  scale_fill_manual(values = jco) + 
  #geom_violin(alpha=0.4, position = position_dodge(width = .75),
  #            size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = F, #outlier.size = -1, 
               color=c("#2EC4B6", "#E71D36"), lwd=0.8, alpha = 0.3)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("mRNAsi Score")) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=1, # 根据自己的数据调节p value的位置
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
