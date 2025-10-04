#subgroup_cnv
load("./Input/TCGA_Data/output_cnv/TCGA-LUSC_CNV.rdata")
data_CNV <- data
dim(data_CNV)
cnv <- data_CNV[,-1]
cnv <- cnv[,c(6,1:5)]
##提取肿瘤样品
tumor_seg <- cnv[substr(cnv$Sample,14,15)=="01",]
names(tumor_seg) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
#tumor_seg$Sample <- substr(tumor_seg$Sample,1,16)
cmoic <- readRDS("./Output/Data/cmoic.RDS")
C1.sample <- cmoic$clust.res[cmoic$clust.res$clust == 1 ,"samID"]
C2.sample <- cmoic$clust.res[cmoic$clust.res$clust == 2,"samID"]
tumor_seg <- tumor_seg[substr(tumor_seg$Sample,1,16) %in% union(C1.sample,C2.sample),]
C1.seg.cnv <- tumor_seg[substr(tumor_seg$Sample,1,16) %in% C1.sample,]
C2.seg.cnv <- tumor_seg[substr(tumor_seg$Sample,1,16) %in% C2.sample,]
write.table(C1.seg.cnv, file="./Input/cnv_gastic/tumor.C1.seg.txt", sep="\t", row.names=F, quote = F)
write.table(C2.seg.cnv, file="./Input/cnv_gastic/tumor.C2.seg.txt", sep="\t", row.names=F, quote = F)
write.table(tumor_seg, file="./Input/cnv_gastic/tumor.seg.txt", sep="\t", row.names=F, quote = F)


# Create a chromosomes reference objects function
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg38"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#str(chrom)




pdf("./Output/Figure/LUSC_copy_number_gistic_score.pdf",12,5)
# Import gistic2 results read gistic output file
scores <- read.table("./Input/cnv_gastic/561299_all/scores.gistic", sep="\t",header=T,stringsAsFactors = F)
head(scores)
unique(scores$Chromosome)
#把染色体名从阿拉伯数字改为“chr1”、“chrX”的形式
scores[scores$Chromosome==23, "Chromosome"] <- "X"
scores[scores$Chromosome==24, "Chromosome"] <- "Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))

# Important step for accurate length to match back to continual chrom loci
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score) - 0.1, max(scores$G.score) + 0.1)
title <- paste0("TCGA LUSC overall copy number gistic score", " ", "n=", 341)

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.15,ylim[2]-0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
dev.off()


pdf("./Output/Figure/LUSC_copy_number_percentage.pdf",12,7)
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency<- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency<- scores.del$frequency * -100

# copy number percentage plot
 seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim<- c(-90,90)
title=paste0("LUSC overall copy number percentage"," ","n=",341)

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gain/loss percentage in cohort", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .5))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(-80,80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
dev.off()


pdf("./Output/Figure/CLUSTER_cnv.scores.gistic.pdf",15,12)
par(mfrow=c(2,1), mar = par()$mar + c(3,0,0,3))


### c1 ###
scores <- read.table("./Input/cnv_gastic/561300_c1/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("CMOIC1 copy number gistic score"," ","n=",171)

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### c1 ###
scores <- read.table("./Input/cnv_gastic/561301_c2/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("CMOIC2 copy number gistic score"," ","n=",170)

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
dev.off()


pdf("CLUSTER_cnv.frequence.pdf",15,12)
par(mfrow=c(2,1), mar = par()$mar + c(3,0,0,3))

### ESCC ###
scores <- read.table("./Input/cnv_gastic/561300_c1/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("ESCC, n=",171)

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


### EAC ###
scores <- read.table("./Input/cnv_gastic/561301_c2/scores.gistic", sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency <- scores.del$frequency * -100
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$frequency)-0.1,max(scores$frequency)+0.1)
title=paste0("CMOIC2, n=",170)

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "Frequency", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+10,ylim[2]-10)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


dev.off()




# 读取GISTIC结果
# arm-level
armCNV <- read.delim("./Input/cnv_gastic/561299_all/broad_values_by_arm.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# focal-level
geneCNV <- read.delim("./Input/cnv_gastic/561299_all/all_data_by_genes.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 读取样本分组信息
C1.sample <- colnames(armCNV)[substr(colnames(armCNV),1,16)%in% C1.sample]
C2.sample <- colnames(armCNV)[substr(colnames(armCNV),1,16)%in% C2.sample]
C1.sub <- data.frame(
  "sample" = C1.sample,
  "cluster"= rep("C1",length(C1.sample))
)

C2.sub <- data.frame(
  "sample" = C2.sample,
  "cluster"= rep("C2",length(C2.sample))
)
  
subt <- rbind(C1.sub,C2.sub)
# 设置扩增和缺失的阈值
cna.cutoff <- 0 # 这里认为>0则为扩增，<0则为缺失

# 循环计算不同组别的扩增和缺失比例
outTab.arm <- NULL
for (i in rownames(armCNV)) {
  message(i)
  
  # 取出当前臂的结果
  tmp <- data.frame(cna = as.numeric(armCNV[i,]),
                    group = subt$cluster,
                    stringsAsFactors = F)
  
  # 初始化loss和amp
  loss <- amp <- tmp
  loss$cna <- ifelse(loss$cna < cna.cutoff,"LOSS","Others") # 当cna小于阈值则loss，正常或扩增认为是others
  amp$cna <- ifelse(amp$cna > -cna.cutoff,"AMP","Others") # 当cna大于阈值则amplification，正常或缺失认为是others
  
  # 构建数据框
  loss.dt <- as.data.frame.array(table(loss$cna,loss$group)) 
  amp.dt <- as.data.frame.array(table(amp$cna,amp$group))
  
  # 计算扩增或者缺失的比例，以及对应在两组间的p值
  if(!is.element("AMP",rownames(amp.dt))) { # 如果在该臂中不存在扩增，则比例记为0，p值记为空
    amp.pct <- c(0,0)
    amp.p <- NA
  } else {
    amp.pct <- as.numeric(amp.dt[1,]/colSums(amp.dt))
    amp.p <- fisher.test(amp.dt)$p.value
  }
  
  if(!is.element("LOSS",rownames(loss.dt))) { # 如果在该臂中不存在缺失，则比例记为0，p值记为空
    loss.pct <- c(0,0)
    loss.p <- NA
  } else {
    loss.pct <- as.numeric(loss.dt[1,]/colSums(loss.dt))
    loss.p <- fisher.test(loss.dt)$p.value
  }
  
  outTab.arm <- rbind.data.frame(outTab.arm,
                                 data.frame(arm = i,
                                            loss.C1 = loss.pct[1],
                                            loss.C2 = loss.pct[2],
                                            p.loss = loss.p,
                                            amp.C1 = amp.pct[1],
                                            amp.C2 = amp.pct[2],
                                            p.amp = amp.p,
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
}

# 输出到文件
write.table(outTab.arm, file = "./Output/Data/output_percentage and fisher test of cnv in two groups.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(reshape)
library(ggplot2)
library(patchwork)
wide <- outTab.arm[,1:4]
df <- as.data.frame(melt(wide, measure.vars = c("loss.C1", "loss.C2")))
df$variable <- gsub("loss.","",df$variable,fixed = T)
df$arm <- factor(df$arm, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22p","22q"))

p_bot <- ggplot(df, aes(arm, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
  xlab(NULL) + ylab("Deletion frequency") +
  
  scale_fill_manual(values = c("#2EC4B6", "#E71D36")) + # 可以用这行修改配色
  annotate(geom="text", # 标星*
           x = 1:nrow(wide),
           y = ifelse(wide[,2] > wide[,3],
                      wide[,2] + 0.03, # 微调*的位置
                      wide[,3] + 0.03),
           size = 10, angle = 0, fontface = "bold",
           label = ifelse(wide[,4] < 0.05,"*",""),
           color = "red") + # 星的颜色
  
  theme_bw() + 
  theme(axis.ticks.y = element_line(size = 0.2),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12, color = "black"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.position = "none",
        #panel.border = element_blank(),
        legend.title = element_blank()) +
  scale_y_reverse()

p_bot

wide <- outTab.arm[,c(1,5:7)]
df <- as.data.frame(melt(wide, measure.vars = c("amp.C1", "amp.C2")))
df$variable <- gsub("amp.","",df$variable,fixed = T)
df$arm <- factor(df$arm, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22p","22q"))

p_top <- ggplot(df, aes(arm, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
  xlab(NULL) + ylab("Amplification frequency") +
  scale_fill_manual(values = c("#2EC4B6", "#E71D36")) + # 可以用这行修改配色 
  annotate(geom="text", # 标星*
           x = 1:nrow(wide),
           y = ifelse(wide[,2] > wide[,3],
                      wide[,2] + 0.01, # 微调*的位置
                      wide[,3] + 0.01),
           size = 10, angle = 0, fontface = "bold",
           label = ifelse(wide[,4] < 0.05,"*",""),
           color = "red") +
  
  theme_bw() + 
  theme(axis.ticks.y = element_line(size = 0.2),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.position = "top",
        #panel.border = element_blank(),
        legend.title = element_blank())

p_top


pal <- p_top + p_bot +
  plot_layout(ncol = 1,widths = c(6,6), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.4, 'cm'))
pal

# 输出pdf文件
ggsave("./Output/Figure/arm-level scnv frequency barplot.pdf", width = 12, height = 10)

focal.level <- factor(geneCNV$Cytoband, levels = unique(geneCNV$Cytoband))
focal.cna <- apply(geneCNV[,setdiff(colnames(geneCNV), c("Cytoband","Gene ID"))], 2, 
                   function(x) tapply(x, INDEX=factor(geneCNV$Cytoband), FUN=mean, na.rm=TRUE)) 
focal.cna <- as.data.frame(focal.cna)
focal.cna <- focal.cna[-1,]
focal.cna <- focal.cna[levels(focal.level),] # 按照染色体臂顺序排列

# 保存到文件
write.table(focal.cna, file = "./Output/Data/output_focal-level cnv.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

tmp1 <- as.data.frame(focal.cna[,C1.sample]) # 取出类1的SCNA
tmp1$scna <- as.numeric(rowSums(tmp1)) # 对SCNA取行和，表示该亚型的在染色体上的SCNA水平
tmp1$loci <- 1:nrow(tmp1)
lw1 <- loess(scna ~ loci,data=tmp1,span = 0.02) # 直接使用行和会让折线变得很尖锐，所以这里我采用loess平滑来表现SCNA的趋势（原文并没有采用）

tmp2 <- as.data.frame(focal.cna[,C2.sample]) # 取出类2的SCNA
tmp2$scna <- as.numeric(rowSums(tmp2)) # 对SCNA取行和，表示该亚型的在染色体上的SCNA水平
tmp2$loci <- 1:nrow(tmp2)
lw2 <- loess(scna ~ loci,data=tmp2,span = 0.02) # 直接使用行和会让折线变得很尖锐，所以这里我采用loess平滑来表现SCNA的趋势（原文并没有采用）


ylim <- pretty(c(min(c(range(lw1$fitted),range(lw2$fitted))),max(c(range(lw1$fitted),range(lw2$fitted))))) # 得到比较“温和”的y轴区间
xlab <- rownames(focal.cna) # 得到具体的染色体focal
xlab <- sapply(strsplit(xlab,".",fixed = T),"[",1) # 只取focal中“.”之前的字符串部分
xlab <- substr(xlab, 1, nchar(xlab) - 2); txt <- unique(xlab) # 去掉字符串末尾的2个字符，以获得带有q或者p的染色体臂字符串，此时取独特的字符作为x轴的标签
xlab <- factor(xlab, levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22p","22q")) # 赋予因子
txt <- txt[1:36]
pdf("./Output/Figure/focal-scna line chart.pdf", width = 17, height = 8)
par(bty="o", mgp = c(2,0.3,0), mar = c(3.1,3.1,2.1,3.1),tcl=-.25, font.main = 3, las = 1)
plot(1:730,     
     lw1$fitted,
     ylim = range(ylim),
     ylab = "Focal SCNA level",
     xlab = "",
     xaxt = "n",
     yaxs = "i", # 移除y轴和图像间的空隙
     xaxs = "i", # 移除x轴和图像间的空隙
     type = "l",
     cex = 1.5,
     col = "#DD492E")
lines(1:730,
      lw2$fitted,
      col = "#40548A", #蓝色
      cex = 1.5)
axis(side = 1, at = unique(cumsum(table(xlab))), labels = txt, cex.axis = 0.8, gap.axis = -1) # gap.axis为负数强制标记出所有的x轴标签
grid() # 添加网格
legend("topright", legend = c("c1","c2"), lty = 1, col = c("#DD492E","#40548A"), cex = 1, bty = "n")

invisible(dev.off())
