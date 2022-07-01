library(DESeq2)
# library(dplyr)
# library(stringr)
# library(ggplot2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(rstatix)
#library(VennDiagram)
#library(apeglm)
library(ggpubr) 
normalize_length = "TRUE"
select_length = "TRUE"
length_cutoff = 10000
# the cutoff of high expression level in ctrl
cut <- 20
site_cut <- 500
args = commandArgs(T)
treated = args[1]  # treated sample's mark
species = args[2]
if (length(args) > 3) {
  spikein="T"
  # size factor = 1/scale factor
  sizefactor = 1 / as.numeric(args[3:length(args)])
  print("sizefactor")
  print(sizefactor)
} else{
  spikein="F"
} 

######### functions ############
volcano <- function(df, outname, ylimit=FALSE, xlimit=c(FALSE)){
  pval_threshold <- 0.05
  logfc_threshold <- 0
  high_ex <- as.factor(abs(df$log2FoldChange) >=logfc_threshold & df$padj <=pval_threshold) 
  significance_up <-df %>% filter(log2FoldChange >=logfc_threshold & padj <= pval_threshold) %>% nrow()
  significance_down <-df %>% filter(log2FoldChange <= -logfc_threshold & padj <= pval_threshold) %>% nrow()
  volcano = ggplot(data=df, 
                   aes(x=log2FoldChange, y=-log10(padj), 
                       colour=high_ex)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_classic() +
    # geom_vline(xintercept = logfc_threshold) +
    # geom_vline(xintercept = -logfc_threshold) +
    # (yintercept = -log10(pval_threshold)) +
    scale_color_manual(values = c("gray","red"))+
    xlab("log2 fold change") + ylab("-log10 FDR") +
    labs(title = sub("(.*).pdf", "\\1", outname), 
         subtitle = paste0("n=", nrow(df), "   significance_up=", significance_up, "   significance_down=", significance_down)) +
    theme(legend.position = 'none') 
    #volcano = volcano + geom_text_repel(aes(label=ifelse(padj <= 0.05 & abs(log2FoldChange) >= logfc_threshold,name2, '')), max.overlaps = 15)
  if (min(df$log2FoldChange, na.rm = TRUE) >= 0) {
    volcano = volcano + xlim(-1, max(df$log2FoldChange, na.rm = TRUE)) + geom_point(color="red")
  }
  if (ylimit != FALSE) {
    volcano = volcano + coord_cartesian(ylim=c(0, ylimit))
  }
  if (xlimit[1] != FALSE) {
    volcano = volcano + xlim(xlimit)
  }
  ggsave(outname, width = 1800,height = 1800, dpi = 300, units = "px", plot = volcano)
}

boxplot_n <- function(counts_temp, n){
  counts_temp_wider = counts_temp %>% pivot_wider(names_from = samples, values_from = counts) %>% dplyr::select(-ID)
  error_site = sapply(counts_temp_wider,function(x) boxplot.stats(x)$stats[5])
  ymax = max(error_site)
  stat.test <- counts_temp  %>% wilcox_test(counts ~ samples) %>% adjust_pvalue() %>% add_significance("p.adj") %>% add_x_position(x = "samples")
  stat.test$y.position = seq(ymax*1.1, ymax*1.1 + (ymax*0.1) * (nrow(stat.test) - 1), ymax*0.1)
  p <- ggplot(counts_temp, aes(x = as.factor(samples), y = counts)) + ylab("Normalized Counts") + xlab("") +
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = samples))+
    geom_boxplot(width = 0.5,  outlier.colour = NA, aes(color = samples)) +
    theme_classic()+theme(legend.position = "none") + labs(subtitle = paste0("n=", nrow(counts_temp) / length(unique(counts_temp$samples))))
  p + coord_cartesian(ylim =c(0, ymax*1.3)) + stat_pvalue_manual(stat.test, label = "p.adj.signif",tip.length = 0)
  ggsave(paste0(n, "_boxplot.pdf"), dpi = 300, width = 1300, height = 2000, units = "px")
}


boxplot_2_counts<- function(counts, res, n){
  # plot with counts ----sense
  counts_temp = counts %>% dplyr::select(contains("_ave")) %>% rownames_to_column("ID") %>%  filter(ID %in% res$ID) %>%  pivot_longer(cols = -ID, values_to = "counts", names_to = "samples") 
  boxplot_n(counts_temp, n = paste0("gene_", n))
  # plot with counts ----antisense
  counts_temp = counts %>% dplyr::select(contains("_ave")) %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% 
    mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>%  filter(ID %in% res$ID) %>%  pivot_longer(cols = -ID, values_to = "counts", names_to = "samples")
  boxplot_n(counts_temp, n = paste0("PROMPTs_", n))
}



Ref <- read.table(paste0("~/reference/",species,"/RefSeq_all.txt"),sep = "\t")
Ref <- Ref %>% dplyr::select(transcriptID = V2, chr = V3, strand = V4, start = V5, end = V6, geneID = V13) %>% mutate(length = end - start) %>% 
  filter(!str_detect(string = chr, pattern = "_"))
Ref_N <- Ref  %>% filter(str_detect(transcriptID, "N"))
Ref_disctincted <- distinct(Ref_N, geneID, length, .keep_all = TRUE)
# read information of all transcripts
info <- read.csv(paste0("~/reference/",species,"/RefSeq_all.txt"), sep = "\t")
info_id <- info %>% dplyr::select(ID = name, name2, txStart, txEnd) %>% mutate(length = txEnd - txStart) 

{control = "CTRL"
  database <- read.csv('database.csv',header = T,check.names = FALSE)
  ##arrange the database, same as the range of featurecounts
  database_rep <- database %>% group_by(condition) %>% mutate(su = length(name)) %>% mutate(rep = 1:unique(su)) %>% 
    dplyr::select(-su) %>% mutate(condition_rep = paste(condition,"rep",rep,sep = "_")) %>% mutate(name = condition_rep)
  data <- read.table('featureCounts.txt',header = T, quote = '\t', skip = 1) 
  data[,7:ncol(data)] <- lapply(data[,7:ncol(data)], as.integer)
  Len <- data %>% dplyr::select(transcriptID = Geneid, Length)
  sampleNmaes <- as.factor(database_rep$condition_rep)
  countData <-as.matrix(data[,7:ncol(data)])
  colnames(countData) <- sampleNmaes
  rownames(countData) <- data$Geneid
  rownames(database) <- sampleNmaes
  database$condition <- as.factor(database$condition)
  database$condition <- relevel(database$condition,ref=control)
  ##DESeq2
  dds <- DESeqDataSetFromMatrix(countData, colData = database, design = ~condition)
  mcols(dds)$basepairs <- Len[,"Length"]
  if (spikein=="T") {
    # normalized with spike-in's size factor
    sizeFactors(dds) = sizefactor
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)
  }else{
    dds <- DESeq(dds)
    }
  vsd = vst(dds, blind=FALSE)
  vst.dds <- vst(dds)
  vst <- assay(vst.dds)
}
# distance between samples
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- database_rep$condition_rep
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p = pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
ggsave("distance.pdf", plot = p)


# for res analyze
uniq_condition = unique(database$condition)
ifelse("res" %in% uniq_condition, resc<-"T", resc<-"F")

vst <- vst %>%as.data.frame() %>% mutate(CTR_ave = rowMeans(dplyr::select(.,contains(control)))) %>% 
  mutate(`treated_ave` = rowMeans(dplyr::select(.,contains(treated))))
counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% mutate(CTR_ave = rowMeans(dplyr::select(.,contains(control)))) %>% 
  mutate(`treated_ave` = rowMeans(dplyr::select(.,contains(treated))))
if (resc == "T") {
  vst <- vst %>% mutate(`RES_ave` = rowMeans(dplyr::select(.,contains("res"))))
  counts <- counts %>% mutate(`RES_ave` = rowMeans(dplyr::select(.,contains("res"))))
}

# calculate the ctrl's genes normalized counts per kb
ctrl_count <- counts(dds, normalized = TRUE) %>% as.data.frame()
ctrl_count <- ctrl_count %>% dplyr::select(contains("CTRL")) %>% mutate(ctrl_ave = rowMeans(dplyr::select(ctrl_count,contains("CTRL")) )) %>% 
  dplyr::select(ctrl_ave) %>% rownames_to_column("transcriptID") %>% mutate(ID = transcriptID) 


dir.create('Differential')
setwd('Differential')

differentialData <- list()

med <-median(vst$CTR_ave)


##get TSS distance PN df 
{
  ##filter antisense TSS site and filter TSS site at 0-500 bp
  NM = read.table(paste0("~/reference/",species,"/RefSeq.bed"),sep = "\t") %>%  filter(!str_detect(V1,"_"))
  RT = read.table("../merge.bed",sep = "\t") %>% filter(!str_detect(V1,"_")) %>% mutate(V4=gsub(".*_(N.*)","\\1",V4))
  NM_p <- NM %>% 
    filter(V6 == "+") %>% 
    dplyr::select(V4,V2)
  NM_m = NM %>% 
    filter(V6 == "-") %>% 
    dplyr::select(V4,V3)
  RT_p <-RT %>% 
    filter(V6 == "-") %>% 
    dplyr::select(V4,V3)
  RT_m <-RT %>% 
    filter(V6 == "+") %>% 
    dplyr::select(V4,V2)
  Positive <- left_join(NM_p,RT_p,by = "V4") %>% mutate(site = V2 - V3) %>% dplyr::select(V4,site)
  Negavite <- left_join(NM_m,RT_m,by="V4") %>% mutate(site = V2-V3) %>%  dplyr::select(V4,site)
  PN <- rbind(Positive,Negavite) %>%  filter(site < site_cut & site > 0)  # %>% filter(V4 %in% name_list)
  PN_list <- PN$V4  ##the list of filtered transcript id
}
#plot the pca plot of all samples
plotPCA(vst.dds) + theme_classic()
ggsave("PCA.pdf", width = 1800, height = 1800, dpi = 300, units = "px" )

#plot volcano plot of antisense transcsripts
#function of plot volcanno plot
volcano_as <- function(df, outname, color=c("gray","red")){
  pval_threshold <- 0.05
  logfc_threshold <- 1
  high_ex <- as.factor(abs(df$log2FoldChange) >=logfc_threshold & df$padj <=pval_threshold) 
  significance_up <-df %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% nrow()
  significance_down <-df %>% filter(log2FoldChange <= -1 & padj <= 0.05) %>% nrow()
  volcano = ggplot(data=df, 
                   aes(x=log2FoldChange, y=-log10(padj), 
                       colour=high_ex)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_classic() + 
    # geom_vline(xintercept = logfc_threshold) +
    # geom_vline(xintercept = -logfc_threshold) +
    # geom_hline(yintercept = -log10(pval_threshold)) +
    xlim(c(-10, 10)) + 
    scale_color_manual(values = color)+
    xlab("log2 fold change") + ylab("-log10 FDR") +
    labs(title = sub("(.*).pdf", "\\1", outname), 
         subtitle = paste0("n=", nrow(df), "   significance_up=", significance_up, "   significance_down=", significance_down)) +
    theme(legend.position = 'none')
  #geom_text_repel(aes(label=ifelse(padj < 10**100 & abs(log2FoldChange) >= 4,SYMBOL, '')), max.overlaps = 30)
  ggsave(outname, width = 2000,height = 2000, dpi = 300, units = "px", plot = volcano)
}
# antisense valcanno plot iaa vs ctrl
as_df <- results(dds, contrast = c("condition",treated,"CTRL")) %>%as.data.frame() %>% 
  rownames_to_column("transcriptID") %>% filter(str_detect(transcriptID, "RT_"))
volcano_as(as_df, "antisense_treated_vs_CTRL.pdf")

if (resc == "T") {
  # antisense volcanno plot res vs iaa
  as_df <- results(dds, contrast = c("condition", "res",treated)) %>%as.data.frame() %>% 
    rownames_to_column("transcriptID") %>% filter(str_detect(transcriptID, "RT_"))
  volcano_as(as_df, "antisense_res_vs_treated.pdf")
  
  # antisense volcanno plot res vs ctrl
  as_df <- results(dds, contrast = c("condition","res","CTRL")) %>%as.data.frame() %>% 
    rownames_to_column("transcriptID") %>% filter(str_detect(transcriptID, "RT_"))
  volcano_as(as_df, "antisense_res_vs_CTRL.pdf")
}

##dplyr::select sense appear with AS
#get as list with AS significance upregulate
namelist = results(dds,contrast = c("condition",treated,"CTRL")) %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(startsWith(ID,"RT_")) %>%
  mutate(ID=gsub("RT_(N.*)","\\1", ID)) %>% filter(log2FoldChange >=1 & padj <=0.05) %>% distinct(ID, .keep_all = TRUE)
namelist <- namelist$ID

# select the distance of antisense TSS 0-300bp
{
  ##focus on 0-300bp
  PN2 <- rbind(Positive,Negavite)  %>% distinct(V4, .keep_all = TRUE) %>% filter(site >=-2000 &site<= 2000) %>% filter(V4 %in% namelist)
  PN2_300<-PN2 %>% filter(site>0 & site < site_cut) %>% distinct(V4, .keep_all = TRUE)
  list_300 <- PN2_300$V4
  write.table(PN2_300$V4,"list300.txt",quote = FALSE, row.names = FALSE,col.names = FALSE)
  x <- ggplot(PN2,aes(site))+
    geom_histogram(binwidth = 10)+
    theme(legend.position = 'none') +
    xlab("distance(bp)")+ 
    xlim(c(-2000, 2000))+
    labs(title = "The distribution of antisense TSS sits" ,subtitle = paste0("n=", nrow(PN2), "   0-", site_cut, "(bp)=", nrow(PN2_300)))+
    theme_classic() + theme()
  ggsave('RT_TSS_site.pdf',width = 2400,height = 1800, dpi = 300, units = "px", plot = x)
}



#get the result that only have significance antisense appear

res_300 <- results(dds,contrast = c("condition",treated,"CTRL")) %>% as.data.frame() %>% rownames_to_column("ID") %>%
  filter(ID %in% list_300)  %>% left_join(info_id, by = "ID") %>% left_join(ctrl_count, by = "ID") %>%
  mutate(ctrl_ave = ctrl_ave*1000/length) %>% distinct(ID, .keep_all = TRUE)
res_anti <-  counts %>% dplyr::select(contains("_ave")) %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% 
  mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>% filter(ID %in% list_300) %>% pivot_longer(cols = -ID, values_to = "log", names_to = "samples") %>% 
  mutate(log = log2(log + 1)) %>% distinct(ID, .keep_all = TRUE)
ggplot(res_anti, aes(x = as.factor(samples), y = log, color = samples)) + ylab("log2(normalized counts + 1)") + xlab("") + stat_boxplot(geom = "errorbar",width=0.2) +
  geom_boxplot(width = 0.5,  outlier.colour = NA)+ theme_classic()+theme(legend.position = "none")
ggsave("antisense_boxplot.pdf", dpi = 300)

# analyze data with quantile of ctrl expression level 
# output with quantile volcanol and quantile gene list
dir.create("quantile")
setwd("quantile")
quan <- quantile(res_300$ctrl_ave)
n = 0
quan_res = res_300 %>% mutate(quantile = ntile(ctrl_ave, 4) / 4) 

for (i in {1:4}) {
  res = quan_res %>% filter(quantile == i/4)
  write.table(res$ID ,paste0(n, "_quantile.txt"),quote = FALSE, col.names = FALSE, row.names = FALSE)
  volcano(res, paste0(n,"~",n+0.25,"_sense.pdf"), 50,c(-4, 4))
  res_anti = results(dds,contrast = c("condition",treated,"CTRL")) %>% as.data.frame() %>% 
    rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>% filter(ID %in% res$ID)%>% 
    left_join(info_id, by = "ID") %>% distinct(ID, .keep_all = TRUE)
  volcano(res_anti, paste0("PROMPTs_", n,"~",n+0.25,"_volcanno.pdf"), 150, c(-5,15))
  write.csv(res, paste0(n, "_quantile.csv"))
  boxplot_2_counts(counts, res, n)
  n = n +0.25
}


{
  #  plot for 0~0.75
  n="0~0.75"
  res = quan_res %>% filter(quantile %in% c(0.25, 0.5, 0.75))
  volcano(res, paste0(0,"~", 0.75,"_sense.pdf"))
  boxplot_2_counts(counts, res, n)
  res_anti = results(dds,contrast = c("condition",treated,"CTRL")) %>% as.data.frame() %>% 
    rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>% filter(ID %in% res$ID)%>% 
    left_join(info_id, by = "ID") %>% distinct(ID, .keep_all = TRUE)
  volcano(res_anti, paste0("PROMPTs_", 0,"~",0.75,"_volcanno.pdf"))
  # plot for all PROMPTs
  n="all"
  res_anti = results(dds,contrast = c("condition",treated,"CTRL")) %>% as.data.frame() %>% 
    rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>%  distinct(ID, .keep_all = TRUE)
  boxplot_2_counts(counts, res = res_anti, n)
}



setwd("..")

#for res conditon
if (resc == "T") {
  res_resc <- results(dds,contrast = c("condition","res","CTRL")) %>% as.data.frame() %>% rownames_to_column("ID") %>%
    filter(ID %in% list_300)  %>% left_join(info_id, by = "ID") %>% left_join(ctrl_count, by = "ID") %>% mutate(ctrl_ave = ctrl_ave*1000/length) %>% distinct(ID, .keep_all = TRUE)
  res_resc_as <- results(dds,contrast = c("condition","res","CTRL")) %>% as.data.frame() %>% rownames_to_column("ID") %>%
    filter(str_detect(ID, "RT_")) %>% mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>% left_join(info_id, by = "ID") %>% distinct(ID, .keep_all = TRUE)
  volcano(res_resc %>% filter(ctrl_ave >= cut), paste0("resc_",cut ,"(normalized counts per kb)_sense.pdf"))
}


# n quantile boxplot
quan_box <- function(n, df, outname){
  df1 <- df %>% mutate(quantile = ntile(ctrl_ave, n) / n) %>% 
    group_by(quantile) %>% mutate(ctrl = quantile) %>% arrange(quantile)
  box <- ggplot(df1, aes(x=as.factor(quantile), y=log2FoldChange, color = quantile))+
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = quantile)) + geom_boxplot(outlier.colour = NA) + 
    theme_classic() +  xlab("quantile of the ctrl's expression level(normalized counts per kb)")+
    theme(legend.position = "none")+
    ylab("log2(treated/ctrl)")  # + (yintercept = 0, color = "grey20") 
  ggsave(outname, dpi = 300)
  # plot with expression level 
  box <- ggplot(df1, aes(x=as.factor(quantile), y=ctrl_ave, color = quantile))+
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = quantile)) +  geom_boxplot(outlier.colour = NA) + 
    theme_classic() +  xlab("quantile of the ctrl's expression level(normalized counts per kb)")+
    theme(legend.position = "none") + ylim(0,80) 
  ggsave(paste0("gene_", outname), dpi = 300)
  #save quantile 4-10 for heatmap plot
  df2<- df1 %>% filter(quantile >=0.4)
  write.table(df2$ID,file = "quan4_10.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

quan_box(4, res_300, "as_box4_plot.pdf")
quan_box(10, res_300, "as_box10_plot.pdf")





write.table(res_300$ID, file = "list_300.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# write transcriptID separated with strands 
Ref_300 <- Ref %>% filter(transcriptID %in% res_300$ID)
Ref_300_fw <- filter(Ref_300, strand == "+")
Ref_300_rev <- filter(Ref_300, strand == "-")
write.table(Ref_300_fw$transcriptID, file = "list300_fw.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Ref_300_rev$transcriptID, file = "list300_rev.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


# plot MA for sense transcripts appearence with antisense
MA_res <- res_300 %>% dplyr::select(ctrl_ave, log2FoldChange) %>% mutate(padj = res_300$padj <= 0.05) 

plotMA(MA_res, colLine = "grey40", colNonSig = "gray60", colSig = "red", xlab = "CTRL expression level(normalized counts per kb)", ylab = "log2(treated/CTRL)",
       main = paste0("MA plot of transcripts with antisense appeared\nn=",nrow(MA_res), 
                     "   significance up=" ,nrow(filter(MA_res, padj == "TRUE" & log2FoldChange >0)),
                     "   significance down=" ,nrow(filter(MA_res, padj == "TRUE" & log2FoldChange <0))))
ggsave("S_withAS_MAplot.pdf", width = 2600,height = 1600, dpi = 300, units = "px")


## quantile of gene expression level  for PROMPTS

gene_quan = res_300 %>% mutate(quantile = ntile(ctrl_ave, 4)) %>% dplyr::select(ID, quantile)
anti_counts <-  counts %>% dplyr::select(contains("_ave")) %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_")) %>% 
  mutate(ID = gsub("RT_(.*)", "\\1", ID)) %>% filter(ID %in% list_300) %>% merge(gene_quan, by = "ID")

anti_counts_longer = anti_counts %>% 
  pivot_longer(values_to = "Normalized Counts", cols = c(CTR_ave, treated_ave), names_to = "Samples")
if (resc == "T") {
  anti_counts_longer = anti_counts_longer %>% dplyr::select(-RES_ave) 
} 


box <- ggplot(anti_counts, aes(x=as.factor(quantile), y=treated_ave, color = quantile))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = quantile)) +  geom_boxplot(outlier.colour = NA) + 
  theme_classic() +  xlab("quantile of the ctrl's expression level(normalized counts per kb)")+
  theme(legend.position = "none") + ylim(0,500)
box
box <- ggplot(anti_counts, aes(x=as.factor(quantile), y=CTR_ave, color = quantile))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = quantile)) +  geom_boxplot(outlier.colour = NA) + 
  theme_classic() +  xlab("quantile of gene expression level")+
  theme(legend.position = "none") + ylim(0,500)
box

box <- ggplot(anti_counts_longer, aes(x=as.factor(quantile), y=`Normalized Counts`, color = Samples))+
  geom_boxplot(outlier.colour = NA) + 
  theme_classic() +  xlab("quantile of gene expression level")+ 
  scale_color_manual(values = c("blue","red"))+
  ylim(0,450)
box
ggsave("quan_anti.pdf", dpi = 600, width = 5.5, height = 4, plot = box)
