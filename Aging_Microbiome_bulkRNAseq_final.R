library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(scales)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(genefilter)
library(geneplotter)
library(plyr)
library(VennDiagram)
library(devtools)
library(forcats)
library(tidyr)
library(variancePartition)
library(tibble)
library(reshape2)
library(dplyr)
set.seed(42)

## Summary
# Samples
## Mice were exposed to sequential antibiotic treatment followed by colonization with either microbiome from same-aged mice (isochronic microbiome-iMB) or with the microbiome from 8 week old mice (young microbiome-yMB). 
## Colonic and Ileum tissue was collected from each mice at week 40 during treatment and at week 120 (final) which is define as Timepoint.

# load Ileum dataset
# load Merge counts and metadata
setwd("/Users/joana/Documents/IKMB/RNAseq/Ageing_final/AgingMiceBulkRNAIlleum/")
mergedCounts <- read.delim("Ileum_mergedCounts.txt")

rownames(mergedCounts)<-mergedCounts$GeneIDs
mergedCounts$GeneIDs=NULL

infotable<-read.delim("Ileum_metaData.txt")


SI<-subset(infotable, infotable$Organ %in% 'Si8')
SI$SampleID <- paste0("X", SI$SampleID)
count_data_SI<-mergedCounts[, as.character(SI$SampleID)]

write.csv2(count_data_SI, file='SI_Ileum_mergedCounts.csv')
write.csv2(SI, file='SI_Ileum_metaData.csv')

SI_treat<-subset(SI, !SI$Treatment %in% 'control')
count_data_SI_treat<-mergedCounts[, as.character(SI_treat$SampleID)]


## perform PCA
# use DEseq2 to transform data
dds_counts=DESeqDataSetFromMatrix(countData = count_data_SI_treat, colData =SI_treat, design = ~ 1 )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized=TRUE)
rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

write.csv(rlog.norm.counts, file= 'SI_Ileum_Norm_counts_treat.csv')
#rlog.norm.counts<-read.csv('SI_Ileum_Norm_counts.csv')

pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
pc_df$SampleID <- SI_treat$SampleID
pc_df$Timepoint <- SI_treat$Timepoint
pc_df$MouseID <- SI_treat$MouseID
pc_df$Treatment <-SI_treat$Treatment

pc_df <- count_data_SI_treat %>%
  colSums() %>%
  data.frame("SampleID" = names(.), "Reads" = .) %>%
  merge(pc_df, by = "SampleID")

pc_df <- count_data_SI_treat %>%
  '>' (0) %>%
  colSums() %>%
  data.frame("SampleID" = names(.), "Genes" = .) %>%
  merge(pc_df, by = "SampleID")


write.csv(pc_df, file= 'PC_Ileum_treat.csv')

eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)


#52 MouseID potential outlier


pdf("pca_plot_ageing_Ileum_treat.pdf", width = 5, height = 3.5)
P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, label=MouseID, color=Treatment, shape=Timepoint)) +
  geom_point(size=4) + geom_text(hjust=-0.5, nudge_y = 0.4) + xlab("PC1 (30.83% variance)") + 
  ylab("PC3 (7.82% variance)")
P <- P + scale_color_manual(values = c( "#D55E00","#0072B2",'grey'))
P <- P + scale_size_area(max_size=4)
P <- P + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + theme_bw()
P
dev.off()


# Perform variancepartition
#prepare tables
isexpr <- rowSums(fpm(dds_counts) > 1) >= 0.5 * ncol(dds_counts)
quantvst <- vst(dds_counts[isexpr, ])
quantvst <- assay(quantvst)

#compute variance partition
SI_treat$outlier <- ifelse(SI_treat$MouseID == 52, "y", "n")
pc_df$outlier <-SI_treat$outlier
SI_treat$Reads <- pc_df$Reads
SI_treat$Genes <- pc_df$Genes

form <- ~ Timepoint + Treatment + outlier + Reads + Genes

info_SI_treat<- SI_treat[, c(12, 13, 15,16,17 )]

varPart <- fitExtractVarPartModel(quantvst, form,info_SI_treat )

#export results
varMean <- colMeans(varPart)
varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
varPart_results <- melt(varPart)
varPart_results <- subset(varPart_results, varPart_results$variable != "Residuals")

write.table(varPart_sorted, "variance_partition_Ileum_simple_DESeq_batch.txt", sep = '\t', quote = FALSE)

varMean <- as.data.frame(varMean)
varMean$Parameter <- rownames(varMean)
varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]

write.table(varMean, "variance_partition_mean_Ileum_DESeq_batch.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Violin plot
varPart_sorted <- add_column(varPart_sorted, "ID" = "DESeq")

varPart_sorted$ID<-rownames(varPart_sorted)
no_factors <- melt(varPart_sorted, id.vars = "ID")
no_factors$variable <- factor(no_factors$variable)
no_factors$value <- as.numeric(no_factors$value)


pdf("Variance_partition_violin_Ileum_DESeq.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()

no_factors2<-subset(no_factors, !no_factors$variable %in% 'outlier')

pdf("Variance_partition_violin_Ileum_DESeq_nooutlier.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors2, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()

## correlation between principal components and variables
form <- ~ Timepoint + Treatment + outlier + Reads + Genes+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8

Corr_mat <- canCorPairs(form, pc_df) %>%
  .[paste0(rep("PC", 8), 1:8), c( "Genes", "Reads", 'Timepoint','Treatment', 'outlier' )] %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2))

Corr_mat$Components <- factor(Corr_mat$Components, levels = paste0(rep("PC", 8), 1:8) %>% rev())
Corr_mat$variable <- factor(Corr_mat$variable, levels = c(,"Genes", "Reads", 'Timepoint','Treatment', 'outlier'))

var_names <- c("Timepoint" = "Timepoint", 
               "Treatment" = "Treatment",
               "outlier" = "outlier",
               "Reads" = "Reads",
               'Genes'='Genes')

pdf("PC_correlation_pca_metadata.pdf", width = 7, height = 4)
ggplot(Corr_mat, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black", size = 3) +
  scale_x_discrete(labels = as_labeller(var_names)) +
  scale_size(range = c(1, 11), guide = "none") +
  scale_fill_gradient(low = "white", high = "#2171B5") +
  scale_color_gradient(low = "white", high = "#2171B5") +
  labs(fill = "Correlation", color = "Correlation", size = "Correlation", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        panel.border = element_rect(color = "grey"),
        axis.text = element_text(size = 11)) +
  coord_equal()
dev.off()


## remove outlier and repeat
SI_treat_no<-subset(SI_treat, !SI_treat$outlier %in% 'y')
count_data_SI_treat_no<-mergedCounts[, as.character(SI_treat_no$SampleID)]


## perform PCA
# use DEseq2 to transform data
dds_counts=DESeqDataSetFromMatrix(countData = count_data_SI_treat_no, colData =SI_treat_no, design = ~ 1 )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized=TRUE)
rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

write.csv(rlog.norm.counts, file= 'SI_Ileum_Norm_counts_treat_noOutlier.csv')
#rlog.norm.counts<-read.csv('SI_Ileum_Norm_counts.csv')

pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
pc_df$SampleID <- SI_treat_no$SampleID
pc_df$Timepoint <- SI_treat_no$Timepoint
pc_df$MouseID <- SI_treat_no$MouseID
pc_df$Treatment <-SI_treat_no$Treatment

pc_df <- count_data_SI_treat_no %>%
  colSums() %>%
  data.frame("SampleID" = names(.), "Reads" = .) %>%
  merge(pc_df, by = "SampleID")

pc_df <- count_data_SI_treat_no %>%
  '>' (0) %>%
  colSums() %>%
  data.frame("SampleID" = names(.), "Genes" = .) %>%
  merge(pc_df, by = "SampleID")


write.csv(pc_df, file= 'PC_Ileum_treat_noOutlier.csv')

eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)
eigs[4] / sum(eigs)


pdf("pca_plot_ageing_Ileum_treat_noOutlier_final.pdf", width = 5, height = 3.5)
P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC4, label=MouseID, fill=Treatment, shape=Timepoint)) +
  geom_point(size=4) + geom_text(hjust=-0.5, nudge_y = 0.4) + xlab("PC1 (30.83% variance)") + 
  ylab("PC4 (7.16% variance)")
P <- P + scale_fill_manual(values = c("#D55E00", "#0072B2")) 
P <- P + scale_shape_manual(values = c(21, 24, 22, 25))
P <- P + scale_size_area(max_size=4)
P <- P + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + theme_bw()
P
dev.off()


# Perform variancepartition
#prepare tables
isexpr <- rowSums(fpm(dds_counts) > 1) >= 0.5 * ncol(dds_counts)
quantvst <- vst(dds_counts[isexpr, ])
quantvst <- assay(quantvst)

#compute variance partition
SI_treat_no$Reads <- pc_df$Reads
SI_treat_no$Genes <- pc_df$Genes

form <- ~ Timepoint + Treatment  + Reads + Genes

info_SI_treat<- SI_treat_no[, c(12, 13, 16,17 )]

varPart <- fitExtractVarPartModel(quantvst, form,info_SI_treat )

#export results
varMean <- colMeans(varPart)
varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
varPart_results <- melt(varPart)
varPart_results <- subset(varPart_results, varPart_results$variable != "Residuals")

write.table(varPart_sorted, "variance_partition_Ileum_simple_DESeq_batch_noOutlier.txt", sep = '\t', quote = FALSE)

varMean <- as.data.frame(varMean)
varMean$Parameter <- rownames(varMean)
varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]

write.table(varMean, "variance_partition_mean_Ileum_DESeq_batch_noOutlier.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Violin plot
varPart_sorted$ID<-rownames(varPart_sorted)
no_factors <- melt(varPart_sorted, id.vars = "ID")
no_factors$variable <- factor(no_factors$variable)
no_factors$value <- as.numeric(no_factors$value)


pdf("Variance_partition_violin_Ileum_DESeq_noOutlier.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()


## correlation between principal components and variables
form <- ~ Timepoint + Treatment +  Reads + Genes+ PC1 + PC2 + PC3 + PC4 + PC5 

Corr_mat <- canCorPairs(form, pc_df) %>%
  .[paste0(rep("PC", 5), 1:5), c( "Genes", "Reads", 'Timepoint','Treatment' )] %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2))

Corr_mat$Components <- factor(Corr_mat$Components, levels = paste0(rep("PC", 5), 1:5) %>% rev())
Corr_mat$variable <- factor(Corr_mat$variable, levels = c("Genes", "Reads", 'Timepoint','Treatment'))

var_names <- c("Timepoint" = "Timepoint", 
               "Treatment" = "Treatment",
               "Reads" = "Reads",
               'Genes'='Genes')

pdf("PC_correlation_pca_metadata_noOutlier.pdf", width = 7, height = 4)
ggplot(Corr_mat, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black", size = 3) +
  scale_x_discrete(labels = as_labeller(var_names)) +
  scale_size(range = c(1, 11), guide = "none") +
  scale_fill_gradient(low = "white", high = "#2171B5") +
  scale_color_gradient(low = "white", high = "#2171B5") +
  labs(fill = "Correlation", color = "Correlation", size = "Correlation", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        panel.border = element_rect(color = "grey"),
        axis.text = element_text(size = 11)) +
  coord_equal()
dev.off()

## selection of genes contributing to Treatment differences 
sig_Treatment<-subset(varPart_sorted, varPart_sorted$Treatment > 0.10)

#Remove genes effected by Timepoint, Genes or Reads
sig_Timepoint<-subset(varPart_sorted, varPart_sorted$Timepoint > 0.10)
sig_Reads<-subset(varPart_sorted, varPart_sorted$Reads > 0.10)
sig_Genes<-subset(varPart_sorted, varPart_sorted$Genes > 0.10)

remove<-rbind(sig_Timepoint,sig_Reads, sig_Genes )
remove_genes <- remove %>%
  distinct(ID) %>%
  pull(ID) 

Treatment_Genes_only<- subset(sig_Treatment, !rownames(sig_Treatment) %in% remove_genes)

## preform enrichment analyis
# GO enrichment
markers_entrez <- bitr(Treatment_Genes_only$ID,
                       fromType = "ENSEMBL",
                       toType="ENTREZID",
                       OrgDb=org.Mm.eg.db)
Treatment_Genes_only$ENSEMBL<-Treatment_Genes_only$ID
Treatment_Genes_only<-merge(Treatment_Genes_only, markers_entrez,by.x='ENSEMBL')

markers_entrez <- bitr(Treatment_Genes_only$ID,
                       fromType = "ENSEMBL",
                       toType="SYMBOL",
                       OrgDb=org.Mm.eg.db)
Treatment_Genes_only$ENSEMBL<-Treatment_Genes_only$ID
Treatment_Genes_only<-merge(Treatment_Genes_only, markers_entrez,by.x='ENSEMBL')

# is the gene up or downregulated in yMB?
dds_counts=DESeqDataSetFromMatrix(countData = count_data_SI_treat_no, colData =SI_treat_no, design = ~  Treatment )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
dds_norm=DESeq(dds_counts)
dds_normr=results(dds_norm)

write.csv(dds_normr, file="DE_Ileum_8wVSageing_noOutlier.csv")

dds_normr$ENSEMBL<-row.names(dds_normr)
Down_comp<-as.data.frame(subset(dds_normr, log2FoldChange < 0))
Up_comp<-as.data.frame(subset(dds_normr, log2FoldChange > 0))

Down_comp$signal<-c('Upregulated in yMB')
Up_comp$signal<-c('Downregulated in yMB')
total<-rbind(Down_comp,Up_comp)

total_signal<- total[, c(7, 8 )]

Treatment_Genes_only_signal<-merge(Treatment_Genes_only, total_signal, by='ENSEMBL')

# Compare clusterProfiler
CompGO <- as.data.frame(compareCluster(ENTREZID~signal, data=Treatment_Genes_only_signal,OrgDb='org.Mm.eg.db', fun="enrichGO", ont = "BP"))

ENTREZID<-as.data.frame(CompGO$geneID)
nmax <- max(stringr::str_count(ENTREZID$`CompGO$geneID`, "\\/")) + 1
ENTREZID2<-tidyr::separate(ENTREZID, `CompGO$geneID`, paste0("col", seq_len(nmax)), sep = "\\/", fill = "right")
rows<- nrow(ENTREZID2)
GENEID <- as.data.frame(matrix(mapIds(org.Mm.eg.db, as.character(unlist(ENTREZID2)), "SYMBOL","ENTREZID"), rows))
col<-colnames(GENEID)
GENEID2<-unite(GENEID, col='GENEID', col, sep=',')
GENEID2$GENEID<-gsub("NULL,","",as.character(GENEID2$GENEID))
GENEID2$GENEID<-gsub(",NULL","",as.character(GENEID2$GENEID))
CompGO$GENEID<-GENEID2$GENEID
write.csv2(CompGO, file='GOComparison_Ileum_Treatment_noOutlier_WithSignal.csv')

CompGO<-read.csv2('GOComparison_Colon_TreatmentNOTime_WithSignal_interest.csv')

CompGO %>% group_by(Cluster) %>% top_n(n= 10, wt = -p.adjust) -> terms
CompGOresults_sub<-subset(CompGO, CompGO$p.adjust <0.05)
CompGOresults_sub2<-subset(CompGOresults_sub, CompGOresults_sub$Count>=3)
#tmp <- CompGOresults_sub2[CompGOresults_sub2$Description %in% terms$Description,]
CompGOresults_sub2$Description <- ifelse(nchar(CompGOresults_sub2$Description)>60,
                                         paste(substr(CompGOresults_sub2$Description, 1, 60),"[...]",sep=""),
                                         CompGOresults_sub2$Description)
CompGOresults_sub2$Description <- factor(CompGOresults_sub2$Description,levels=unique(CompGOresults_sub2$Description))
CompGOresults_sub4<-CompGOresults_sub2 %>% group_by(Cluster) %>% top_n(n= 20, wt = -p.adjust)

CompGOresults_sub4$Cluster2 = with(CompGOresults_sub4, reorder(Cluster))


pdf("GOComparison_TreatmentnoOutlier_WithSignal_all.pdf", width = 8, height =6)
ggplot(CompGOresults_sub4, aes(x = p.adjust, y = Description,size = Count, fill=signal)) +
  geom_point(shape=21) +
  scale_fill_manual(values=c("cornflowerblue", "brown1"))+
  ylab(NULL) +
  theme_bw() +
  #facet_wrap(~ Cluster, ncol = 2, scales = "free") +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## colon
count_data<-read.csv2('/Users/joana/Documents/IKMB/RNAseq/Ageing_final/count_table_ageing.csv')
col_data<-read.csv2('/Users/joana/Documents/IKMB/RNAseq/Ageing_final/infotable_ageing.csv')

table(col_data$Tissue)

rownames(count_data)<-count_data$Geneid
count_data$Geneid=NULL
count_data$gene_name=NULL

Colon<-subset(col_data, col_data$Tissue %in% 'colon')
count_data_Colon<-count_data[, as.character(Colon$Sample)]

Colon2 <- Colon %>%
  separate(Name, into = c("Timepoint", "Treatment", "MouseID", 'Tissue'), sep = "_")

Colon<- Colon2[, c(1, 2,3,4,5 )]

Colon$Treatment<-Colon2$mb
write.csv2(count_data_Colon, file='Colon_mergedCounts.csv')
write.csv2(Colon, file='Colon_metaData.csv')


Colon_treat<-subset(Colon, !Colon$Treatment %in% 'control')
count_data_Colon_treat<-count_data_Colon[, as.character(Colon_treat$Sample)]


## perform PCA
# use DEseq2 to transform data
dds_counts=DESeqDataSetFromMatrix(countData = count_data_Colon_treat, colData =Colon_treat, design = ~ 1 )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized=TRUE)
rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

write.csv(rlog.norm.counts, file= 'Colon_Norm_counts_treat.csv')
#rlog.norm.counts<-read.csv('SI_Ileum_Norm_counts.csv')

pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
pc_df$Sample <- Colon_treat$Sample
pc_df$Timepoint <- Colon_treat$Timepoint
pc_df$MouseID <- Colon_treat$MouseID
pc_df$Treatment <-Colon_treat$Treatment

pc_df <- count_data_Colon_treat %>%
  colSums() %>%
  data.frame("Sample" = names(.), "Reads" = .) %>%
  merge(pc_df, by = "Sample")

pc_df <- count_data_Colon_treat %>%
  '>' (0) %>%
  colSums() %>%
  data.frame("Sample" = names(.), "Genes" = .) %>%
  merge(pc_df, by = "Sample")


write.csv(pc_df, file= 'PC_Colon_treat.csv')

eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)


#52 MouseID potential outlier


pdf("pca_plot_ageing_Colon_treat.pdf", width = 5, height = 3.5)
P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, label=MouseID, color=Treatment, shape=Timepoint)) +
  geom_point(size=4) + geom_text(hjust=-0.5, nudge_y = 0.4) + xlab("PC1 (29.97% variance)") + 
  ylab("PC2 (11.36% variance)")
P <- P + scale_color_manual(values = c( "#D55E00","#0072B2",'grey'))
P <- P + scale_size_area(max_size=4)
P <- P + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + theme_bw()
P
dev.off()


# Perform variancepartition
#prepare tables
isexpr <- rowSums(fpm(dds_counts) > 1) >= 0.5 * ncol(dds_counts)
quantvst <- vst(dds_counts[isexpr, ])
quantvst <- assay(quantvst)

#compute variance partition
Colon_treat$outlier <- ifelse(Colon_treat$MouseID == 52, "y", "n")
pc_df$outlier <-Colon_treat$outlier
Colon_treat$Reads <- pc_df$Reads
Colon_treat$Genes <- pc_df$Genes

form <- ~ Timepoint + Treatment + outlier + Reads + Genes

info_Colon_treat<- Colon_treat[, c(2, 3, 6,7 , 8)]

varPart <- fitExtractVarPartModel(quantvst, form,info_Colon_treat )

#export results
varMean <- colMeans(varPart)
varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
varPart_results <- melt(varPart)
varPart_results <- subset(varPart_results, varPart_results$variable != "Residuals")

write.table(varPart_sorted, "variance_partition_Colon_simple_DESeq_batch.txt", sep = '\t', quote = FALSE)

varMean <- as.data.frame(varMean)
varMean$Parameter <- rownames(varMean)
varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]

write.table(varMean, "variance_partition_mean_Colon_DESeq_batch.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Violin plot
varPart_sorted$ID<-rownames(varPart_sorted)
no_factors <- melt(varPart_sorted, id.vars = "ID")
no_factors$variable <- factor(no_factors$variable)
no_factors$value <- as.numeric(no_factors$value)


pdf("Variance_partition_violin_Colon_DESeq.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()

no_factors2<-subset(no_factors, !no_factors$variable %in% 'outlier')

pdf("Variance_partition_violin_Colon_DESeq_nooutlier.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors2, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()

## correlation between principal components and variables
form <- ~ Timepoint + Treatment + outlier + Reads + Genes+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8

Corr_mat <- canCorPairs(form, pc_df) %>%
  .[paste0(rep("PC", 8), 1:8), c( "Genes", "Reads", 'Timepoint','Treatment', 'outlier' )] %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2))

Corr_mat$Components <- factor(Corr_mat$Components, levels = paste0(rep("PC", 8), 1:8) %>% rev())
Corr_mat$variable <- factor(Corr_mat$variable, levels = c("Genes", "Reads", 'Timepoint','Treatment', 'outlier'))

var_names <- c("Timepoint" = "Timepoint", 
               "Treatment" = "Treatment",
               "outlier" = "outlier",
               "Reads" = "Reads",
               'Genes'='Genes')

pdf("PC_correlation_pca_metadata_Colon.pdf", width = 7, height = 4)
ggplot(Corr_mat, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black", size = 3) +
  scale_x_discrete(labels = as_labeller(var_names)) +
  scale_size(range = c(1, 11), guide = "none") +
  scale_fill_gradient(low = "white", high = "#2171B5") +
  scale_color_gradient(low = "white", high = "#2171B5") +
  labs(fill = "Correlation", color = "Correlation", size = "Correlation", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        panel.border = element_rect(color = "grey"),
        axis.text = element_text(size = 11)) +
  coord_equal()
dev.off()


## remove outlier and repeat
Colon_treat_no<-subset(Colon_treat, !Colon_treat$outlier %in% 'y')
count_data_Colon_treat_no<-count_data[, as.character(Colon_treat_no$Sample)]


## perform PCA
# use DEseq2 to transform data
dds_counts=DESeqDataSetFromMatrix(countData = count_data_Colon_treat_no, colData =Colon_treat_no, design = ~ 1 )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized=TRUE)
rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

write.csv(rlog.norm.counts, file= 'Colon_Norm_counts_treat_noOutlier.csv')
#rlog.norm.counts<-read.csv('SI_Ileum_Norm_counts.csv')

pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
pc_df$Sample <- Colon_treat_no$Sample
pc_df$Timepoint <- Colon_treat_no$Timepoint
pc_df$MouseID <- Colon_treat_no$MouseID
pc_df$Treatment <-Colon_treat_no$Treatment

pc_df <- count_data_Colon_treat_no %>%
  colSums() %>%
  data.frame("Sample" = names(.), "Reads" = .) %>%
  merge(pc_df, by = "Sample")

pc_df <- count_data_Colon_treat_no %>%
  '>' (0) %>%
  colSums() %>%
  data.frame("Sample" = names(.), "Genes" = .) %>%
  merge(pc_df, by = "Sample")


write.csv(pc_df, file= 'PC_Colon_treat_noOutlier.csv')

eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)
eigs[4] / sum(eigs)


pdf("pca_plot_ageing_Colon_treat_noOutlier_final.pdf", width = 5, height = 3.5)
P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC3, label=MouseID, fill=Treatment, shape=Timepoint)) +
  geom_point(size=4) + geom_text(hjust=-0.5, nudge_y = 0.4) + xlab("PC1 (33.75% variance)") + 
  ylab("PC3 (8.80% variance)")
P <- P + scale_fill_manual(values = c("#D55E00", "#0072B2")) 
P <- P + scale_shape_manual(values = c(21, 24, 22, 25))
P <- P + scale_size_area(max_size=4)
P <- P + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + theme_bw()
P
dev.off()


# Perform variancepartition
#prepare tables
isexpr <- rowSums(fpm(dds_counts) > 1) >= 0.5 * ncol(dds_counts)
quantvst <- vst(dds_counts[isexpr, ])
quantvst <- assay(quantvst)

#compute variance partition
Colon_treat_no$Reads <- pc_df$Reads
Colon_treat_no$Genes <- pc_df$Genes

form <- ~ Timepoint + Treatment  + Reads + Genes

info_Colon_treat<- Colon_treat_no[, c(2, 3, 7,8 )]

varPart <- fitExtractVarPartModel(quantvst, form,info_Colon_treat )

#export results
varMean <- colMeans(varPart)
varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
varPart_results <- melt(varPart)
varPart_results <- subset(varPart_results, varPart_results$variable != "Residuals")

write.table(varPart_sorted, "variance_partition_Colon_simple_DESeq_batch_noOutlier.txt", sep = '\t', quote = FALSE)

varMean <- as.data.frame(varMean)
varMean$Parameter <- rownames(varMean)
varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]

write.table(varMean, "variance_partition_mean_Colon_DESeq_batch_noOutlier.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Violin plot
varPart_sorted$ID<-rownames(varPart_sorted)
no_factors <- melt(varPart_sorted, id.vars = "ID")
no_factors$variable <- factor(no_factors$variable)
no_factors$value <- as.numeric(no_factors$value)


pdf("Variance_partition_violin_Colon_DESeq_noOutlier.pdf", width = 3, height = 2.5)
p <- ggplot(no_factors, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.07) +
  labs(y = "Variance explained") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"), 
        axis.title = element_text(size = 10)) +
  theme(strip.placement = "outside", 
        strip.background  = element_blank(), 
        panel.border = element_blank())
plot(p)
dev.off()


## correlation between principal components and variables
form <- ~ Timepoint + Treatment +  Reads + Genes+ PC1 + PC2 + PC3 + PC4 + PC5 

Corr_mat <- canCorPairs(form, pc_df) %>%
  .[paste0(rep("PC", 5), 1:5), c( "Genes", "Reads", 'Timepoint','Treatment' )] %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2))

Corr_mat$Components <- factor(Corr_mat$Components, levels = paste0(rep("PC", 5), 1:5) %>% rev())
Corr_mat$variable <- factor(Corr_mat$variable, levels = c("Genes", "Reads", 'Timepoint','Treatment'))

var_names <- c("Timepoint" = "Timepoint", 
               "Treatment" = "Treatment",
               "Reads" = "Reads",
               'Genes'='Genes')

pdf("PC_correlation_pca_metadata_noOutlier_Colon.pdf", width = 7, height = 4)
ggplot(Corr_mat, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black", size = 3) +
  scale_x_discrete(labels = as_labeller(var_names)) +
  scale_size(range = c(1, 11), guide = "none") +
  scale_fill_gradient(low = "white", high = "#2171B5") +
  scale_color_gradient(low = "white", high = "#2171B5") +
  labs(fill = "Correlation", color = "Correlation", size = "Correlation", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        panel.border = element_rect(color = "grey"),
        axis.text = element_text(size = 11)) +
  coord_equal()
dev.off()

## selection of genes contributing to Treatment differences 
sig_Treatment<-subset(varPart_sorted, varPart_sorted$Treatment > 0.10)

#Remove genes effected by Timepoint, Genes or Reads
sig_Timepoint<-subset(varPart_sorted, varPart_sorted$Timepoint > 0.10)
sig_Reads<-subset(varPart_sorted, varPart_sorted$Reads > 0.10)
sig_Genes<-subset(varPart_sorted, varPart_sorted$Genes > 0.10)

remove<-rbind(sig_Timepoint,sig_Reads, sig_Genes )
remove_genes <- remove %>%
  distinct(ID) %>%
  pull(ID) 

Treatment_Genes_only<- subset(sig_Treatment, !rownames(sig_Treatment) %in% remove_genes)

## preform enrichment analyis
# GO enrichment
markers_entrez <- bitr(Treatment_Genes_only$ID,
                       fromType = "ENSEMBL",
                       toType="ENTREZID",
                       OrgDb=org.Mm.eg.db)
Treatment_Genes_only$ENSEMBL<-Treatment_Genes_only$ID
Treatment_Genes_only<-merge(Treatment_Genes_only, markers_entrez,by.x='ENSEMBL')

markers_entrez <- bitr(Treatment_Genes_only$ID,
                       fromType = "ENSEMBL",
                       toType="SYMBOL",
                       OrgDb=org.Mm.eg.db)
Treatment_Genes_only$ENSEMBL<-Treatment_Genes_only$ID
Treatment_Genes_only<-merge(Treatment_Genes_only, markers_entrez,by.x='ENSEMBL')

# is the gene up or downregulated in yMB?
dds_counts=DESeqDataSetFromMatrix(countData = count_data_Colon_treat_no, colData =Colon_treat_no, design = ~  Treatment )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
dds_norm=DESeq(dds_counts)
dds_normr=results(dds_norm)

write.csv(dds_normr, file="DE_Colon_8wVSageing_noOutlier.csv")

dds_normr$ENSEMBL<-row.names(dds_normr)
Down_comp<-as.data.frame(subset(dds_normr, log2FoldChange < 0))
Up_comp<-as.data.frame(subset(dds_normr, log2FoldChange > 0))

Down_comp$signal<-c('Upregulated in yMB')
Up_comp$signal<-c('Downregulated in yMB')
total<-rbind(Down_comp,Up_comp)

total_signal<- total[, c(7, 8 )]

Treatment_Genes_only_signal<-merge(Treatment_Genes_only, total_signal, by='ENSEMBL')

# Compare clusterProfiler
CompGO <- as.data.frame(compareCluster(ENTREZID~signal, data=Treatment_Genes_only_signal,OrgDb='org.Mm.eg.db', fun="enrichGO", ont = "BP"))

ENTREZID<-as.data.frame(CompGO$geneID)
nmax <- max(stringr::str_count(ENTREZID$`CompGO$geneID`, "\\/")) + 1
ENTREZID2<-tidyr::separate(ENTREZID, `CompGO$geneID`, paste0("col", seq_len(nmax)), sep = "\\/", fill = "right")
rows<- nrow(ENTREZID2)
GENEID <- as.data.frame(matrix(mapIds(org.Mm.eg.db, as.character(unlist(ENTREZID2)), "SYMBOL","ENTREZID"), rows))
col<-colnames(GENEID)
GENEID2<-unite(GENEID, col='GENEID', col, sep=',')
GENEID2$GENEID<-gsub("NULL,","",as.character(GENEID2$GENEID))
GENEID2$GENEID<-gsub(",NULL","",as.character(GENEID2$GENEID))
CompGO$GENEID<-GENEID2$GENEID
write.csv2(CompGO, file='GOComparison_Colon_Treatment_noOutlier_WithSignal.csv')

CompGO<-read.csv2('GOComparison_Colon_TreatmentNOTime_WithSignal_interest.csv')

CompGO %>% group_by(Cluster) %>% top_n(n= 10, wt = -p.adjust) -> terms
CompGOresults_sub<-subset(CompGO, CompGO$p.adjust <0.05)
CompGOresults_sub2<-subset(CompGOresults_sub, CompGOresults_sub$Count>=3)
#tmp <- CompGOresults_sub2[CompGOresults_sub2$Description %in% terms$Description,]
CompGOresults_sub2$Description <- ifelse(nchar(CompGOresults_sub2$Description)>60,
                                         paste(substr(CompGOresults_sub2$Description, 1, 60),"[...]",sep=""),
                                         CompGOresults_sub2$Description)
CompGOresults_sub2$Description <- factor(CompGOresults_sub2$Description,levels=unique(CompGOresults_sub2$Description))
CompGOresults_sub4<-CompGOresults_sub2 %>% group_by(Cluster) %>% top_n(n= 20, wt = -p.adjust)

CompGOresults_sub4$Cluster2 = with(CompGOresults_sub4, reorder(Cluster))


pdf("GOComparison_TreatmentnoOutlier_WithSignal_all_Colon.pdf", width = 8, height =6)
ggplot(CompGOresults_sub4, aes(x = p.adjust, y = Description,size = Count, fill=signal)) +
  geom_point(shape=21) +
  scale_fill_manual(values=c("cornflowerblue", "brown1"))+
  ylab(NULL) +
  theme_bw() +
  #facet_wrap(~ Cluster, ncol = 2, scales = "free") +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



