library("DESeq2")
library("pheatmap")
library("tidyverse")
library("dplyr")
library("tidyr")
library("stringr")
library("biomartr")
library("data.table")

#==========================================
# DESeq2 analysis
#==========================================

# read in the feature counts matrix; adjust path as necessary

y <- read.table("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/fC_Bstricta.txt", sep = "\t", header = T)
x <- y[,c(7:48)] # select only the columns with the count data


# set columns names in the data to match the row names in the meta data
colnames(x)<- c("A10_1", "A10_2", "A10_3", "A26_1", "A26_2", "A26_3", "A5_1", "A5_2", "A5_3",
                "B19_1", "B19_2", "B19_3", "B28_1", "B28_2", "B28_3", "B30_1", "B30_2", "B30_3",
                "C10_1", "C10_2", "C10_3", "C11_1", "C11_2", "C11_3", "C22_1", "C22_2", "C22_3",
                "D14_1", "D14_2", "D14_3", "D16_1", "D16_2", "D16_3", "D20_1", "D20_2", "D20_3",
                "E13_1", "E13_2", "E13_3", "E7_1", "E7_2", "E7_3")
rownames(x) <- y[,1]

# create a prelimiary metadata file

genotype = c(rep("Bstricta", 42))
population = as.factor(c(rep("2553", 9), rep("2710", 9), rep("2890", 9), rep("3133", 9), rep("3342", 6)))
band = c(rep("montane", 18),rep("subapline", 9), rep("alpine", 15))
elevation <- c(rep(2553, 9), rep(2710, 9), rep(2890, 9), rep(3133, 9), rep(3342, 6))
aspect = c(rep(140.19, 9), rep(218.65, 9), rep(66.8, 9), rep(26.5, 9), rep(161.5, 6))
cardinal = c(rep("SE", 9), rep("W", 9), rep("NE", 18), rep("SW", 6))
slope = c(rep(11.6, 9), rep(2.8, 9), rep(6.3, 9), rep(14.9, 9), rep(16.4, 6))
water = as.factor(c(rep("moderate", 9), rep("low", 9), rep("high", 18), rep("moderate", 6)))

Bs_metadata <- data.frame(genotype, population, band, elevation, aspect, cardinal, slope, water)
rownames(Bs_metadata) <- c("A10_1", "A10_2", "A10_3", "A26_1", "A26_2", "A26_3", "A5_1", "A5_2", "A5_3",
                           "B19_1", "B19_2", "B19_3", "B28_1", "B28_2", "B28_3", "B30_1", "B30_2", "B30_3",
                           "C10_1", "C10_2", "C10_3", "C11_1", "C11_2", "C11_3", "C22_1", "C22_2", "C22_3",
                           "D14_1", "D14_2", "D14_3", "D16_1", "D16_2", "D16_3", "D20_1", "D20_2", "D20_3",
                           "E13_1", "E13_2", "E13_3", "E7_1", "E7_2", "E7_3") 



# create DESeq object
dds_Bs <- DESeqDataSetFromMatrix(countData = x, colData = Bs_metadata, design = ~population)

# collapse technical replicates

dds_Bs$sample <- c("A10", "A10", "A10", "A26", "A26", "A26", "A5", "A5", "A5",
                   "B19", "B19", "B19", "B28", "B28", "B28", "B30", "B30", "B30",
                   "C10", "C10", "C10", "C11", "C11", "C11", "C22", "C22", "C22",
                   "D14", "D14", "D14", "D16", "D16", "D16", "D20", "D20", "D20",
                   "E13", "E13", "E13", "E7", "E7", "E7")

dds_Bs <- collapseReplicates(dds_Bs, dds_Bs$sample)

# use this metadata dataframe — this one has the correct dimensions for downstream math
col_geno <- c(rep("Bstricta", 14))
col_pop <- as.factor(c(rep("2553", 3), rep("2710", 3), rep("2890", 3), rep("3133", 3), rep("3342", 2)))
col_band <- as.factor(c(rep("montane",6),rep("subapline", 3), rep("alpine", 5)))
col_elevation <- as.factor(c(rep(2553, 3), rep(2710, 3), rep(2890,3), rep(3133, 3), rep(3342, 2)))
col_aspect <- c(rep(140.19, 3), rep(218.65, 3), rep(66.8,3), rep(26.5, 3), rep(161.5, 2))
col_cardinal <- as.factor(c(rep("SE", 3), rep("W", 3), rep("NE", 6), rep("SW", 2)))
col_slope <- c(rep(11.6, 3), rep(2.8, 3), rep(6.3, 3), rep(14.9, 3), rep(16.4, 2))
col_water <- as.factor(c(rep(3, 3), rep(1, 3), rep(5, 6), rep(3, 2))) #1 = low; 3 = moderate; 5 = high
col_solar <- c(rep(16.73, 3), rep(16.44, 3), rep(15.65, 3), rep(14.6, 5)) #collected on 07/26/2024
col_mxvpd <- c(rep(12.09, 3), rep(11.25, 3), rep(9.46, 3), rep(8.75, 5)) #collected on 07/26/2024
col_temp <- c(rep(3.41, 3), rep(2.10, 3), rep(1.46, 3), rep(0.92, 5)) #collected on 07/26/2024

col_meta <- data.frame(col_geno, col_pop, col_band, col_elevation, col_aspect, col_cardinal, col_slope, col_water, col_solar, col_mxvpd, col_temp)
rownames(col_meta) <- c("A10", "A26", "A5", "B19", "B28", "B30", "C10", "C11", "C22", "D14", "D16", "D20", "E13", "E7") 

# remove genes with counts < 15 in 75% of samples (14 samples * 0.75 = 10.5...rounding up) to reduce variability
# OG dataset
nrow(dds_Bs) # 27416

dds75 <- dds_Bs[rowSums(counts(dds_Bs) >=15) >=11,]
nrow(dds75) # 19290

dds75 <- estimateSizeFactors(dds75)
sizeFactors(dds75)

# normalize to read counts per million (TPM)

Bs_normalized <- counts(dds75, normalized = TRUE)
Bs_normalized$z_scores <- scale(assay(vsd))

#write.csv(Bs_normalized, "Bs_normalized.txt")
# log transform
# DESeq2 uses variance stabilizing transformation (VST)

vsd_Bs <- vst(dds75, blind = TRUE)

# extract the matrix
vsd_mat_Bs <- assay(vsd_Bs) #transformed values
#calculate z-scores
z_scores <- scale(assay(vsd_Bs))
#write.csv(z_scores, "z.scores.csv")

# calculate the pairwise correlation values
vsd_cor_Bs <- cor(vsd_mat_Bs)

# heatmap!

library(RColorBrewer)
library(viridis)
 
heat.colors <- viridis_pal(option = "D")
df <- as.data.frame(colData(dds75) [,c("population", "sample")])

pheatmap(vsd_cor_Bs, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, color = viridis::mako(n = 5), main = "pairwise correlation of variance stablized gene expression values", annotation_col = df)


# PCA
plotPCA(vsd_Bs, intgroup = c("population")) +  scale_color_manual(values = c("#FF3300","#FF9933", "#660066",  "#669900", "#9999FF"), name = "Population") +
  theme_minimal() + theme(axis.title.x = element_text(size = 12), 
                          axis.title.y = element_text(size = 12), 
                          axis.text.x = element_text(size = 12), 
                          axis.text.y = element_text(size = 12), 
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 12))

#==============
# run DE analysis
#==============
# use the fit that has the smallest absolute median residual; this is the best test to use.
#DE_obj <- DESeq(dds75, test = "Wald", fitType = "mean")
#DE_obj <- DESeq(dds75, test = "Wald", fitType = "parametric")
DE_obj <- DESeq(dds75, test = "Wald", fitType = "local") # this one!

residual <- log(mcols(DE_obj)$dispGeneEst) - log(mcols(DE_obj)$dispFit)
abs(median(residual)) #smallest absolute median residual of all the tests

# check dispersion estimates
plotDispEsts(DE_obj)

# results for the pairwise populations
results_ab <- results(DE_obj, contrast = c("population", "2553", "2710"))
results_ac <- results(DE_obj, contrast = c("population", "2553", "2890"))
results_ad <- results(DE_obj, contrast = c("population", "2553", "3133"))
results_ae <- results(DE_obj, contrast = c("population", "2553", "3342"))
results_bc <- results(DE_obj, contrast = c("population", "2710", "2890"))
results_bd <- results(DE_obj, contrast = c("population", "2710", "3133"))
results_be <- results(DE_obj, contrast = c("population", "2710", "3342"))
results_cd <- results(DE_obj, contrast = c("population", "2890", "3133"))
results_ce <- results(DE_obj, contrast = c("population", "2890", "3342"))
results_de <- results(DE_obj, contrast = c("population", "3133", "3342"))

# plot MA plots
par(mfrow=c(2,5))
plotMA(results_ab, ylim=c(-8,8), alpha = 0.01, main = "A to B", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_ac, ylim=c(-8,8), alpha = 0.01, main = "A to C", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_ad, ylim=c(-8,8), alpha = 0.01, main = "A to D", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_ae, ylim=c(-8,8), alpha = 0.01, main = "A to E", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_bc, ylim=c(-8,8), alpha = 0.01, main = "B to C", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_bd, ylim=c(-8,8), alpha = 0.01, main = "B to D", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_be, ylim=c(-8,8), alpha = 0.01, main = "B to E", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_cd, ylim=c(-8,8), alpha = 0.01, main = "C to D", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_ce, ylim=c(-8,8), alpha = 0.01, main = "C to E", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)
plotMA(results_de, ylim=c(-8,8), alpha = 0.01, main = "D to E", colNonSig = "lightsteelblue", colSig = "mediumvioletred", cex = .75)


# Arranging by adjusted p-value
df_A_B <- data.frame(results_ab)
df_A_C <- data.frame(results_ac)
df_A_D <- data.frame(results_ad)
df_A_E <- data.frame(results_ae)
df_B_C <- data.frame(results_bc)
df_B_D <- data.frame(results_bd)
df_B_E <- data.frame(results_be)
df_C_D <- data.frame(results_cd)
df_C_E <- data.frame(results_ce)
df_D_E <- data.frame(results_de)

rownames(df_A_B) <- gsub("gene:", "", rownames(df_A_B))
rownames(df_A_C) <- gsub("gene:", "", rownames(df_A_C))
rownames(df_A_D) <- gsub("gene:", "", rownames(df_A_D))
rownames(df_A_E) <- gsub("gene:", "", rownames(df_A_E))
rownames(df_B_C) <- gsub("gene:", "", rownames(df_B_C))
rownames(df_B_D) <- gsub("gene:", "", rownames(df_B_D))
rownames(df_B_E) <- gsub("gene:", "", rownames(df_B_E))
rownames(df_C_D) <- gsub("gene:", "", rownames(df_C_D))
rownames(df_C_E) <- gsub("gene:", "", rownames(df_C_E))
rownames(df_D_E) <- gsub("gene:", "", rownames(df_D_E))

A_B_sig <- subset(df_A_B, padj < 0.01)
A_C_sig <- subset(df_A_C, padj < 0.01)
A_D_sig <- subset(df_A_D, padj < 0.01)
A_E_sig <- subset(df_A_E, padj < 0.01)
B_C_sig <- subset(df_B_C, padj < 0.01)
B_D_sig <- subset(df_B_D, padj < 0.01)
B_E_sig <- subset(df_B_E, padj < 0.01)
C_D_sig <- subset(df_C_D, padj < 0.01)
C_E_sig <- subset(df_C_E, padj < 0.01)
D_E_sig <- subset(df_D_E, padj < 0.01)

A_B_sig <- A_B_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
A_C_sig <- A_C_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
A_D_sig <- A_D_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
A_E_sig <- A_E_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
B_C_sig <- B_C_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
B_D_sig <- B_D_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
B_E_sig <- B_E_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
C_D_sig <- C_D_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
C_E_sig <- C_E_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")
D_E_sig <- D_E_sig %>% arrange(padj) %>% rownames_to_column(var = "ensembl_gene_id")



degs <- c(A_B_sig$ensembl_gene_id, A_C_sig$ensembl_gene_id, A_D_sig$ensembl_gene_id, A_E_sig$ensembl_gene_id, 
              B_C_sig$ensembl_gene_id, B_D_sig$ensembl_gene_id, B_E_sig$ensembl_gene_id, 
              C_D_sig$ensembl_gene_id, C_E_sig$ensembl_gene_id, 
              D_E_sig$ensembl_gene_id) 

degs <- unique(degs) # length = 775

write.csv(degs, "DEGs.csv", row.names = FALSE)

# determining number of up/down regulated genes in each pairwise comparision
nrow(A_B_sig) #152
length(which(A_B_sig$log2FoldChange < 0)) # 52
length(which(A_B_sig$log2FoldChange > 0)) # 100

nrow(A_C_sig) #207
length(which(A_C_sig$log2FoldChange < 0)) # 70
length(which(A_C_sig$log2FoldChange > 0)) # 137

nrow(A_D_sig) #222
length(which(A_D_sig$log2FoldChange < 0)) # 54
length(which(A_D_sig$log2FoldChange > 0)) # 168

nrow(A_E_sig) #734
length(which(A_E_sig$log2FoldChange < 0)) # 239
length(which(A_E_sig$log2FoldChange > 0)) # 495

nrow(B_C_sig) #165
length(which(B_C_sig$log2FoldChange < 0)) # 81
length(which(B_C_sig$log2FoldChange > 0)) # 84

nrow(B_D_sig) #177
length(which(B_D_sig$log2FoldChange < 0)) # 80
length(which(B_D_sig$log2FoldChange > 0)) # 97

nrow(B_E_sig) #437
length(which(B_E_sig$log2FoldChange < 0)) # 178
length(which(B_E_sig$log2FoldChange > 0)) # 259

nrow(C_D_sig) #140
length(which(C_D_sig$log2FoldChange < 0)) # 60
length(which(C_D_sig$log2FoldChange > 0)) # 80

nrow(C_E_sig) #253
length(which(C_E_sig$log2FoldChange < 0)) # 11
length(which(C_E_sig$log2FoldChange > 0)) # 142

nrow(D_E_sig) #333
length(which(D_E_sig$log2FoldChange < 0)) # 141
length(which(D_E_sig$log2FoldChange > 0)) # 192

#==========================
## visualizing the results

#### reformat the results as data frames ####
ab_df <- as.data.frame(results_ab)
ab_df <- ab_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

ac_df <- as.data.frame(results_ac)
ac_df <- ac_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

ad_df <- as.data.frame(results_ad)
ad_df <- ad_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

ae_df <- as.data.frame(results_ae)
ae_df <- ae_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

bc_df <- as.data.frame(results_bc)
bc_df <- bc_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

bd_df <- as.data.frame(results_bd)
bd_df <- bd_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

be_df <- as.data.frame(results_be)
be_df <- be_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

cd_df <- as.data.frame(results_cd)
cd_df <- cd_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

ce_df <- as.data.frame(results_ce)
ce_df <- ce_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

de_df <- as.data.frame(results_de)
de_df <- de_df %>%
  rownames_to_column(var = "geneID") %>%
  mutate(threshold = padj < 0.01)

#### volcano plots ####
ab_volcano <- ggplot(ab_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2553 and 2710") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

ac_volcano <- ggplot(ac_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2553 and 2890") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

ad_volcano <- ggplot(ad_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2553 and 3133") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

ae_volcano <- ggplot(ae_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2553 and 3342") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

bc_volcano <- ggplot(bc_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2710 and 2890") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

bd_volcano <- ggplot(bd_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2710 and 3133") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

be_volcano <- ggplot(be_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2710 and 3342") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

cd_volcano <- ggplot(cd_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2890 and 3133") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

ce_volcano <- ggplot(ce_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("2890 and 3342") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

de_volcano <- ggplot(de_df) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ylim(c(0,15)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("3133 and 3342") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15)))

#### altogether now ####
library(gridExtra)

grid.arrange(ab_volcano, ac_volcano, ad_volcano, ae_volcano,
             bc_volcano, bd_volcano, be_volcano,
             cd_volcano, ce_volcano,
             de_volcano, ncol = 5)

