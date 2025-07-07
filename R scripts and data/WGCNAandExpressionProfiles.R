# ==================================
# WGCNA
# tutorial part 1: https://www.youtube.com/watch?v=gYE59uEMXT4
# tutorial part 2: https://www.youtube.com/watch?v=mzXIxjPr_Mc
# ==================================

library(WGCNA)
library(DESeq2)
library(tidyverse)
library(devtools)
#install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
library("biomartr")


#y <- read.table("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/featCounts_mx.txt", sep = "\t", header = T)
#y <- read.table("/Users/busch_lab/Dropbox/Shelby_Ch3/Analysis/working_data/featCounts_mx.txt", sep = "\t", header = T)
y <- read.table("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/fC_Bstricta.txt", sep = "\t", header = T)

genotype = c(rep("Bstricta", 42))
population = as.factor(c(rep("A", 9), rep("B", 9), rep("C", 9), rep("D", 9), rep("E", 6)))
band = c(rep("montane", 18),rep("subapline", 9), rep("alpine", 15))
#elevation <- as.factor(c(rep(2553, 9), rep(2710, 9), rep(2890, 9), rep(3133, 9), rep(3342, 6)))
aspect = c(rep(140.19, 9), rep(218.65, 9), rep(66.8, 9), rep(26.5, 9), rep(161.5, 6))
cardinal = c(rep("SE", 9), rep("W", 9), rep("NE", 18), rep("SW", 6))
slope = c(rep(11.6, 9), rep(2.8, 9), rep(6.3, 9), rep(14.9, 9), rep(16.4, 6))
water = as.factor(c(rep("moderate", 9), rep("low", 9), rep("high", 18), rep("moderate", 6)))


Bs_metadata <- data.frame(genotype, population, band, aspect, cardinal, slope, water)#,elevation)
rownames(Bs_metadata) <- c("A10_1", "A10_2", "A10_3", "A26_1", "A26_2", "A26_3", "A5_1", "A5_2", "A5_3",
                           "B19_1", "B19_2", "B19_3", "B28_1", "B28_2", "B28_3", "B30_1", "B30_2", "B30_3",
                           "C10_1", "C10_2", "C10_3", "C11_1", "C11_2", "C11_3", "C22_1", "C22_2", "C22_3",
                           "D14_1", "D14_2", "D14_3", "D16_1", "D16_2", "D16_3", "D20_1", "D20_2", "D20_3",
                           "E13_1", "E13_2", "E13_3", "E7_1", "E7_2", "E7_3") 

x <- y[,c(7:48)] # select only the columns with the count data


# set columns names in the data to match the row names in the meta data
colnames(x)<- c("A10_1", "A10_2", "A10_3", "A26_1", "A26_2", "A26_3", "A5_1", "A5_2", "A5_3",
                "B19_1", "B19_2", "B19_3", "B28_1", "B28_2", "B28_3", "B30_1", "B30_2", "B30_3",
                "C10_1", "C10_2", "C10_3", "C11_1", "C11_2", "C11_3", "C22_1", "C22_2", "C22_3",
                "D14_1", "D14_2", "D14_3", "D16_1", "D16_2", "D16_3", "D20_1", "D20_2", "D20_3",
                "E13_1", "E13_2", "E13_3", "E7_1", "E7_2", "E7_3")
rownames(x) <- y[,1]

# create DESeq object
dds_Bs <- DESeqDataSetFromMatrix(countData = x, colData = Bs_metadata, design = ~population)

# collapse technical replicates

dds_Bs$sample <- c("A10", "A10", "A10", "A26", "A26", "A26", "A5", "A5", "A5",
                   "B19", "B19", "B19", "B28", "B28", "B28", "B30", "B30", "B30",
                   "C10", "C10", "C10", "C11", "C11", "C11", "C22", "C22", "C22",
                   "D14", "D14", "D14", "D16", "D16", "D16", "D20", "D20", "D20",
                   "E13", "E13", "E13", "E7", "E7", "E7")

dds_Bs2 <- collapseReplicates(dds_Bs, dds_Bs$sample)

# use this metadata dataframe — this one has the correct dimensions for downstream math
col_geno <- c(rep("Bstricta", 14))
col_pop <- as.factor(c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3), rep("E", 2)))
col_band <- as.factor(c(rep("montane",6),rep("subapline", 3), rep("alpine", 5)))
col_elevation <- as.factor(c(rep(2553, 3), rep(2710, 3), rep(2890,3), rep(3133, 3), rep(3342, 2)))
col_aspect <- c(rep(140.19, 3), rep(218.65, 3), rep(66.8,3), rep(26.5, 3), rep(161.5, 2))
col_cardinal <- as.factor(c(rep("SE", 3), rep("W", 3), rep("NE", 6), rep("SW", 2)))
col_slope <- c(rep(11.6, 3), rep(2.8, 3), rep(6.3, 3), rep(14.9, 3), rep(16.4, 2))
col_water <- as.factor(c(rep(3, 3), rep(1, 3), rep(5, 6), rep(3, 2))) #1 = low; 3 = moderate; 5 = high
#col_solar <- c(rep(16.73, 3), rep(16.44, 3), rep(15.65, 3), rep(14.6, 5)) #collected on 07/26/2024
radiation <- c(rep(23.57, 3), rep(23.27, 3), rep(21.94, 3), rep(21.9, 3), rep(20.85, 2))#collected on 10/09/2024 (sloped solar irradiance, 800m resolution; June, July, and August values averaged. MJ/m2/day)
#col_mxvpd <- c(rep(12.09, 3), rep(11.25, 3), rep(9.46, 3), rep(8.75, 5)) #collected on 07/26/2024
mxvpd <- c(rep(23.22, 3), rep(22.16, 3), rep(20.16, 3), rep(17.00, 3), rep(14.58, 2))#collected on 10/09/2024 (800m resolution)
#col_temp <- c(rep(3.41, 3), rep(2.10, 3), rep(1.46, 3), rep(0.92, 5)) #collected on 07/26/2024
temp_800m <- c(rep(14.34, 3), rep(13.07, 3), rep(12.81, 3), rep(12.16, 3), rep(8.99, 2))#collected on 10/09/2024 

col_meta <- data.frame(col_geno, col_pop, col_band, col_elevation, col_aspect, col_cardinal, col_slope, col_water, radiation, mxvpd, temp_800m)
#col_meta <- data.frame(col_geno, col_pop,  col_aspect, col_slope, col_water, radiation, mxvpd, temp_800m)
rownames(col_meta) <- c("A10", "A26", "A5", "B19", "B28", "B30", "C10", "C11", "C22", "D14", "D16", "D20", "E13", "E7") 

# remove genes with counts < 15 in 75% of samples (14 samples * 0.75 = 10.5...rounding up)
# OG dataset
nrow(dds_Bs2) # 32833

dds75 <- dds_Bs2[rowSums(counts(dds_Bs2) >=15) >=11,]
nrow(dds75) # 16168

# normalize with function vst() from the DESeq2 package; transpose so that it's in the correct orientation for downstream analyses

dds_norm <- vst(dds75)

norm_counts <- assay(dds_norm) %>%
  t()

write.csv(norm_counts, "wgcna_normalized_bstricta.csv")

# ================================
# NETWORK CONTSTRUCTION
# ================================

# choose soft threshold
powers <- c(c(1:20), seq(from = 22, to=50, by=2))

sft <- pickSoftThreshold(norm_counts, powerVector = powers, networkType = "unsigned", verbose = 5)

# you want to find a power that maximizes R^2 and minimizes mean connectivity
sft.res <- sft$fitIndices

# visualizing the powers
# R^2
r2 <- ggplot(sft.res, aes(x = Power, y = SFT.R.sq, label = Power)) +
  geom_point()+
  geom_text(nudge_y = 0.03) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "scale free topology model fit, signed R^2") +
  ggtitle("Scale independence") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

# mean connectivity
meank <- ggplot(sft.res, aes(x = Power, y = mean.k., label = Power)) +
  geom_point()+
  geom_text( nudge_y = 1, nudge_x = 0.5) +
  labs(x = "Power", y = " mean connectivity") + 
  ggtitle("Mean connectivity") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

library(gridExtra)
grid.arrange(r2, meank, nrow = 2)
# hmmm maybe 10?

power <- sft$powerEstimate # 16.

# convert the matrix to numeric

norm_counts[] <- sapply(norm_counts, as.numeric)

# to avoid confusion later on, temporarily assign the
# built in cor() function to temp_cor and the WGCNA cor funtion to cor

temp_cor <- cor
cor <- WGCNA::cor

# using the blockwise approach because it is less computationally expensive
# (compared to the stepwise approach)

bwnet <- blockwiseModules(norm_counts,
                          maxBlockSize = 20000,
                          power = power,
                          TOMType = "unsigned",
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          set.seed(1234),
                          verbose = 5) #adjust blocksize based on RAM available - 16GB can handle about 20000; this will process all the genes in one block
# 47 modules found!

bwnet_signed <- blockwiseModules(norm_counts,
                          maxBlockSize = 20000,
                          power = power,
                          TOMType = "signed",
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          set.seed(1234),
                          verbose = 5)

# once the blockwise module function is done, you can reassign the temp_cor() function back to cor()


cor <- temp_cor

# ==============================
# DETERMINE MODULE EIGENGENES
# ==============================

module_eigengenes <- bwnet$MEs

head(module_eigengenes)

#write.csv(module_eigengenes, "module_eigengenes_Bstricta.csv")

# look to see the number of genes assigned to each module
genespermod <- table(bwnet$colors)
write.table(genespermod, "geneNumbersinMod_Bs.csv")


# plot the dendrogram with module colors both before and after merging

plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors), 
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.03)

plotDendroAndColors(bwnet$dendrograms[[1]], 
                    bwnet$colors,
                    "Module colors",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.03)

# neat.

#=====================
# correlations of WGNCA modules to environment
# using PC1 from the correlations_PCA.R script
#=====================
library(tidyverse)

module_eigengenes <- read.csv("module_eigengenes_Bstricta.csv")
module_eigengenes <-column_to_rownames(module_eigengenes, var = "X")

norm_counts <- read.csv("wgcna_normalized_bstricta.csv")
norm_counts <- column_to_rownames(norm_counts, var = "X")

nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)

pc1 <- read.csv("pca_eigenvals.csv")[,1:2]
pc1 <- column_to_rownames(pc1, var = "X")


library(CorLevelPlot)
library(WGCNA)

heatmap.data <- merge(module_eigengenes, pc1, by = "row.names")
heatmap.data <- heatmap.data %>%
  column_to_rownames(var = "Row.names")

colnames(heatmap.data) <- gsub("ME", "", colnames(heatmap.data))


PC1_plot <- CorLevelPlot(heatmap.data,
                         x = names(heatmap.data)[48],
                         y = names(heatmap.data)[1:47],
                         col = c("mediumpurple4", "mediumpurple1", "white",  "darkorange1", "darkorange4"),
                         cexLabX = 0,
                         fontLabY = 1,
                         titleX = "PC1",
                         titleY = "WGCNA module",
                         rotTitleY = 90,
                         posColKey = "right",
                         posLab = "none",
                         cexCorval = 1.25,
                         colCorval = "black") 

# different colors


PC1_plot <- CorLevelPlot(heatmap.data,
                         x = names(heatmap.data)[48],
                         y = names(heatmap.data)[1:47],
                         col =c("#006C00","#1F9622","#74B975", "#B0D8B0", "#E3F0E3",
                                "#FAE8F2","#F3BFDE","#E38FC3","#CD55A5","#A60B7F"),
                        cexLabX = 0,
                         fontLabY = 1,
                         titleX = "PC1",
                         titleY = "WGCNA module",
                         rotTitleY = 90,
                         posColKey = "right",
                         posLab = "none",
                         cexCorval = 1.25,
                         colCorval = "black") 


## add colored squares to symbolize the modules
library(ggplot2)

colors <- names(heatmap.data)[1:47]

# Create sample data for square coordinates and colors
data <- data.frame(
  x = c(1),
  y = c(1:47),
  color = colors
)

# Create the plot using geom_rect()
ggplot(data, aes(xmin = x - 0.1, xmax = x + 0.1 , ymin = y -0.5, ymax = y+0.5, fill = color)) + 
  geom_rect() +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# correlation dataframe
z <- cor(module_eigengenes, pc1)
z_p<- corPvalueStudent(z, nSamples)


colnames(z) <- c("PC1_r")
colnames(z_p) <- c("PC1_pvals")
elevation.stats <- merge(z, z_p, by=0, all = T)
elevation.stats <- elevation.stats %>% column_to_rownames("Row.names")
siggies <- c(which(elevation.stats$PC1_pvals < 0.01))
siggies <- sort(unique(siggies)) 
elevation.stats <- elevation.stats %>% rownames_to_column(var = "moduleColor")
elevation.colors <- elevation.stats$moduleColor[siggies]
elevation.mod <- gsub("ME", "", elevation.colors)

r_p <- elevation.stats[siggies,]
rownames(r_p) <- r_p$moduleColor
r_p <- r_p[,2:3]

rownames(r_p) <-gsub("ME", "", rownames(r_p))

r_p <- round(r_p, 4)
colnames(r_p) <- c("weighted Pearson correlation (r)", "p-value")

pos.cor <- r_p[c(which(r_p$`weighted Pearson correlation (r)` > 0)),]
negative.cor <- r_p[c(which(r_p$`weighted Pearson correlation (r)` < 0)),]

write.csv(r_p, "wgcna_moduleCorrelationresults_bstricta.csv")

#=============================
# Extract genes from modules of interest!
#=============================

# the gene-module mapping is in bwnet$colors; 
# extract this out into its own dataframe and then we can subset out the modules we're interested in

allgenes <- as.data.frame(bwnet$colors)

allgenes <- rownames_to_column(allgenes, var = "geneID")
allgenes$geneID <- gsub("gene.", "", allgenes$geneID)
#write.csv(allgenes, "moduleMappedGenes_bstricta.csv")

for(i in elevation.mod){
  mod <- paste0("mod_",i)
  assign(mod, allgenes %>% filter(bwnet$colors == i))
}

#if reading in the csv file use:

allgenes <- read.csv("moduleMappedGenes_bstricta.csv")[,2:3]

for(i in elevation.mod){
  mod <- paste0("mod_",i)
  assign(mod, allgenes %>% filter(bwnet.colors == i))
}

# universal
all_gene_mods <- data.frame(NULL)

dlist <- lapply ( ls(patt='mod_'), get)

for (i in 1:length(dlist)) {
  temp <- dlist[[i]]
  all_gene_mods <- rbind(all_gene_mods, temp)
}

#all_gene_mods <- rownames_to_column(all_gene_mods, var = "geneID")

mod.list <- vector("list", length = length(elevation.mod))

for(i in elevation.mod){
  mod <-  allgenes %>% filter(bwnet$colors == i)
  mod.list[[i]] <- mod
}

mod.genes <- do.call(rbind, mod.list)

degs <- as.data.frame(read.csv("DEGs_Bstricta.csv")[,2])

colnames(degs) <- "geneID"

degs_mods <- inner_join(mod.genes, degs, by = "geneID")
degs_mods$geneID <- gsub(".v1.2", "", degs_mods$geneID)

write.csv(degs_mods, "module_DEGs_bstricta.csv", row.names = FALSE)

#use this file for functional enrichment analysis :)

#### Additional/optional stuff ####
exclude <- filter(allgenes, !`bwnet$colors` %in% rownames(r_p))
hubGenes <- chooseTopHubInEachModule(norm_counts, bwnet$colors, omitColors = exclude$`bwnet$colors`, power = power, type = "unsigned")

hubGenes <- as.data.frame(hubGenes) %>% rownames_to_column(var = "module")
hubGenes$hubGenes <- gsub("gene:", "", hubGenes$hubGenes)


oneHub <- chooseOneHubInEachModule(norm_counts, bwnet$colors, omitColors = exclude$`bwnet$colors`, power = power, type = "unsigned")

#=====================
# plotting the "eigengene" expression profiles of the PC1-correlated modules
#======================

genes=t(module_eigengenes)

A <- genes[,1:3] %>% as.data.frame() %>% mutate_at(c('A10', 'A26', 'A5'), as.numeric)
B <- genes[,4:6] %>% as.data.frame() %>% mutate_at(c('B19', 'B28', 'B30'), as.numeric)
C <- genes[,7:9] %>% as.data.frame() %>% mutate_at(c('C10', 'C11', 'C22'), as.numeric)
D <- genes[,10:12] %>% as.data.frame() %>% mutate_at(c('D14', 'D16', 'D20'), as.numeric)
E <- genes[,13:14] %>% as.data.frame() %>% mutate_at(c('E13', 'E7'), as.numeric)

A <- as.data.frame(rowMeans(A, na.rm = T))
colnames(A) <- "A"
B <- as.data.frame(rowMeans(B, na.rm = T))
colnames(B) <- "B"
C <- as.data.frame(rowMeans(C, na.rm = T))
colnames(C) <- "C"
D <- as.data.frame(rowMeans(D, na.rm = T))
colnames(D) <- "D"
E <- as.data.frame(rowMeans(E, na.rm = T))
colnames(E) <- "E"

pops <- cbind(A, cbind(B, cbind(C, cbind(D, E))))
rownames(pops) <- gsub("ME", "", rownames(pops))

ele.big <- as.data.frame(pops)
elevation.mod <- as.array(elevation.mod)
ele.big <- ele.big[rownames(ele.big) %in% elevation.mod,]
ele.big <- as.data.frame(t(ele.big))
ele.big2 <- gather(ele.big, key = "module", value = "expression")
ele.big2$population <- rep(c("A", "B", "C", "D", "E"), 9)

env = read.csv('raw_metadata.csv', header=T, stringsAsFactors = F)
env = env[,c(2,5:13)]
env = env[c(1,4,7,10,13),]
env[,1]=c("A","B","C","D","E")

n <-9 
env2 <- do.call("rbind", replicate( 
  n, env, simplify = FALSE)) 

ele.big2=cbind(ele.big2, env2)
ele.big2 = group_by(ele.big2, module)

ele.colors <- c(cyan = "cyan", tan = "tan", steelblue = "steelblue", 
                darkslateblue= "darkslateblue", lightyellow = "gold",
                blue = "blue", green = "green", lightcyan1 = "lightcyan1")

elev <- ggplot(data = ele.big2[ele.big2$module=="cyan",], aes(x = col_elevation, y = expression)) +
  geom_line(color = "cyan")+
  geom_line(data = ele.big2[ele.big2$module=="tan",], color = "tan")+
  geom_line(data = ele.big2[ele.big2$module=="steelblue",], color = "steelblue")+
  geom_line(data = ele.big2[ele.big2$module=="lightcyan1",], color = "lightcyan3")+
  geom_line(data = ele.big2[ele.big2$module=="darkslateblue",], color = "darkslateblue") + 
  geom_line(data = ele.big2[ele.big2$module=="lightyellow",], color = "gold")+ 
  geom_line(data = ele.big2[ele.big2$module=="blue",], color = "blue")+ 
  geom_line(data = ele.big2[ele.big2$module=="green",], color = "green")+
  xlab("Elevation (m)")+ylab("Module expression")+ theme_bw()+theme(legend.position="none")

