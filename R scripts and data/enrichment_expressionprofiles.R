library(tidyverse)
library(readxl)
library(reshape2)
# read in the annotation file downloaded from https://phytozome-next.jgi.doe.gov/info/Bstricta_v1_2
anfile <- read.delim("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/Ch3_Analysis/Bstricta_278_v1.2.annotation_info.txt", sep ="\t")

# read in the datafile containing the list of genes that are significantly correlated with elevation/PC1
genes <- read.delim("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/Ch3_Analysis/module_DEGs_bstricta.txt", sep = "\t")

# join the two datafiles
at_goes <- inner_join(genes, anfile[,c(2,10:11)], by= c("geneid" = "locusName"), unmatched = "drop")

# extract the arabidopsis "best hits" 
at_genes <- at_goes$Best.hit.arabi.name

#write.csv(at_genes, "at_genes.csv", row.names = F) # read this file into Metascape.org


## after metascape ##

go_genes <- read_excel("./Metascape_Bstricta_05302025/metascape_result.xlsx", sheet = "Annotation")

# match the At geneIDs with the Boechera stricta geneIDs

bs_go <- as.data.frame(right_join(go_genes, at_goes, by = c("hits" = "Best.hit.arabi.name")))

# top 11 enriched categories
a <- as.data.frame(bs_go[bs_go$`GO:0007623 circadian rhythm` == "1.0",22])
b <- as.data.frame(bs_go[bs_go$`GO:0009416 response to light stimulus` == "1.0",22])
c <- as.data.frame(bs_go[bs_go$`GO:0009658 chloroplast organization` == "1.0",22])
d <- as.data.frame(bs_go[bs_go$`GO:0010228 vegetative to reproductive pha` == "1.0",22])
e <- as.data.frame(bs_go[bs_go$`GO:0071482 cellular response to light sti` == "1.0",22])
f <- as.data.frame(bs_go[bs_go$`GO:0016119 carotene metabolic process` == "1.0",22])
g <- as.data.frame(bs_go[bs_go$`GO:0009266 response to temperature stimul` == "1.0",22])
h <- as.data.frame(bs_go[bs_go$`GO:0010109 regulation of photosynthesis` == "1.0",22])
i <- as.data.frame(bs_go[bs_go$`GO:1901617 organic hydroxy compound biosy` == "1.0",22])
j <- as.data.frame(bs_go[bs_go$`GO:0006714 sesquiterpenoid metabolic proc` == "1.0",22])
k <- as.data.frame(bs_go[bs_go$`GO:0030258 lipid modification` == "1.0", 22])


names(a) <- "original_id"
names(b) <- "original_id"
names(c) <- "original_id"
names(d) <- "original_id"
names(e) <- "original_id"
names(f) <- "original_id"
names(g) <- "original_id"
names(h) <- "original_id"
names(i) <- "original_id"
names(j) <- "original_id"
names(k) <- "original_id"

# read in z-scores datafame
#normal.scores <- read.delim("Bs_normalized.txt", sep = ",")

z.scores <- read.csv("z.scores.csv")
z.scores$X <- gsub(".v1.2", "", z.scores$X)
colnames(z.scores)[1] <- "original_id"

# average z-score for each population
pops <- c("2253", "2710", "2890", "3133", "3342")
avg.Z <-data.frame(original_id = z.scores$original_id,
                   A = rowMeans(select(z.scores, A10, A26, A5)),
                   B = rowMeans(select(z.scores, B19, B28, B30)),
                   C = rowMeans(select(z.scores, C10, C11, C22)),
                   D = rowMeans(select(z.scores, D14, D16, D20)),
                   E = rowMeans(select(z.scores, E13, E7)))

elevation <- as.factor(c(2553, 2710, 2890, 3133,3342))

# make expression plots! top 4
#### circadian rhythm ####
a_z <- inner_join(avg.Z, a, by ="original_id")
t_a <- setNames(data.frame(t(a_z[,-1])), a_z[,1])
t_a <- cbind(t_a, elevation)

x <- t_a
x <- melt(x, id.vars="elevation")
colnames(x)[2:3] <- c("gene ID", "z-score")

a_plot <- ggplot(x, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("circadian rhythm")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

#### response to light stimulus ####
b_z <- inner_join(avg.Z, b, by ="original_id")
t_b <- setNames(data.frame(t(b_z[,-1])), b_z[,1])
t_b <- cbind(t_b, elevation)

x <- t_b
x <- melt(x, id.vars="elevation")
colnames(x)[2:3] <- c("gene ID", "z-score")

b_plot <- ggplot(x, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("response to light stimulus")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

#### chloroplast organization ####
c_z <- inner_join(avg.Z, c, by ="original_id")
t_c <- setNames(data.frame(t(c_z[,-1])), c_z[,1])
t_c <- cbind(t_c, elevation)

x <- t_c
x <- melt(x, id.vars="elevation")
colnames(x)[2:3] <- c("gene ID", "z-score")

c_plot <- ggplot(x, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("chloroplast organization")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

#### vegetative to reproductive phase transition of meristem ####
d_z <- inner_join(avg.Z, d, by ="original_id")
t_d <- setNames(data.frame(t(d_z[,-1])), d_z[,1])
t_d <- cbind(t_d, elevation)

x <- t_d
x <- melt(x, id.vars="elevation")
colnames(x)[2:3] <- c("gene ID", "z-score")

d_plot <- ggplot(x, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("vegetative to reproductive phase \ntransition of meristem")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

### compound plot! ####

library(ggpubr)

ggarrange(a_plot, b_plot, c_plot, d_plot, ncol = 2, nrow = 2)
