library(readxl)
library(tidyverse)
library(reshape2)

#go_genes <- read_excel("/Users/busch_lab/Dropbox/Shelby_Ch3/Analysis/working_data/Ch3_Analysis/modulespecific_METASCAPERESULTS/metascape_result.xlsx", sheet = "Annotation")
#go_genes <- read_excel("./modulespecific_METASCAPERESULTS/metascape_result.xlsx", sheet = "Annotation")

go_genes <- read.csv("top15.csv")

geneID <- read.delim("Bostr_GOannotations.txt", sep = "\t", header = T)

z <- geneID %>%
  separate_rows(GO, sep = ',')

anno_genes <- left_join(go_genes, z, by = c("TermID" = "GO"))









#circadian <- go_genes[go_genes$`GO:0007623 circadian rhythm` == "1",2]
#carotenoid <- go_genes[go_genes$`ath00906 Carotenoid biosynthesis - Arab` == "1",2]
#photosynthesis <- go_genes[go_genes$`GO:0019684 photosynthesis, light reaction` == "1",2]
#metabolites <- go_genes[go_genes$`ath00999 Biosynthesis of various plant` == "1",2]
#sesquit <- go_genes[go_genes$`GO:0006714 sesquiterpenoid metabolic proc` == "1",2]
#hydroxy <- go_genes[go_genes$`GO:1901617 organic hydroxy compound biosy` == "1",2]
#cellgrowth <- go_genes[go_genes$`GO:0009825 multidimensional cell growth` == "1",2]
#systemprocess <- go_genes[go_genes$`GO:0003008 system process` == "1",2]
#abioStress <- go_genes[go_genes$`GO:0071214 cellular response to abiotic s` == "1",2]
#nGlycan <- go_genes[go_genes$`ath00510 N-Glycan biosynthesis - Arabid` == "1",2]
#lightSTim <- go_genes[go_genes$`GO:0009416 response to light stimulus` == "1",2]

# read in z-scores datafame
z.scores <- read.delim("Bs_normalized.txt", sep = ",")
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

#inner join for each of the top 15 enriched categories and transpose 
elevation <- as.factor(c(2553, 2710, 2890, 3133,3342))

e <- inner_join(avg.Z, anno_genes, by = c("original_id"= "X.locusName"))
t_e <- setNames(data.frame(t(e[,-1])), e[,1])
t_e <- cbind(t_e, elevation)

d <- t_e
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("expression profiles")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#### circadian rhythm ####
circ_z <- inner_join(avg.Z, circadian, by = "original_id")
t_circ <- setNames(data.frame(t(circ_z[,-1])), circ_z[,1])
t_circ <- cbind(t_circ, elevation)

d <- t_circ
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

circ_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  ggtitle("circadian rhythm")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  
  

#### carotenoid biosynthesis ####
carot_z <- inner_join(avg.Z, carotenoid, by = "original_id")
t_caro <- setNames(data.frame(t(carot_z[,-1])), carot_z[,1])
t_caro <- cbind(t_caro, elevation)

d <- t_caro
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

carot_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("carotenoid biosynthesis")

#### photosynthesis ####
photo_z <- inner_join(avg.Z, photosynthesis, by = "original_id")
t_photo <- setNames(data.frame(t(photo_z[,-1])), photo_z[,1])
t_photo <- cbind(t_photo, elevation)

d <- t_photo
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

photo_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("photosynthesis")

#### secondary metabolism ####

metab_z <- inner_join(avg.Z, metabolites, by = "original_id")
t_metab <- setNames(data.frame(t(metab_z[,-1])), metab_z[,1])
t_metab <- cbind(t_metab, elevation)

d <- t_metab
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

metab_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("biosynthesis of secondary metabolites")

#### sesquiterpenoids ####
sesq_z <- inner_join(avg.Z, sesquit, by = "original_id")
t_sesq <- setNames(data.frame(t(sesq_z[,-1])), sesq_z[,1])
t_sesq <- cbind(t_sesq, elevation)

d <- t_sesq
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

sesq_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw() +
  ggtitle("biosynthesis of sesquiterpeniods")

#### hydroxy compounds ####
hydro_z <- inner_join(avg.Z, hydroxy, by = "original_id")
t_hydro <- setNames(data.frame(t(hydro_z[,-1])), hydro_z[,1])
t_hydro <- cbind(t_hydro, elevation)

d <- t_hydro
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

hydro_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("biosynthesis of organic hydroxy compounds")

#### cell growth ####
cell_z <- inner_join(avg.Z, cellgrowth, by = "original_id")
t_cell <- setNames(data.frame(t(cell_z[,-1])), cell_z[,1])
t_cell <- cbind(t_cell, elevation)

d <- t_cell
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

cell_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("multidimensional cell growth")

#### system process ####
syst_z <- inner_join(avg.Z, systemprocess, by = "original_id")
t_syst <- setNames(data.frame(t(syst_z[,-1])), syst_z[,1])
t_syst <- cbind(t_syst, elevation)

d <- t_syst
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

syst_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw() +
  ggtitle("system process")

#### abiotic stress response ####
abio_z <- inner_join(avg.Z, abioStress, by = "original_id")
t_abio <- setNames(data.frame(t(abio_z[,-1])), abio_z[,1])
t_abio <- cbind(t_abio, elevation)

d <- t_abio
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

abio_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw() +
  ggtitle("abiotic stress response")

#### n-glycan biosynthesis ####
nGly_z <- inner_join(avg.Z, nGlycan, by = "original_id")
t_nGly <- setNames(data.frame(t(nGly_z[,-1])), nGly_z[,1])
t_nGly <- cbind(t_nGly, elevation)

d <- t_nGly
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

nGly_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw()+
  ggtitle("n-Glycan biosynthesis")


#### response to radiation ####
lite_z <- inner_join(avg.Z, lightSTim, by = "original_id")
t_lite <- setNames(data.frame(t(lite_z[,-1])), lite_z[,1])
t_lite <- cbind(t_lite, elevation)

d <- t_lite
d <- melt(d, id.vars="elevation")
colnames(d)[2:3] <- c("gene ID", "z-score")

lite_plot <- ggplot(d, aes(elevation,`z-score`, col=`gene ID`)) + 
  geom_point() + 
  geom_line(aes(group = `gene ID`))+
  stat_smooth() +
  theme_bw() +
  ggtitle("response to radiation")

######

library(ggpubr)
ggarrange(circ_plot, carot_plot, photo_plot, metab_plot,
          sesq_plot, hydro_plot, cell_plot, syst_plot,
          abio_plot, nGly_plot,lite_plot, ncol = 6, nrow = 2)

library(gridExtra)
grid.arrange(circ_plot, carot_plot, photo_plot, metab_plot,
             sesq_plot, hydro_plot, cell_plot, syst_plot,
             abio_plot, nGly_plot,lite_plot,)


##### Big GO table ####
cr <- circadian
cr$go <- rep("Circadian rhythm")
car <- carotenoid
car$go <- rep("Carotenoid biosynthesis - Arabidopsis thaliana")
pho <- photosynthesis
pho$go <- rep("Photosynthesis, light reaction")
secmet <- metabolites
secmet$go <- rep("Biosynthesis of various plant secondary metabolites")
sesq <- sesquit
sesq$go <- rep("Sesquiterpenoid metabolic process")
hydro <- hydroxy
hydro$go <- rep("Organic hydroxy compound biosynthesis")
cell <- cellgrowth
cell$go <- rep("Multidimensional cell growth")
system$go <- rep("System process")
abio <- abioStress
abio$go <- rep("Cellular response to abiotic stress")
nGly <- nGlycan
nGly$go <- rep("N-Glycan biosynthesis - Arabidopsis thaliana")
light <- lightSTim
light$go <- rep("Response to light stimulus")

big_go <- rbind(cr, car, pho,secmet,sesq,hydro, cell, system, abio, nGly, light)

better_go <- aggregate(go~original_id,data=big_go, FUN = function(X) paste(unique(X), collapse=", "))

write(better_go, "better_go.csv")