### script for Supplemental Table 1 ###

library(tidyverse)

anfile <- read.delim("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/Ch3_Analysis/Bstricta_278_v1.2.annotation_info.txt", sep ="\t")

# read in the datafile containing the list of genes that are significantly correlated with elevation/PC1
genes <- read.delim("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/Ch3_Analysis/module_DEGs_bstricta.txt", sep = "\t")

# join the two datafiles
at_goes <- inner_join(genes, anfile[,c(2,10:11)], by= c("geneid" = "locusName"), unmatched = "drop")

go_genes <- read_excel("./Metascape_Bstricta_05302025/metascape_result.xlsx", sheet = "Annotation")

# match the At geneIDs with the Boechera stricta geneIDs

bs_go <- as.data.frame(right_join(go_genes, at_goes, by = c("hits" = "Best.hit.arabi.name")))

# top 11 enriched categories
a <- as.data.frame(bs_go[bs_go$`GO:0007623 circadian rhythm` == "1.0",c(22, 1, 5, 23)])
b <- as.data.frame(bs_go[bs_go$`GO:0009416 response to light stimulus` == "1.0",c(22, 1, 5, 23)])
c <- as.data.frame(bs_go[bs_go$`GO:0009658 chloroplast organization` == "1.0",c(22, 1, 5, 23)])
d <- as.data.frame(bs_go[bs_go$`GO:0010228 vegetative to reproductive pha` == "1.0",c(22, 1, 5, 23)])
e <- as.data.frame(bs_go[bs_go$`GO:0071482 cellular response to light sti` == "1.0",c(22, 1, 5, 23)])
f <- as.data.frame(bs_go[bs_go$`GO:0016119 carotene metabolic process` == "1.0",c(22, 1, 5, 23)])
g <- as.data.frame(bs_go[bs_go$`GO:0009266 response to temperature stimul` == "1.0",c(22, 1, 5, 23)])
h <- as.data.frame(bs_go[bs_go$`GO:0010109 regulation of photosynthesis` == "1.0",c(22, 1, 5, 23)])
i <- as.data.frame(bs_go[bs_go$`GO:1901617 organic hydroxy compound biosy` == "1.0",c(22, 1, 5, 23)])
j <- as.data.frame(bs_go[bs_go$`GO:0006714 sesquiterpenoid metabolic proc` == "1.0",c(22, 1, 5, 23)])
k <- as.data.frame(bs_go[bs_go$`GO:0030258 lipid modification` == "1.0", c(22, 1, 5, 23)])

write.csv(a, "a.csv")
write.csv(b, "b.csv")
write.csv(c, "c.csv")
write.csv(d, "d.csv")
write.csv(e, "e.csv")
write.csv(f, "f.csv")
write.csv(g, "g.csv")
write.csv(h, "h.csv")
write.csv(i, "i.csv")
write.csv(j, "j.csv")
write.csv(k, "k.csv")
