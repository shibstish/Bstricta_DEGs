library(adegenet)
library(ade4)
library(vcfR)
#install.packages("hierfstat")
library(hierfstat)
x <- read.vcfR("/Users/shelbyt/Library/CloudStorage/Dropbox/Shelby_Ch3/Analysis/working_data/bstricta_G.vcf")
y <- vcfR2genind(x)
pops <- as.factor(c(rep("2553", 3), rep("2710", 2), rep("2890", 2), rep("3133", 3), rep("3342", 2)))
z <- genind2hierfstat(y,pop=pops)

basic <- basic.stats(z)

per_site <- basic$perloc$Fst
write.csv(per_site, "hierFstat_sitewise_fst.csv")

overall_fst <- wc(z)
global <- overall_fst$FST
inbreeding <- overall_fst$FIS
fst <- data.frame(global = global, inbreeding = inbreeding)
write.csv(fst, "hierFstat_global_fst.csv")

pairwise_fst <- pairwise.neifst(z,diploid=TRUE)
rownames(pairwise_fst) <- c("2553", "2710", "2890", "3133", "3342")
colnames(pairwise_fst) <- c("2553", "2710", "2890", "3133", "3342")

write.csv(pairwise_fst, "hierFstat_pairwise_fst.csv")