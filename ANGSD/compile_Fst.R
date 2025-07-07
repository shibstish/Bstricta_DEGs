library(dplyr)
setwd("/scratch/dd66718/diversity/angsd/FST")

# Read in all *.print files of per site Fst values
file_list <- list.files(pattern = "*\\.print$")

# Combine all of the files into a single df
combined_data <- do.call(rbind, lapply(file_list, function(file) {
				 data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
				 colnames(data) <- c("chr", "pos", "A", "B")
				 return(data)
}))

# Calculate average weighted Fst for each chr pos locus
avg_fst <- combined_data %>% group_by(chr, pos) %>% summarize(average_fst = mean(A/B, na.rm = T))

write.table(avg_fst, "average_fst_locus.txt", row.names = FALSE, quote = FALSE, sep = "\t")

