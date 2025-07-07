# Load the commandArgs function to capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure that a file base name is provided
if (length(args) == 0) {
  stop("No input file provided. Usage: Rscript print_fst.R <file_base>")
}

# Capture the input file base from the argument
file_base <- args[1]

# Construct the full file name by appending ".print"
file_name <- paste0(file_base, ".print")

# Read the data (assuming tab-delimited)
data <- read.table(file_name, header = FALSE, sep = "\t")

# Calculate the ratio of sum of A (column 3) and B (column 4)
sum_unweighted <- sum(data$V3, na.rm = T)
sum_weighted <- sum(data$V4, na.rm = T)
ratio <- sum_unweighted / sum_weighted

# Calculate variance for A (column 3) and B (column 4)
var_unweighted <- var(data$V3, na.rm = T)
var_weighted <- var(data$V4, na.rm = T)

output_file <- paste0(file_base, ".txt")
sink(output_file)

cat("Ratio of sum of unweighted to weighted Fst: ", ratio, "\n")
cat("Variance of unweighted Fst: ", var_unweighted, "\n")
cat("Variance of weighted Fst: ", var_weighted, "\n")

sink()
