library(BiocManager)
library(GEOquery)
library(limma)
install.packages("RobustRankAggreg")
library(RobustRankAggreg)


# Read the configuration file containing paths to top1000 files
working_dir <- "C:\\Users\\itaia\\OneDrive\\Project_TMED"
setwd(working_dir)
current_dir <- getwd()
print(current_dir)
config_file <- "Configs_top1000.txt"
top1000_paths <- read.delim(config_file, header = FALSE, stringsAsFactors = FALSE)

# Create an empty data frame to store the combined results
# combined_table <- data.frame(matrix(ncol = 0, nrow = 0))

de_results <- list()
# Loop through each top1000 file and combine the data
for (i in 1:length(top1000_paths)) {
  # Read the top1000 file
  top1000_file <- top1000_paths[[i]]
  print(top1000_file)
  # Create a temporary list to store the data for each file
  temp_list <- list()
  
  # Read each top1000 file and store the data
  for (j in 1:length(top1000_file)) {
    temp_file <- top1000_file[j]
    top1000_data <- read.delim(temp_file, header = TRUE, stringsAsFactors = FALSE)
    print(top1000_data)
    de_results[[j]] <- top1000_data
  }
}
# Combine the DEG results from each dataset into a single matrix
glist <- do.call(cbind, de_results)

result <- aggregateRanks(glist = de_results, N = 1000, method = "stuart")
# Run RobustRankAggreg on the DEG matrix
rra_result <- RRA(de_matrix)

# Extract the integrated DEGs based on the result
integrated_degs <- rownames(rra_result$aggregatedRank)
print(integrated_degs)






write.table(combined_table, file ="dup_top1000_combined_table.tsv", row.names = FALSE, sep = "\t")
# Remove duplicated genes, if any
combined_table <- unique(combined_table)

# Print the combined table to check the results
print(combined_table)
write.table(combined_table, file ="top1000_combined_table.tsv", row.names = FALSE, sep = "\t")


##Ranking
# Compute the maximum AveExpr
max_ave_expr <- max(combined_table$AveExpr)
# Apply the significance cutoff values
filtered_table <- subset(combined_table, adj.P.Val < 0.1 & P.Value <= 0.1 & (AveExpr/max_ave_expr)*100 >= 30 & abs(logFC) > 1)

# Rank the genes based on the statistical parameters
ranked_table <- filtered_table[order(filtered_table$adj.P.Val), ]

# Select the top 100 significant genes
top_50_genes <- head(ranked_table, 50)

# Print the top 50 significant genes
print(top_50_genes)
write.table(combined_table, file ="combined_top50.tsv", row.names = FALSE, sep = "\t")
