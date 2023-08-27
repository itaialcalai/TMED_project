# full analysis of the data
library(BiocManager)
library(GEOquery)
library(limma)

run_study_analysis <- function(study_name, study_group_string, study_format, study_num, study_country, study_tissue) {
  tryCatch(
    {
    gset <- getGEO(study_name, GSEMatrix = TRUE, AnnotGPL = TRUE)
    if (length(gset) > 1) idx <- grep(study_format, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    
    # Group membership for all samples
    gsms <- study_group_string
    sml <- strsplit(gsms, split = "")[[1]]
    
    # filter out excluded samples (marked as "X")
    sel <- which(sml != "X")
    sml <- sml[sel]
    gset <- gset[ ,sel]
    
    # log2 transformation
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
    LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    if (LogC) {
      ex[which(ex <= 0)] <- NaN
      exprs(gset) <- log2(ex)
    }
    
    # assign samples to groups and set up the design matrix
    gs <- factor(sml)
    groups <- make.names(c("T2D", "normal"))
    levels(gs) <- groups
    gset$group <- gs
    design <- model.matrix(~group + 0, gset)
    colnames(design) <- levels(gs)
    
    gset <- gset[complete.cases(exprs(gset)), ]  # Skip missing values
    
    fit <- lmFit(gset, design)  # Fit linear model
    
    # Set up contrasts of interest and recalculate model coefficients
    cts <- paste(groups[1], groups[2], sep = "-")
    cont.matrix <- makeContrasts(contrasts = cts, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    
    # Compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
    
    # check columns
    colnames(tT)
    # Subset and write the table
    tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "AveExpr", "t", "B", "logFC", "Gene.symbol", "Gene.ID", "Gene.title"))
    write.table(tT, file = paste0(study_name, "_output_table.tsv"), row.names = FALSE, sep = "\t")
    
    # Write the top 1000 genes to a separate file
    tT_top1000 <- head(tT, 1000)
    write.table(tT_top1000, file = paste0(study_name, "_top1000_genes.tsv"), row.names = FALSE, sep = "\t")
    
    # Define the genes to be plotted
    control_genes <- c("TMED9", "PDE4B", "NFKB1", "SRR", "GAPDH")
    
    # Subset the table for the genes of interest, including NaN values and duplicates
    genes_subset <- subset(tT, Gene.symbol %in% control_genes)
    
    # Output the table of selected genes
    write.table(genes_subset, file = paste0(study_name, "_control_genes.tsv"), row.names = FALSE, sep = "\t")
    
    # Subset the table for the genes of interest, excluding NaN values and duplicates
    plot_genes_subset <- subset(tT, Gene.symbol %in% control_genes & is.finite(AveExpr) & is.finite(adj.P.Val) & is.finite(logFC) & !duplicated(Gene.symbol))
    
    
    # Visulize
    # volcano plot
    
    # Plot the scatter plot only if there are valid genes to plot
    if (nrow(plot_genes_subset) > 0) {
      pdf(paste0(study_name, "_scatter_plot.pdf"))
      plot(plot_genes_subset$logFC, -log10(plot_genes_subset$adj.P.Val),
           xlim = range(plot_genes_subset$logFC, na.rm = TRUE), ylim = range(-log10(plot_genes_subset$adj.P.Val), na.rm = TRUE),
           xlab = "Log Fold Change (LogFC)", ylab = "-log10(Adjusted p-value)",
           main = paste("Scatter Plot of Control Genes Expression",
                        "\n", study_name,
                        "\nCountry:", study_country,
                        " Tissue:", study_tissue,
                        " Samples:", study_num),
           pch = 16, cex = abs(plot_genes_subset$AveExpr) / max(plot_genes_subset$AveExpr),
           col = ifelse(plot_genes_subset$Gene.symbol == "TMED9", "red",
                        ifelse(plot_genes_subset$Gene.symbol == "GAPDH", "blue", "green")))
      
      # Add gene labels
      text(plot_genes_subset$logFC, -log10(plot_genes_subset$adj.P.Val),
           labels = plot_genes_subset$Gene.symbol, pos = 3)
    } else {
      print("No valid genes found for plotting.")
    }
    dev.off()
    },
    error=function(e) {
      print(paste("In study ",study_name ))
      print(paste("Acured Error:", e))
    }
  )
}

# Create a list to store study objects
studies <- list()

# Define a class to represent a study
Study <- setRefClass(
  "Study",
  fields = list(
    name = "character",
    country = "character",
    tissue = "character",
    n = "integer",
    group_string = "character",
    file_format = "character"
  ),
  methods = list(
    initialize = function(name, country, tissue, n, group_string, file_format) {
      name <<- name
      country <<- country
      tissue <<- tissue
      n <<- n
      group_string <<- group_string
      file_format <<- file_format
    }
  )
)

# Function to read the config file and populate the studies list
read_config_file <- function(filename) {
  studies <- list()
  
  con <- file(filename, "r")
  while (length(line <- readLines(con, n = 1)) > 0) {
    line <- trimws(line)
    
    if (nchar(line) > 0) {
      elements <- strsplit(line, ",")[[1]]
      
      if (length(elements) == 6) {
        name <- elements[1]
        country <- elements[2]
        tissue <- elements[3]
        n <- as.integer(elements[4])
        group_string <- elements[5]
        file_format <- elements[6]
        
        study <- Study$new(name, country, tissue, n, group_string, file_format)
        studies <- c(studies, study)
      } else {
        cat(paste("Ignoring invalid line:", line, "\n"))
      }
    }
  }
  close(con)
  
  return(studies)
}
working_dir <- "C:\\Users\\itaia\\OneDrive\\Project_TMED"
setwd(working_dir)
print("Reading config file")
studies <- read_config_file("Configs.txt")
num_studies <- length(studies)
counter <- 0
print(paste("Number of studies:", num_studies))

print("Performing analysis")
for (study in studies) {
  counter <- counter + 1
  print(paste("Running analysis for study:", study$name))
  print(paste("Country:", study$country))
  print(paste("Tissue:", study$tissue))
  print(paste("Number of Samples:", study$n))
  name <- study$name
  group_string <- study$group_string
  format <- study$file_format
  tissue <- study$tissue
  num <- study$n
  country <- study$country
  
  tryCatch(
    {
      # Call the function from the loaded script and pass arguments
      run_study_analysis(name, group_string, format, num, country, tissue)
      print(paste("Finished running analysis for", study$name))
      print(paste("Remaining studies:", num_studies - counter))
    },
    error = function(e) {
      print(paste("Error occurred for", study$name, "- Skipping to the next study."))
      print(paste("Error message:", conditionMessage(e)))
    }
  )
}
print("Finished running all analysis")
