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
      
      # Subset to tables
      tT <- subset(tT, select = c("adj.P.Val", "P.Value", "AveExpr", "t", "B", "logFC", "Gene.symbol", "Gene.ID", "Gene.title"))
      # top 1000
      tT_top1000 <- head(tT, 1000)
      # Define Control Genes
      control_genes <- c("TMED9", "PDE4B", "NFKB1", "SRR", "GAPDH")
      # Subset the table for the genes of interest, including NaN values and duplicates
      genes_subset <- subset(tT, Gene.symbol %in% control_genes)
      # summarize test results as "up", "down" or "not expressed" by P val threshold and LFC
      dT <- decideTests(fit2, adjust.method="fdr", p.value=0.08, lfc=0)
      # Extract the indices of up regulated and down regulated genes from dT matrix
      upregulated_indices <- which(dT == 1)
      downregulated_indices <- which(dT == -1)
      # Extract the gene symbols corresponding to the indices
      gene_symbols <- tT$Gene.symbol
      # Obtain the top 5 up regulated genes
      top_up_genes <- gene_symbols[upregulated_indices[1:3]]
      # Obtain the top 5 down regulated genes
      top_down_genes <- gene_symbols[downregulated_indices[1:3]]
      # Define TMED family
      TMED_genes <- c("TMED1", "TMED2", "TMED3", "TMED4", "TMED5", "TMED6", "TMED7", "TMED9", "TMED10", "TMED10P1")
      # subset TMEDs
      TMED_genes_subset <- subset(tT, Gene.symbol %in% TMED_genes)
      # Create a vector of gene symbols for all genes in the analysis
      gene_symbols <- tT$Gene.symbol
      
      # Create an empty status table with gene symbols as the first column
      status_table <- data.frame(Gene.symbol = gene_symbols, Status = 0, stringsAsFactors = FALSE)
      
      # Set the status to 1 for TMED genes and control genes
      status_table$Status[gene_symbols %in% c(TMED_genes, control_genes)] <- 1
      
      
      # Output tables
      write.table(tT, file = paste0(study_name, "_output_table.tsv"), row.names = FALSE, sep = "\t")
      # Write the top 1000 genes to a separate file
      write.table(tT_top1000, file = paste0(study_name, "_top1000_genes.tsv"), row.names = FALSE, sep = "\t")
      # write the table of control genes
      write.table(genes_subset, file = paste0(study_name, "_control_genes.tsv"), row.names = FALSE, sep = "\t")
      # write the TMED family
      write.table(TMED_genes_subset, file = paste0(study_name, "_TMED_genes.tsv"), row.names = FALSE, sep = "\t")
      
      
      
      # Visualize
      # Open a PDF file for output
      pdf(paste0(study_name, "_plots.pdf"))
      
      # Create title
      general_title <- paste(
        "\n\nStudy:", study_name,
        "\n\nCountry:", study_country,
        "\n\nTissue:", study_tissue,
        "\n\nSamples:", study_num
      )
      top_up_text <- paste("\n\nSignificance cutoff by FDR < 0.08 (adj.P Value):\nTop up:", paste(top_up_genes, collapse = ", "))
      top_down_text <- paste("\nTop down:", paste(top_down_genes, collapse = ", "))
      full_title <- paste(general_title, top_up_text, top_down_text)
      
      # Set the plot margins
      par(mar = c(5, 4, 4, 2) + 0.1)
      
      # Plot the title as text without axes
      plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      text(x = 0.5, y = 0.5, full_title, cex = 1.2, font = 2, adj = 0.5)
      
      
      # volcano plot of full genes and list of top up and top down
      colnames(fit2) # list contrast names
      ct <- 1        # choose contrast of interest
      # Generate the title with top up and down genes:
      title_text <- paste("Volcano plot:", colnames(fit2)[ct])
      # Calculate the number of genes with non-zero status for the contrast of interest
      highlight_genes <- length(which(status_table$Status != 0))
      
      # Create a vector of '+' signs for all genes
      names_vector <- rep('+', nrow(fit2))
      
      # Generate the volcano plot using the modified parameters
      volcanoplot(fit2, coef = ct, main = title_text, pch = 20,
                  highlight = highlight_genes, names = names_vector)
      #add legend
      legend("topright", legend = "TMEDs and Control in color", xpd = TRUE)
      
      
      # volcano for TMED family and control
      unified_genes <- c(control_genes, TMED_genes)
      # Subset the table for the genes of interest, excluding NaN values and duplicates
      plot_genes_subset <- subset(tT, Gene.symbol %in% unified_genes & is.finite(AveExpr) & is.finite(P.Value) & is.finite(logFC) & !duplicated(Gene.symbol))
      # Define the shape based on the adj Pval condition
      shape <- ifelse(plot_genes_subset$adj.P.Val < 0.08, 24, 16)
      # Plot the data with different shape for significant, and color for control
      plot(plot_genes_subset$logFC, -log10(plot_genes_subset$P.Value),
           xlim = range(plot_genes_subset$logFC, na.rm = TRUE), ylim = range(-log10(plot_genes_subset$P.Value), na.rm = TRUE),
           xlab = "Log Fold Change (LogFC)", ylab = "-log10(p-value)",
           main = paste(title_text, "\nControl and TMED Genes"),
           pch = shape,
           col = ifelse(plot_genes_subset$Gene.symbol %in% c("PDE4B", "NFKB1", "SRR"), "red",
                        ifelse(plot_genes_subset$Gene.symbol == "GAPDH", "blue", "green")))
      
      
      # Add gene labels
      text(plot_genes_subset$logFC, -log10(plot_genes_subset$P.Value),
           labels = plot_genes_subset$Gene.symbol, pos = 3)
      # Add legend
      legend("topright", legend = "FDR Significant", pch = 24, xpd = TRUE)
      
      
      # MD plot (log fold change vs mean log expression)
      # Generate the MD plot using the modified parameters
      plotMD(fit2, column = ct, status = status_table$Status, legend = FALSE, pch = 20, cex = 1, main = "MD plot (log fold change vs mean log expression)")
      abline(h = 0)
      # add legend
      legend("topright", legend = "TMEDs and Control in color", xpd = TRUE)
      
      
      # Create MD  subset MD plot
      plot(plot_genes_subset$AveExpr, plot_genes_subset$logFC, 
           xlab = "Average Log Expression", ylab = "LogFC",
           main = "MD Plot\nControl and TMED Genes",
           pch = shape, col = ifelse(plot_genes_subset$Gene.symbol %in% c("PDE4B", "NFKB1", "SRR"), "red",
                                     ifelse(plot_genes_subset$Gene.symbol == "GAPDH", "blue", "green")))
      
      # Add gene symbols as labels to the points
      text(plot_genes_subset$AveExpr, plot_genes_subset$logFC,
           labels = plot_genes_subset$Gene.symbol, pos = 3, cex = 0.8, col = "black")
      
      # Add horizontal and vertical lines at logFC = 0
      abline(h = 0, lty = 2)
      abline(v = mean(plot_genes_subset$AveExpr), lty = 2)
      # Add legend
      legend("topright", legend = "FDR Significant", pch = 24, xpd = TRUE)
      
      
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
studies <- read_config_file("Configs_test.txt")
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
