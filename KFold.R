run_Kfold_analysis <- function(study_name, study_group_string, study_format, study_num, study_country, study_tissue, threshold) {
  tryCatch(
    {
      library(BiocManager)
      library(GEOquery)
      library(limma)
      library(boot)
      
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
      gset <- gset[, sel]
      
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
      
      # Perform k-fold cross-validation
      k <- 10  # Specify the number of folds
      
      cv_fit <- cv.glm(data = exprs(gset), design, K = k)  # Perform k-fold cross-validation
      print(str(cv_fit))
      
      # Access the cross-validation results
      cv_error <- cv_fit$delta  # Get the cross-validation errors
      cv_indices <- cv_fit$which  # Get the indices for each fold
      
      # Loop over each fold and assess the performance
      for (fold in 1:k) {
        train_indices <- cv_indices[[fold]]$train  # Indices of samples in the training set
        test_indices <- cv_indices[[fold]]$test  # Indices of samples in the test set
        
        train_data <- exprs(gset)[, train_indices]  # Training data
        train_design <- design[train_indices, ]  # Training design matrix
        
        test_data <- exprs(gset)[, test_indices]  # Test data
        test_design <- design[test_indices, ]  # Test design matrix
        
        # Fit linear model on the training set
        train_fit <- lmFit(train_data, train_design)
        
        # Set up contrasts of interest and recalculate model coefficients on the training set
        train_cts <- paste(groups[1], groups[2], sep = "-")
        train_cont.matrix <- makeContrasts(contrasts = train_cts, levels = train_design)
        train_fit2 <- contrasts.fit(train_fit, train_cont.matrix)
        
        # Compute statistics and table of top significant genes on the training set
        train_fit2 <- eBayes(train_fit2, 0.01)
        train_tT <- topTable(train_fit2, adjust = "fdr", sort.by = "B", number = Inf)
        
        # Evaluate performance on the test set
        test_fit <- lmFit(test_data, test_design)
        test_fit2 <- contrasts.fit(test_fit, train_cont.matrix)  # Use the same contrasts from the training set
        test_tT <- topTable(test_fit2, adjust = "fdr", sort.by = "B", number = Inf)
        
        # ... Perform any desired analysis or evaluation on the training and test results ...
        # Initialize a set or list to store the selected genes that appear as significant in multiple folds
        selected_genes <- c()
        
        # Perform k-fold cross-validation
        for (fold in 1:k) {
          # ... previous code for fold-specific analysis ...
          
          # Compute the list of top significant genes (DEGs) for the current fold
          top_genes <- tT$Gene.symbol
          
          # Check if the top genes are already present in the selected_genes set/list
          for (gene in top_genes) {
            if (gene %in% selected_genes) {
              # Increment a counter or keep track of the fold index where the gene appears
              # You can use a counter or a vector to keep track of the fold index, depending on your requirements
              # For simplicity, let's assume we are using a counter variable `gene_count`
              gene_count[gene] <- gene_count[gene] + 1
            } else {
              # Add the gene to the selected_genes set/list and initialize its counter
              selected_genes <- c(selected_genes, gene)
              gene_count[gene] <- 1
            }
          }
        }
        
        # Determine the genes that consistently appeared as significant across multiple folds
        consistently_significant_genes <- names(gene_count[gene_count >= threshold])
        
        # Print the consistently significant genes
        print(consistently_significant_genes)
        write.table(consistently_significant_genes, file =paste0(study_name,"_topKfold.tsv"), row.names = FALSE, sep = "\t")
        
      }
      
    },
    error=function(e) {
      print(paste("In study ", study_name))
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
study <- studies[1]
name <- study$name
group_string <- study$group_string
format <- study$file_format
tissue <- study$tissue
num <- study$n
country <- study$country
run_Kfold_analysis(name, group_string, format, num, country, tissue)

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
      run_Kfold_analysis(name, group_string, format, num, country, tissue)
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