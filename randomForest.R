install.packages("randomForest")

run_RR_analysis <- function(studies) {
  tryCatch(
    {
      library(GEOquery)
      library(randomForest)
      
      # Create lists to store expression data and labels from each dataset
      expression_list <- list()
      labels_list <- list()
      
      # Download and store datasets
      for (study in studies) {
        gset <- getGEO(study$name, GSEMatrix = TRUE, AnnotGPL = TRUE)
        if (length(gset) > 1) idx <- grep(study$file_format, attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # Make proper column names to match toptable
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # Group membership for all samples
        gsms <- study$group_string
        sml <- strsplit(gsms, split = "")[[1]]
        
        # Filter out excluded samples (marked as "X")
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[, sel]
        
        # Assign labels to the samples based on 'sml'
        labels <- ifelse(sml == "1", "healthy", ifelse(sml == "0", "sick", NA))
        
        expression_list[[study$name]] <- exprs(gset)
        labels_list[[study$name]] <- labels
      }
      print("1\n")
      for (i in 1:length(expression_list)) {
        print(ncol(expression_list[[i]]))
      }
      
      
      # Combine expression data and labels from all datasets into a single data frame
      all_expression_data <- do.call(rbind, expression_list)
      all_labels <- unlist(labels_list)
      
      print("2\n")
      # Split the data into training and testing sets
      set.seed(123)
      train_indices <- sample(1:nrow(all_expression_data), 0.7 * nrow(all_expression_data))
      train_data <- all_expression_data[train_indices, ]
      train_labels <- all_labels[train_indices]
      test_data <- all_expression_data[-train_indices, ]
      test_labels <- all_labels[-train_indices]
      
      print("3\n")
      
      # Train the random forest model
      rf_model <- randomForest(train_data, train_labels, ntree = 100)
      
      # Predict labels for the test data
      predictions <- predict(rf_model, test_data)
      
      # Evaluate model performance
      accuracy <- sum(predictions == test_labels) / length(test_labels)
      print(paste("Accuracy:", accuracy))
      
      print("4\n")
      
      # Determine importance of a specific gene
      # Replace 'gene_name' with the actual name of the gene of interest
      gene_importance <- importance(rf_model)[, "gene_name"]
      print(gene_importance)
      
      # Perform other analysis as needed
      
    },
    error = function(e) {
      # Error handling code
      print("An error occurred.")
      print(e)
    }
  )
}
# Create a list to store study objects
studies <- list()
working_dir <- "C:\\Users\\itaia\\OneDrive\\Project_TMED"
setwd(working_dir)
print("Reading config file")
studies <- read_config_file("Configs_test.txt")
num_studies <- length(studies)
counter <- 0
print(paste("Number of studies:", num_studies))

print("Performing analysis")

run_RR_analysis(studies)
