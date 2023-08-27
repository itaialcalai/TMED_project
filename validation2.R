run_DEGs_analysis <- function(studies) {
  tryCatch(
    {
      library(BiocManager)
      library(GEOquery)
      library(randomForest)
      

      # Perform differential expression analysis for each dataset
      de_results <- list()
      
      # Download and store datasets
      for (i in 1:length(studies)) {
        study <- studies[[i+1]]
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
        
        print("111")
        
        # Select differentially expressed genes based on a threshold (e.g., adjusted p-value or fold change)
        de_results[[i]] <- decideTests(fit2)

      }
      print(de_results)
      # Combine the DEG results from each dataset into a single matrix
      de_matrix <- do.call(cbind, de_results)
      
      # Run RobustRankAggreg on the DEG matrix
      rra_result <- RRA(de_matrix)
      
      # Extract the integrated DEGs based on the result
      integrated_degs <- rownames(rra_result$aggregatedRank)
      print(integrated_degs)
     
      
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

run_DEGs_analysis(studies)
