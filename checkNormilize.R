# full analysis of the data
library(BiocManager)
library(GEOquery)
library(limma)

run_check_norm <- function(study_name, study_group_string, study_format, study_num, study_country, study_tissue) {
  gsms <- NULL
  sml <- NULL
  gset <- NULL
  ex <- NULL
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
      
      # normalize
      exprs(gset) <- normalizeBetweenArrays(exprs(gset))

      
      # assign samples to groups and set up the design matrix
      gs <- factor(sml)
      groups <- make.names(c("T2D", "normal"))
      levels(gs) <- groups
      gset$group <- gs
      design <- model.matrix(~group + 0, gset)
      colnames(design) <- levels(gs)
      
      gset <- gset[complete.cases(exprs(gset)), ]  # Skip missing values
      
      # fit <- lmFit(gset, design)  # Fit linear model
      
      # Set up contrasts of interest and recalculate model coefficients
      # cts <- paste(groups[1], groups[2], sep = "-")
      # cont.matrix <- makeContrasts(contrasts = cts, levels = design)
      # fit2 <- contrasts.fit(fit, cont.matrix)
      
      # Compute statistics and table of top significant genes
      # fit2 <- eBayes(fit2, 0.01)
      # tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
      
      # Plot before normalization
      pdf_file <- paste0(study_name, "_after_normalization_boxplot.pdf")
      pdf(pdf_file)
      # box-and-whisker plot
      ex <- exprs(gset)
      ord <- order(gs)  # order samples by group
      palette(c("#1B9E77", "#7570B3"))
      par(mar=c(7,4,2,1))
      title <- paste ("AFter Scale Normalization of Median Absolute Values", "/", annotation(gset), sep ="")
      boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
      legend("topleft", groups, fill=palette(), bty="n")
      # Close the PDF file
      
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
      run_check_norm(name, group_string, format, num, country, tissue)
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
