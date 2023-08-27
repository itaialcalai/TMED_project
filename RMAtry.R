# install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
biocLite("affy")  # Methods for Affymetrix Oligonucleotide Arrays
BiocManager::install("oligo")
biocLite("hgu133a.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
biocLite("hgu133acdf")

library(GEOquery)
library(affy)
library(hgu133a.db)
library(hgu133acdf)
library(data.table)
library(dplyr)
library(tidyverse)
library(oligo)
library(stringr)

run_study_analysis <- function(study_name, study_group_string, study_format, study_num, study_country, study_tissue) {
  tryCatch(
    {
      # get supplimentary files
      getGEOSuppFiles(study_name)

      untar(paste0(study_name,"_RAW.tar"), exdir = "data")
      cels = list.files("data/", pattern = "CEL|cel")
      # sometiles, it is 'CEL', you need to check it first
      sapply(paste("data", cels, sep = "/"), gunzip)
      
      cels = list.files("data/", pattern = "CEL|cel")
      # sometiles, it is 'CEL', you need to check it first
      
      # Set working directory for normalization
      setwd(paste0("./GEOdata/", study_name))
      
      ### for Affy package
      raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
      # perform RMA normalization (log2)
      data.rma.norm = affy::rma(raw.data)
      # Get the expression estimates for each array
      rma = exprs(data.rma.norm)
      
      ### for oligo package
      rawData <- read.celfiles(cels)
      # perform RMA normalization (log2)
      data.rma.norm = rma(rawData)
      # Get the expression estimates for each array
      rma = exprs(data.rma.norm)
      
      
      #---------------------------------- Annotation -------------------------------------#
      tt = cbind(row.names(rma), rma)
      colnames(tt) = c("ProbID", sub(".cel", "", colnames(rma), ignore.case = TRUE))
      
      #### to merge with gene symbols
      # if featureData has no column called Gene Symbol:
      gse <- getGEO(study_name,GSEMatrix=TRUE)
      featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
      #featureData <-  setDT(featureData)[,c("V1","Gene Symbol","V2"):= sapply(str_split(featureData$gene_assignment, " // ",  n = 3), `[`, 2)]
      featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
      #featureData$ID <-  as.character(featureData$ID)
      GSE17714_RMA_data <- merge(tt,featureData, by.x = "ProbID", by.y = "ID")
      
      #colnames(GSE107333_series_matrix_expr_data) <- sub("_exp.*", "", colnames(GSE107333_series_matrix_expr_data))
      
      
      #---------------------------------- Expression matrix Manipulation -------------------------------------#
      # converting gene symbols to row names and keeping rows with max mean in each row
      # removing ID column
      GSE_RMA_data <- select(GSE_RMA_data,-c("ProbID"))
      
      
      copy_GSE_RMA_data <-  GSE_RMA_data
      #GSE16476_series_matrix_expr_data <-  copy_GSE16476_series_matrix_expr_data
      
      cols <- names(GSE_RMA_data)[1:22]
      setDT(GSE_RMA_data)[, (cols) := lapply(.SD, as.character), .SDcols = cols]
      setDT(GSE_RMA_data)[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
      
      # Since there are duplicated gene symbols -
      # calculating means for every row grouping by the gene; and selecting rows with max mean
      GSE_RMA_data <- setDT(GSE4_RMA_data)[, .SD[which.max(rowMeans(.SD))], by=`Gene Symbol`]
      
      # Removing if there are any NA's in Gene Symbols
      GSE4_RMA_data <- GSE_RMA_data[!(is.na(GSE17714_RMA_data$`Gene Symbol`)),]
      
      # converting column to rownames
      GSE_RMA_data <- GSE_RMA_data %>% column_to_rownames(var="Gene Symbol")
      
      # Saving data
      save(GSE_RMA_data,file = paste0(study_name, "_RMA_data.RData"))
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
