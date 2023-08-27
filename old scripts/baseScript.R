R.version
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")
library(GEOquery)
library(limma)

gset <- getGEO("GSE23343", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group membership for all samples
gsms <- "11111110000000000"
sml <- strsplit(gsms, split = "")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("T2D","normal"))
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
tT <- topTable(fit2, adjust = "fdr", sort.by = "adj.P.Val", number = Inf)

# check columns
colnames(tT)
# Subset and write the table
tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "AveExpr", "t", "B", "logFC", "Gene.symbol", "Gene.ID", "Gene.title"))
write.table(tT, file = "GSE23343_output_table.tsv", row.names = FALSE, sep = "\t")

# Write the top 1000 genes to a separate file
tT_top1000 <- head(tT_subset, 1000)
write.table(tT_top1000, file = "GSE23343_top1000_genes.tsv", row.names = FALSE, sep = "\t")

# Define the genes to be plotted
genes_to_plot <- c("TMED9", "PDE4B", "NFKB1", "SRR", "GAPDH")

# Subset the table for the genes of interest
genes_subset <- subset(tT, Gene.symbol %in% genes_to_plot)

# Plot the scatter plot
plot(log2(genes_subset$AveExpr), -log10(genes_subset$adj.P.Val), 
     xlim = range(log2(tT$AveExpr)), ylim = range(-log10(tT$adj.P.Val)),
     xlab = "Log Fold Change (LogFC)", ylab = "-log10(Adjusted p-value)",
     main = "Scatter Plot of Gene Expression",
     pch = 16, cex = genes_subset$AveExpr / max(genes_subset$AveExpr),
     col = ifelse(genes_subset$Gene.symbol == "TMED9", "red",
                  ifelse(genes_subset$Gene.symbol == "GAPDH", "blue", "green")))

# Add gene labels
text(log2(genes_subset$AveExpr), -log10(genes_subset$adj.P.Val), 
     labels = genes_subset$Gene.symbol, pos = 3)




