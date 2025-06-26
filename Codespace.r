#Preparation & Heatmap plot
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)

# Create genotype vector
genotype <- c('smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe')

# Create condition vector
condition <- c('fibrosis','fibrosis','fibrosis','fibrosis','normal','normal','normal')

# Create data frame
smoc2_metadata <- data.frame(genotype,condition)

# Assign the row names of the data frame
rownames(smoc2_metadata) <- c('smoc2_fibrosis1','smoc2_fibrosis2','smoc2_fibrosis3','smoc2_fibrosis4','smoc2_normal1','smoc2_normal3','smoc2_normal4')

# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind=TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2) 

# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))

#Plot PCA
# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd_smoc2, intgroup="condition")

#Plot dispersi
# Plot dispersions
plotDispEsts(dds_smoc2)

# Extract the results of the differential expression analysis
smoc2_res <- results(dds_smoc2, 
                contrast = c('condition','fibrosis','normal'), 
                alpha = 0.05)

# Shrink the log2 fold change estimates to be more accurate
smoc2_res <- lfcShrink(dds_smoc2, 
                    contrast =  c('condition', 'fibrosis', 'normal'),
                    res = smoc2_res)

# Get an overview of the results                    
summary(smoc2_res)

# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_res)

# Subset the results to only return the significant genes with p-adjusted values less than 0.05
smoc2_res_sig <- subset(smoc2_res_all, padj < 0.05)

# Create MA plot
plotMA(smoc2_res)

# Generate logical column 
smoc2_res_all <- data.frame(smoc2_res) %>% mutate(threshold = padj < 0.05)
              
# Create the volcano plot
ggplot(smoc2_res_all) + 
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") + 
        theme(legend.position = "none", 
              plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25)))

# Log transform counts for QC
vsd_all <- vst(dds_all, blind = TRUE)

# Create heatmap of sample correlation values
vsd_all %>% 
        assay() %>% #Extract vst
        cor() %>% #Compute pairwise correlation
        pheatmap(annotation = select(all_metadata, c("genotype", "condition")))

# Create the PCA plot for PC1 and PC2 and color by condition       
plotPCA(vsd_all, intgroup = "condition")

# Create the PCA plot for PC1 and PC2 and color by genotype       
plotPCA(vsd_all, intgroup = "genotype")

# Select significant genese with padj < 0.05
smoc2_sig <- subset(res_all, padj < 0.05) %>%
  				data.frame() %>%
  				rownames_to_column(var = "geneID")

# Extract the top 6 genes with padj values
smoc2_sig %>%
	arrange(padj) %>%
	select(geneID, padj) %>%
	head()

# Create MA plot
plotMA(smoc2_res)

# Generate logical column 
smoc2_res_all <- data.frame(smoc2_res) %>% mutate(threshold = padj < 0.05)
              
# Create the volcano plot
ggplot(smoc2_res_all) + 
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") + 
        theme(legend.position = "none", 
              plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25)))