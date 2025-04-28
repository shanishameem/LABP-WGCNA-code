# --------------------------------------------------
#
# Supplementary code for WGCNA
#
# --------------------------------------------------


### Bioconductor and CRAN libraries used

library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggplot2)
library(tximport)
library(ggrepel)
library(EnsDb.Hsapiens.v86)
library(factoextra)
library(dplyr)
library(WGCNA)
library(CorLevelPlot)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

# --------------------------------------------------
# Set working directory
# --------------------------------------------------

getwd()
setwd()
dir <- "~/basic_quant"

# List all directories containing data

samples <- read.table(file.path(dir, "basic_samples.txt"), header = T)

# Obtain a vector for all path names 

files <- file.path(dir, samples$sample_id, "quant.sf")



# --------------------------------------------------
# Building transcript / gene table with Ensembl human index
# --------------------------------------------------

edb <- EnsDb.Hsapiens.v86

k <- keys(edb, keytype = "TXNAME")
tx2gene <- transcripts(edb, columns = c("tx_id", "gene_id"))
tx2gene <- data.frame(
  TXNAME = mcols(tx2gene)$tx_id,  
  GENEID = mcols(tx2gene)$gene_id  
)
head(tx2gene)


txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE) # Note that this ignores the transcript, and uses lengthscaledTPM, i.e. it accounts for bias in the length of transcripts

# Save counts as data

basic_rnaseq_data <- txi$counts %>%
  round() %>%
  data.frame()


## Annotations

annot <- as.data.frame(genes(edb))
annot <- subset(annot, select = c("gene_id", "gene_name"))


# --------------------------------------------------
# Create meta data
# --------------------------------------------------

sampleTable <- data.frame(timepoint = substr(colnames(basic_rnaseq_data), 10,12))

meta <- data.frame(sampleTable, row.names = colnames(txi$counts))

patient <- as.data.frame(substr(row.names(meta),1,8))
names(patient)[names(patient) =="substr(row.names(meta), 1, 8)"] <- "basic_id"
patient_data_link <- as.data.frame(cbind(patient$basic_id, meta$timepoint))
names(patient_data_link)[names(patient_data_link) =="V1"] <- "basic_id"
names(patient_data_link)[names(patient_data_link) =="V2"] <- "timepoint"
patient_data_link$timepoint <- as.factor(patient_data_link$timepoint)

# --------------------------------------------------
# Apply DESeq2 normalization
# --------------------------------------------------

# Match metadata and counts data 

all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

# Create DESeq2Dataset object

dds_wgcna <- DESeqDataSetFromTximport(txi, colData = patient_data_link, design = ~1) # design not specified

# Keep only rows with more than 10 reads

keep <- rowSums(counts(dds_wgcna)) >= 10
dds_wgcna <- ddsTxi[keep,]

# Apply normalization 

dds_wgcna_norm <- vst(dds_wgcna)
wgcna_norm <- assay(dds_wgcna_norm) %>%
  t()

wgcna_norm2 <- assay(dds_wgcna_norm)

# Check outlier genes

gsg <-goodSamplesGenes(wgcna_norm)
summary(gsg)
gsg$allOK

# Cluster of samples

clust_samples <- hclust(d=dist(wgcna_norm), method = "average")
plot(clust_samples)

# Plot RNA seq counts

mmdata <- assay(dds_wgcna_norm)
mmdata <- as.data.frame(mmdata)
mcol_sel = names(mmdata)  
mmdata <-pivot_longer(mmdata,
                       col = all_of(mcol_sel))

p <- mmdata %>%
  ggplot(., aes(x = name, y = value)) +             
  geom_violin() +                                   
  geom_point(alpha = 0.2) +                         
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Samples", y = "RNA Seq Counts") 

p

# --------------------------------------------------
# Determine soft-threshold power
# --------------------------------------------------

# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function

sft = pickSoftThreshold(
  wgcna_norm,            
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

# Plot Scale independence and Mean connectivity

par(mfrow = c(1,2)); 
cex1 = 0.8; 

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red", pos = 3
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# --------------------------------------------------
# WGCNA construction
# --------------------------------------------------

set.seed(12345)
picked_power = 5
temp_cor <- cor       
cor <- WGCNA::cor        

netwk <- blockwiseModules(wgcna_norm,                
                          
                          # == Adjacency Function ==
                          power = picked_power,               
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 30000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file 
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = F,
                          verbose = 3
)
cor <- temp_cor  

# --------------------------------------------------
# Analysis of modules
# --------------------------------------------------

# Save mergred and unmerged modules

mergedColors = netwk$colors
unmergedColors = netwk$unmergedColors

# Plot the dendrogram and the module colors underneath

plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)


# Get Module Eigengenes per cluster

MEs0 <- moduleEigengenes(wgcna_norm, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# --------------------------------------------------
# Clustering of modules
# --------------------------------------------------

ME.dissimilarity = 1-cor(MEs0, use="complete")

METree = hclust(as.dist(ME.dissimilarity), method = "average") 
par(mar = c(0,4,2,0)) 
par(cex = 0.6);
plot(METree)
abline(h=.25, col = "red") 

# --------------------------------------------------
# Relate modules to external traits 
# --------------------------------------------------

# Save module dataframe

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = (netwk$colors)
)


write_delim(module_df,
            file = "gene_modules_final_sft5.txt",
            delim = "\t")

# Get Module Eigengenes per cluster

MEs0 <- moduleEigengenes(wgcna_norm, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other

MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Extract external traits of interest

WGCNA_traits <- read_csv("~/wgcnatraits.csv")
alltraits <- WGCNA_traits[,c(1, 3)]
alltraits <- column_to_rownames(alltraits, "patient")

# Correlation of modules to traits

MEs0_o <- MEs0[order(rownames(MEs0)),]
module.trait.correlation <- cor(MEs0_o, alltraits, use = 'p')
module.trait.correlation.pvals <- corPvalueStudent(module.trait.correlation, nrow(input_mat))

# Will display correlations and their p-values

textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.correlation.pvals, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)

# Display the correlation values within a heatmap plot

par(mar = c(4, 9, 4, 1)) 

labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(alltraits),
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = T,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))





# --------------------------------------------------
# Functional enrichment analysis for modules
# --------------------------------------------------

## BLUE MODULE

blue_mod_genes <- module_df %>%
  filter(colors == "blue")

# Run GO analysis

blue_mod_genes$Entrez = mapIds(org.Hs.eg.db, keys = blue_mod_genes$gene_id, 
                               column = "ENTREZID", keytype = "ENSEMBL", 
                               multiVals = "first") #Ensembl gene IDs to Entrez gene IDs

blue_mod_genes <- blue_mod_genes[!is.na(blue_mod_genes$Entrez),]

blue_mod_genes$Symbol <- mapIds(org.Hs.eg.db, keys = blue_mod_genes$Entrez, 
                                column = "SYMBOL", keytype = "ENTREZID", 
                                multiVals = "first") #Entrez gene IDs to gene symbols

blue_mod_go <- enrichGO(gene = blue_mod_genes$Entrez,
                        keyType = "ENTREZID",
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

dotplot(blue_mod_go, title="Blue module")


## MAGENTA MODULE

magenta_mod_genes <- module_df %>%
  filter(colors == "magenta")

# Run GO analysis

magenta_mod_genes$Entrez = mapIds(org.Hs.eg.db, keys = magenta_mod_genes$gene_id, 
                                  column = "ENTREZID", keytype = "ENSEMBL", 
                                  multiVals = "first") #Ensembl gene IDs to Entrez gene IDs

magenta_mod_genes <- magenta_mod_genes[!is.na(magenta_mod_genes$Entrez),]

magenta_mod_go <- enrichGO(gene = magenta_mod_genes$Entrez,
                           keyType = "ENTREZID",
                           OrgDb = org.Hs.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05,
                           readable = TRUE)

dotplot(magenta_mod_go, title="Magenta module")

## BROWN MODULE

brown_mod_genes <- module_df %>%
  filter(colors == "brown")

# Run GO analysis

brown_mod_genes$Entrez = mapIds(org.Hs.eg.db, keys = brown_mod_genes$gene_id, 
                                column = "ENTREZID", keytype = "ENSEMBL", 
                                multiVals = "first") #Ensembl gene IDs to Entrez gene IDs

brown_mod_genes <- brown_mod_genes[!is.na(brown_mod_genes$Entrez),]

brown_mod_go <- enrichGO(gene = brown_mod_genes$Entrez,
                         #universe = allgenes,
                         keyType = "ENTREZID",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

dotplot(brown_mod_go, title="Brown module")



# --------------------------------------------------
# Module membership & Gene Significance calculation + Hub gene identification
# --------------------------------------------------

# Define TP2 variable

vtp2 = as.data.frame(alltraits$TP2)
names(vtp2) = "TP2"

modNames = (names(MEs0)) 

# Calculate the module membership and the associated p-values

geneModuleMembership = as.data.frame(cor(wgcna_norm, MEs0, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(wgcna_norm)))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# Calculate the gene significance and associated p-values

# Order in ascending to match vtp2 

wgcna_norm_ordered <- wgcna_norm[order(rownames(wgcna_norm)),]


geneTraitSignificance = as.data.frame(cor(wgcna_norm_ordered, vtp2, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(wgcna_norm)))
names(geneTraitSignificance) = paste("GS.", names(vtp2), sep="")
names(GSPvalue) = paste("p.GS.", names(vtp2), sep="")
head(GSPvalue)



# --------------------------------------------------
#  Plot Module membership vs Gene significance
# --------------------------------------------------


module = "MEmagenta"
column = match(module, modNames)
moduleGenes = module_df$colors =="magenta"
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TP2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "magenta")

module = "MEblue"
column = match(module, modNames)
moduleGenes = module_df$colors =="blue"
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TP2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")

module = "MEbrown"
column = match(module, modNames)
moduleGenes = module_df$colors =="brown"
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TP2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "brown")


