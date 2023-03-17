# Importing libraries

library(DESeq2)
library(tidyverse)
library(airway)

# Importing the data
counts_data <- read.csv("counts_data.csv")
head(counts_data)

# Importing metadata
colData <- read.csv("sample_info.csv")


# Match row names n colData to colum names in counts_data
# First, let's make sure that all the values are equal between the two
all(colnames(counts_data) %in% rownames(colData))

# ow we need to ensure they are in the same order
all(colnames(counts_data) == rownames(colData))


# Now we need to construct a DESeq2 object/input
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)
dds


# Pre-filtering by removing rows with low gene counts
# Let's apply a threshold of 10 read for each gene as a signifiacnt count
# we showed that more than 60% of genes are deleted bc of low reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]


# Setting a Factor level: to define untreated as the reference
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# Run DESeq
dds <- DESeq(dds)

# Save the results of DESeq
res <- results(dds)
res
summary(res)

# to change our p-value threshold
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# Generation of MA plot
plotMA(res)

