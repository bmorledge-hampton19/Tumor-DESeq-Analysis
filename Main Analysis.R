###############################################################################################
#####                     Tumor RNA Seq Analysis Based on Size                            #####
###############################################################################################
# This script sets up an environment that is used by other scripts to analyse tumor RNA Seq
# data based on tumor size.  It then runs those scripts.  (See below for more info)

##### Set Up the Environment #####

#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("metaseqR")

#Include libraries
library(TCGAbiolinks)
library(DESeq2)
library(clusterProfiler)
library(SummarizedExperiment)
library(data.table)
library(org.Hs.eg.db)
library(biomaRt)
library(GO.db)
library(ggplot2)
library(edgeR)
library(metaseqR)

#Ensure Parallelization
#NOTE: THIS MAY NEED TO BE ADJUSTED BASED ON HOW MANY CPU CORES ARE AVAILABLE.
library(BiocParallel)
register(SnowParam(7))

#Get tumor types to analyze from csv
tumorTypes = fread("Sample Tumor Info.csv")
tumorIDs = as.character(tumorTypes$Abbreviation)



##### Analysis Scripts #####

### The following scripts run various steps of the tumor analysis and output their results into the "SavedData" folder.
### These scripts are meant to run sequentially, but they do not need to be run all at the same time, 
### because the results are saved at the end of each script.
### However, you must run the above code snippet to set up the proper environment every time you close RStudio.

#This script downloads data from TCGA based on the supplied tumor IDs.
#The results are outputted to the "Clinical" and "GeneExpression" folders.
source("Scripts/Download Data.R")

#This script formats the size and gene expression data, isolating T1, T2, and T3 tumors,
#removing non-unique patient barcodes, and combining size and expression data in one table.
#The results are outputted to the "SizeAndGeneExpression" folder.
source("Scripts/Format Data.R")

#The following three scripts run DESeq analysis, producing RNA expression comparisons between the tumor sizes.
#The results are outputted to the "DESeqObjects" and "DESeqRawData" folders.
source("Scripts/Setup for DESeq.R")
source("Scripts/Run DESeq.R")
source("Scripts/Extract DESeq Results.R")

#The following script runs Limma Voom analysis, producing RNA expression comparisons between the tumor sizes.
#The results are outputted to the "LimmaVoomObjects" and "LimmaVoomRawData" folders.
source("Scripts/Run Limma Voom.R")

#The following script places genes into relevant categories based on how they are associated with tumor growth.
#It also contains functions to create dotplots and heatmaps based on the results,
#but you will have to invoke them yourself as the script never actually utilizes them.
#(This is to save you from being bombarded with an absolutely ridiculous number of plots.)
#The results are outputted to the "DESeqRefinedData" and "LimmaVoomRefinedData" folders for each tumor type.
source("Scripts/Perform Biologically Relevant Analysis.R")

#The following script determines how common the genes are between tumor types for each of the categories created
#by the previous script.  It also produces barplots to display how well-represented each category was per tumor type.
#The results are outputted to the "DESeqRefinedData" and "LimmaVoomRefinedData" folders under the subfolder "ComparisonTables".
source("Scripts/Analyze Across Tumor Types.R")
