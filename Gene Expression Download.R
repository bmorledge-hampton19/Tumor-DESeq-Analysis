#Need to align gene expression to tumor size using patient barcode.
#BUT sometimes multiple samples are taken from one patient.
#Each patient is only associated with at most one tumor size, so those with multiple samples need to be omitted.
#Also, we only care about tumor sizes T1, T2, and T3.
#(Thanks Pete!)
#Then, run DESeq analysis on the data.


#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")

#Include libraries
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(data.table)

#Ensure arallelization
library(BiocParallel)
register(SnowParam(7))

#Get tumor types to analyze from csv
tumorTypes = fread("Sample Tumor Info.csv")
tumorIDs = as.character(tumorTypes$Abbreviation)



##### Download ALL the data!! #####

#Create the necessary directories for saving data.
dir.create("SavedData/GeneExpression")
dir.create("SavedData/Clinical")

for (i in 1:length(tumorIDs)) {
  
  ##### Get Gene Expression data.
  
  #Construct the query and download the raw data.
  query <- GDCquery(project = paste0("TCGA-",tumorIDs[i]),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    experimental.strategy = "RNA-Seq")
  GDCdownload(query)
  data <- GDCprepare(query)
  
  # We want only primary tumor data, and we want it in a nice data.table.
  expressionData = as.data.table(assay(data[,data$shortLetterCode == "TP"]),keep.rownames = TRUE)
  
  # Save the results.
  save(expressionData, file = paste0("SavedData/GeneExpression/TCGA-",tumorIDs[i],"-GE.rda"))

  ##### Get Clinical data.
  
  query <- GDCquery(project = paste0("TCGA-",tumorIDs[i]),
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
  clinical <- as.data.table(GDCprepare_clinic(query, clinical.info = "patient"))
  save(clinical, file = paste0("SavedData/Clinical/TCGA-",tumorIDs[i],"-C.rda"))
    
}



##### Format Data #####

#Create the necessary directory for saving data.
dir.create("SavedData/SizeAndGeneExpression")

for (i in 1:length(tumorIDs)) {
  
  ##### Load in Clinical and Gene Expression data
  
  load(paste0("SavedData/GeneExpression/TCGA-",tumorIDs[i],"-GE.rda"))
  load(paste0("SavedData/Clinical/TCGA-",tumorIDs[i],"-C.rda"))

  ##### Format Clinical Data

  #Ensure that data actually contains a size column
  if ("stage_event_tnm_categories" %in% colnames(clinical)) {
       
    #Separate tumor size and patient barcode from other columns and truncate tumor size to first two characters.
    clinical = clinical[,.(bcr_patient_barcode,stage_event_tnm_categories)]
    clinical[,stage_event_tnm_categories := substr(stage_event_tnm_categories,1,2)]
    
    #Remove any patient barcodes which appear more than once and set the remaining barcodes as key.
    clinical = clinical[isUnique(clinical[,bcr_patient_barcode])]
    setkey(clinical,bcr_patient_barcode)
    
    #Get only the T1, T2, and T3 size tumors.
    clinical = clinical[grepl("T1|T2|T3", clinical[,stage_event_tnm_categories])]
    
  } else {
    clinical = NULL
  }
  
  # Ensure that we have valid size data before proceeding
  if (!is.null(clinical) && clinical[,.N] > 0) {
  
    ##### Format Expression Data
    
    #Rename column
    colnames(expressionData)[colnames(expressionData)=="rn"] = "Gene"
    
    #Transpose expression Data so that patient barcode can be used as a key.
    expressionData = melt(expressionData, id.vars = "Gene", variable.name = "Patient_Barcode", value.name = "Counts")
    expressionData = dcast(expressionData, Patient_Barcode ~ Gene, value.var = "Counts")
    
    #Truncate patient barcode to be consistent with clinical data, then remove duplicated barcodes and set barcode as the key.
    expressionData[,Patient_Barcode := substr(Patient_Barcode,1,12)]
    expressionData = expressionData[isUnique(expressionData[,Patient_Barcode])]
    setkey(expressionData,Patient_Barcode)
    
  
  
    ##### Join the two data sets and save the new data set

    #Join the two data sets
    sizeAndExpressionData = clinical[expressionData, nomatch=0]

    #Save the new data set.
    save(sizeAndExpressionData, file = paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))

  }
  
}



##### Set up for DESeq #####

#Create the necessary directory for saving data.
dir.create("SavedData/DESeqRawData")

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))) {
  
    # Load in Data
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    # Split up the data
    sizeData = sizeAndExpressionData[,.(bcr_patient_barcode,stage_event_tnm_categories)]
    expressionData = sizeAndExpressionData[,!"stage_event_tnm_categories"]
    
    # Re-Transpose expression data to work with DeSeq
    expressionData = melt(expressionData, id.vars = "bcr_patient_barcode", variable.name = "Gene", value.name = "Counts")
    expressionData = dcast(expressionData, Gene ~ bcr_patient_barcode, value.var = "Counts")
    
    # Create the DESeqDataSet object, and save it.
    DESeqData = DESeqDataSetFromMatrix(countData = expressionData, 
                                 colData = sizeData, 
                                 design = ~ stage_event_tnm_categories,
                                 tidy = TRUE)
    
    save(DESeqData, file = paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSD.rda"))
    
  }
  
}

##### Sanity Check #####

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSD.rda"))) {
    
    #Load in DESeq Data and clinical data
    load(paste0("SavedData/DESeqData/TCGA-",tumorIDs[i],"-DSD.rda"))
    load(paste0("SavedData/Clinical/TCGA-",tumorIDs[i],"-C.rda"))
    
    #Get just the patient barcodes from each set.
    colDataBarcodes = colData(DESeqData)$bcr_patient_barcode
    originalBarcodes = clinical[,bcr_patient_barcode]
    
    #Ensure that all patient barcodes in the DESeqData object are present in the original clinical data.
    if (length(intersect(colDataBarcodes,originalBarcodes)) == length(colDataBarcodes)) {
      print(paste0(tumorIDs[i]," tumor type successfully validated."))
    } else {
      print(paste0("Error in ",tumorIDs[i]," tumor type."))
    }
    
  }
  
}

##### Run DESeq #####
for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqData/TCGA-",tumorIDs[i],"-DSD.rda"))) {
  
    #Load in DESeq Data
    load(paste0("SavedData/DESeqData/TCGA-",tumorIDs[i],"-DSD.rda"))
    
    #Run it!
    DESeqResults = DESeq(DESeqData,parallel = TRUE, BPPARAM = SnowParam(7))
    
    #Save the results
    save(DESeqResults, file = paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSR.rda"))
  
  }
  
}

# Analyze Results
# padj cutoff: < 0.05
# increase cutoff: > log(1.5) = ~0.58
# Intersect T1 vs T2 and T1 vs T3
# What genes are consistently upregulated across tumor types?
# How can we group tumor types biologically?  (Tissue type?  Mesenchyme vs epithelial?)

# Create the necessary Directory for saving data.
dir.create("SavedData/DESeqRefinedData")

# A nice function for creating results from different comparisons
createRefinedResultsTable = function(DESeqResults, genes, size1, size2) {
  
  #Get the results for the indicated comparison
  results = results(DESeqResults, contrast = c("stage_event_tnm_categories",size1,size2))
  
  #Convert the results to a data.table and add the column for gene ID's.
  DESeqResultsTable = as.data.table(results)[,Gene := genes]
  setcolorder(DESeqResultsTable,c(7,1:6))
  
  #Remove insignificant changes and small changes.
  DESeqResultsTable = DESeqResultsTable[padj < 0.05 & abs(log2FoldChange) > 0.58]
  
  return(DESeqResultsTable)
  
}

for (i in 1:length(tumorIds)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSR.rda"))) {
    
    #Load in DESeq Results
    load(paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSR.rda"))
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    #Create the comparison tables
    T1vsT2 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T2")
    T1vsT3 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T3")
    T2vsT3 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T2", "T3")
    
    #Create a sub-directory to save the tables in.
    dir.create(paste0("SavedData/DESeqRefinedData/",tumorIDs[i]))
    
    #Save the tables
    save(T1vsT2, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2.rda"))
    save(T1vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT3.rda"))
    save(T2vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T2vsT3.rda"))
    
  }
  
}

# DAVID? Cluster Profiler?  GG Plot?  SEUART?