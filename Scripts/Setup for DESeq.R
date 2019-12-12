##### Set up for DESeq #####

#Create the necessary directory for saving data.
dir.create("SavedData/DESeqObjects")

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))) {
    
    # Load in Data
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    # Split up the data
    sizeData = sizeAndExpressionData[,.(bcr_patient_barcode,stage_event_tnm_categories)]
    expressionData = sizeAndExpressionData[,!"stage_event_tnm_categories"]
    
    # Re-Transpose expression data to work with DeSeq
    expressionData = melt(expressionData, id.vars = "bcr_patient_barcode", 
                          variable.name = "Gene", value.name = "Counts")
    expressionData = dcast(expressionData, Gene ~ bcr_patient_barcode, value.var = "Counts")
    
    # Filter out lowly expressed Genes
    expressionData = expressionData[filterByExpr(expressionData[,-1], group = sizeData[,stage_event_tnm_categories])]
    
    # Convert to a matrix with row names
    expressionData = as.matrix(expressionData, rownames = TRUE)
    
    # Create the DESeqDataSet object, and save it.
    DESeqData = DESeqDataSetFromMatrix(countData = expressionData, 
                                       colData = sizeData, 
                                       design = ~ stage_event_tnm_categories)
    
    save(DESeqData, file = paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))
    
  }
  
}