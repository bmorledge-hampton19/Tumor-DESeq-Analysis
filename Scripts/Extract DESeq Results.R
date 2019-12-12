##### Extract DESeq Results #####

# Create the necessary Directorys for saving data.
dir.create("SavedData/DESeqRawData")
dir.create("SavedData/DESeqRefinedData")

# A nice function for creating results from different comparisons
createResultsTable = function(DESeqResults, genes, size1, size2) {
  
  #Get the results for the indicated comparison
  results = results(DESeqResults, contrast = c("stage_event_tnm_categories",size1,size2), 
                    parallel = TRUE, BPPARAM = SnowParam(7))
  
  #Convert the results to a data.table and add the column for gene ID's.
  DESeqResultsTable = as.data.table(results)[,Gene := genes]
  setcolorder(DESeqResultsTable,c(7,1:6))
  
  return(DESeqResultsTable)
  
}

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSR.rda"))) {
    
    #Load in DESeq Results
    load(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSR.rda"))
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    #Create the comparison tables
    T1vsT2 = createResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T2")
    T1vsT3 = createResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T3")
    T2vsT3 = createResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T2", "T3")
    
    #Intersect T1vsT2 and T1vsT3 to get the genes that are consistently differentially expressed in large tumors.
    setkey(T1vsT2,Gene)
    setkey(T1vsT3,Gene)
    setkey(T2vsT3,Gene)
    
    T1vsT2IntersectT1vsT3 = T1vsT2[T1vsT3, nomatch = 0]
    
    #Format column names
    invisible(sapply(colnames(T1vsT2IntersectT1vsT3)[2:7], function(X)
      setnames(T1vsT2IntersectT1vsT3,old = X, new = paste0("T1vsT2.",X))))
    invisible(sapply(colnames(T1vsT2IntersectT1vsT3)[8:13], function(X)
      setnames(T1vsT2IntersectT1vsT3,old = X, new = paste0("T1vsT3",substr(X,2,1000)))))
    
    #Create a sub-directory to save the tables in.
    dir.create(paste0("SavedData/DESeqRawData/",tumorIDs[i]))
    dir.create(paste0("SavedData/DESeqRefinedData/",tumorIDs[i]))
    
    #Save the tables
    save(T1vsT2, file = paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T1vsT2.rda"))
    save(T1vsT3, file = paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T1vsT3.rda"))
    save(T2vsT3, file = paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T2vsT3.rda"))
    save(T1vsT2IntersectT1vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2_Intersect_T1vsT3.rda"))
    
  }
  
}