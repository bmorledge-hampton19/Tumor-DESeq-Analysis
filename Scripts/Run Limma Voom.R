##### Limma-Voom #####

#Create the necessary directorys for saving data.
dir.create("SavedData/LimmaVoomObjects")
dir.create("SavedData/LimmaVoomRawData")

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))) {
    
    # Load in Data
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    # Split up the data
    sizeData = sizeAndExpressionData[,.(bcr_patient_barcode,stage_event_tnm_categories)]
    expressionData = sizeAndExpressionData[,!"stage_event_tnm_categories"]
    
    # Re-Transpose expression data to traditional "counts" structure.
    expressionData = melt(expressionData, id.vars = "bcr_patient_barcode", 
                          variable.name = "Gene", value.name = "Counts")
    expressionData = dcast(expressionData, Gene ~ bcr_patient_barcode, value.var = "Counts")
    setkey(expressionData,Gene)
    
    # Derive normalization factors.
    dgeList = DGEList(as.matrix(expressionData, rownames = TRUE), 
                      group = as.vector(sizeData[,stage_event_tnm_categories]))
    normalizedData = calcNormFactors(dgeList)
    
    # Filter out low expressed genes
    normalizedAndTrimmedData = normalizedData[filterByExpr(normalizedData,),]
    
    # Apply the voom transformation
    tumorSizes = sizeData[,stage_event_tnm_categories]
    voomModel = model.matrix(~0 + tumorSizes)
    voomResults = voom(normalizedAndTrimmedData,voomModel,plot = FALSE)
    
    # TESTS
    #trimmedByOneFifth = apply(cpm(normalizedData), 1, function(x) sum(x>0.25)>dim(normalizedData)[2]/2)
    #trimmedByOneFifth = normalizedData[trimmedByOneFifth,]
    #voom(trimmedByOneFifth,voomModel,plot = TRUE)
    #
    #trimmedBy1 = apply(cpm(normalizedData), 1, function(x) sum(x>1)>dim(normalizedData)[2]/2)
    #trimmedBy1 = normalizedData[trimmedBy1,]
    #voom(trimmedBy1,voomModel,plot = TRUE)
    
    
    
    # Use Limma to fit the linear model (and save it)
    limmaFit = lmFit(voomResults,voomModel)
    save(limmaFit, file = paste0("SavedData/LimmaVoomObjects/TCGA-",tumorIDs[i],"-LFit.rda"))
    
    # Create the contrasts we want to analyze.
    T1vsT2Contr = makeContrasts(tumorSizesT1 - tumorSizesT2, levels = colnames(coef(limmaFit)))
    T1vsT3Contr = makeContrasts(tumorSizesT1 - tumorSizesT3, levels = colnames(coef(limmaFit)))
    T2vsT3Contr = makeContrasts(tumorSizesT2 - tumorSizesT3, levels = colnames(coef(limmaFit)))
    
    # Get the results (smoothed using empirical Bayes)
    T1vsT2 = topTable((eBayes(contrasts.fit(limmaFit, T1vsT2Contr))), 
                      sort.by = "none", number = Inf)
    T1vsT3 = topTable((eBayes(contrasts.fit(limmaFit, T1vsT3Contr))), 
                      sort.by = "none", number = Inf)
    T2vsT3 = topTable((eBayes(contrasts.fit(limmaFit, T2vsT3Contr))), 
                      sort.by = "none", number = Inf)
    
    # Convert the results to neat and tidy data.tables.
    T1vsT2 = as.data.table(T1vsT2, keep.rownames = "Gene")
    T1vsT3 = as.data.table(T1vsT3, keep.rownames = "Gene")
    T2vsT3 = as.data.table(T2vsT3, keep.rownames = "Gene")
    
    # Set keys
    setkey(T1vsT2,Gene)
    setkey(T1vsT3,Gene)
    setkey(T2vsT3,Gene)
    
    #Standardize column names with DESeq Data.
    setnames(T1vsT2,c("logFC","P.Value","adj.P.Val"),c("log2FoldChange","pvalue","padj"))
    setnames(T1vsT3,c("logFC","P.Value","adj.P.Val"),c("log2FoldChange","pvalue","padj"))
    setnames(T2vsT3,c("logFC","P.Value","adj.P.Val"),c("log2FoldChange","pvalue","padj"))
    
    #Create a sub-directory to save the tables in.
    dir.create(paste0("SavedData/LimmaVoomRawData/",tumorIDs[i]))
    
    #Save the tables
    save(T1vsT2, file = paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T1vsT2.rda"))
    save(T1vsT3, file = paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T1vsT3.rda"))
    save(T2vsT3, file = paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T2vsT3.rda"))
    
  }
  
}