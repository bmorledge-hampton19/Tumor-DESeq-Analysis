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
    
    #Truncate patient barcode to be consistent with clinical data.  
    #Then, remove duplicated barcodes and set barcode as the key.
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