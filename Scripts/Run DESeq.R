##### Run DESeq #####
for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))) {
    
    #Load in DESeq Data
    load(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))
    
    #Run it!
    DESeqResults = DESeq(DESeqData,parallel = TRUE, BPPARAM = SnowParam(7))
    
    #Save the results
    save(DESeqResults, file = paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSR.rda"))
    
  }
  
}