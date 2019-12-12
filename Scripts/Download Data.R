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