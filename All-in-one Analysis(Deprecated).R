###############################################################################################
#####                     Tumor RNA Seq Analysis Based on Size                            #####
###############################################################################################

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

##### Sanity Check #####

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))) {
    
    #Load in DESeq Data and clinical data
    load(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))
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
  if (file.exists(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))) {
  
    #Load in DESeq Data
    load(paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSD.rda"))
    
    #Run it!
    DESeqResults = DESeq(DESeqData,parallel = TRUE, BPPARAM = SnowParam(7))
    
    #Save the results
    save(DESeqResults, file = paste0("SavedData/DESeqObjects/TCGA-",tumorIDs[i],"-DSR.rda"))
  
  }
  
}



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



##### Biologically relevant analysis #####

generateHeatMap = function(genes, geneOntologyCategory, foldChange) {
  
  geneOntologyResults = enrichGO(genes, 'org.Hs.eg.db', 
                                 keyType = "ENSEMBL", ont=geneOntologyCategory, pvalueCutoff = 0.1)
  
  if(missing(foldChange)) {
    heatplot(setReadable(geneOntologyResults, "org.Hs.eg.db", "ENSEMBL"))
  } else {
    genesWithFoldChanges = setNames(foldChange,genes)
    heatplot(setReadable(geneOntologyResults, "org.Hs.eg.db", "ENSEMBL"), foldChange = genesWithFoldChanges)
  }
  
}

generateDotPlot = function(genes, geneOntologyCategory, title) {
  geneOntologyResults = enrichGO(genes, 'org.Hs.eg.db', pvalueCutoff = 0.5, qvalueCutoff = 0.5, 
                                 keyType = "ENSEMBL", ont=geneOntologyCategory)
  if(missing(title)) {
    enrichplot::dotplot(geneOntologyResults)
  } else {
    enrichplot::dotplot(geneOntologyResults) + ggtitle(title)
  }
  
}

generateBiologicallyRelevantResults = function(T1vsT2, T1vsT3, T2vsT3, tumorType, saveFolder) {
  
  print(paste0("Generating results for ",tumorType," in ",saveFolder,"."))
  
  #Find genes that increase in expression at each increase in tumor size.
  consistentGrowth = merge(T1vsT2[log2FoldChange<0],T2vsT3[log2FoldChange<0],
                           suffixes = c("T1vsT2","T2vsT3"))
  consistentGrowth = filterOnCombinedPAdj(consistentGrowth,consistentGrowth[,.(pvalueT1vsT2,pvalueT2vsT3)])
  
  #Find genes that increase in expression from T1 to T3 but were not picked up in "consistentGrowth"
  lateGrowth = T1vsT3[log2FoldChange<0 & padj<0.05]
  lateGrowth = lateGrowth[!(Gene %in% consistentGrowth$Gene)]
  
  #These genes are expressed highly in T2 tumors, but expression is brought back down in T3 tumors.
  earlyGrowth = merge(T1vsT2[log2FoldChange<0],T2vsT3[log2FoldChange>0],
                      suffixes = c("T1vsT2","T2vsT3"))
  earlyGrowth = filterOnCombinedPAdj(earlyGrowth,earlyGrowth[,.(pvalueT1vsT2,pvalueT2vsT3)])
  
  #Give some info...
  print(paste0("Found ",dim(consistentGrowth)[1]," genes in consistentGrowth."))
  print(paste0("Found ",dim(lateGrowth)[1]," genes in lateGrowth."))
  print(paste0("Found ",dim(earlyGrowth)[1]," genes in earlyGrowth."))
  
  #Save "binned" results.
  dir.create(paste0("SavedData/", saveFolder, "/", tumorType))
  save(consistentGrowth, file = paste0("SavedData/", saveFolder, "/",tumorType,"/consistentGrowth.rda"))
  save(lateGrowth, file = paste0("SavedData/", saveFolder, "/",tumorType,"/lateGrowth.rda"))
  save(earlyGrowth, file = paste0("SavedData/", saveFolder, "/",tumorType,"/earlyGrowth.rda"))
  
  #Generate and save gene ontology info.
  dir.create(paste0("SavedData/", saveFolder, "/", tumorType, "/GO"))
  
  if (dim(consistentGrowth)[1] > 0) {
    geneOntologyConsistentGrowth = enrichGO(consistentGrowth$Gene, 'org.Hs.eg.db', 
                                            pvalueCutoff = 0.5, qvalueCutoff = 0.5, 
                                            keyType = "ENSEMBL", ont="BP")
    save(geneOntologyConsistentGrowth,
         file = paste0("SavedData/", saveFolder, "/",tumorType,"/GO/geneOntologyConsistentGrowth.rda"))
  }
  
  if (dim(lateGrowth)[1] > 0) {
    geneOntologyLateGrowth = enrichGO(lateGrowth$Gene, 'org.Hs.eg.db', 
                                            pvalueCutoff = 0.5, qvalueCutoff = 0.5, 
                                            keyType = "ENSEMBL", ont="BP")
    save(geneOntologyLateGrowth,
         file = paste0("SavedData/", saveFolder, "/",tumorType,"/GO/geneOntologyLateGrowth.rda"))
  }
  
  if (dim(earlyGrowth)[1] > 0) {
    geneOntologyEarlyGrowth = enrichGO(earlyGrowth$Gene, 'org.Hs.eg.db', 
                                            pvalueCutoff = 0.5, qvalueCutoff = 0.5, 
                                            keyType = "ENSEMBL", ont="BP")
    save(geneOntologyEarlyGrowth,
         file = paste0("SavedData/", saveFolder, "/",tumorType,"/GO/geneOntologyEarlyGrowth.rda"))
  }
  
}

filterOnCombinedPAdj = function(data, pValues) {
  combinedStatistic = as.data.table(fisher.method(pValues))
  data[,c("combinedPValue","combinedPAdj") := combinedStatistic[,.(p.value,p.adj)]]
  return(data[combinedPAdj<0.05])
}

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRawData/",tumorIDs[i]))) {
    
    #Load in data
    load(paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T1vsT2.rda"))
    load(paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T2vsT3.rda"))
    load(paste0("SavedData/DESeqRawData/",tumorIDs[i],"/T1vsT3.rda"))
    
    generateBiologicallyRelevantResults(T1vsT2,T1vsT3,T2vsT3,tumorIDs[i],"DESeqRefinedData")
    
    #Load in data
    load(paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T1vsT2.rda"))
    load(paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T2vsT3.rda"))
    load(paste0("SavedData/LimmaVoomRawData/",tumorIDs[i],"/T1vsT3.rda"))
    
    dir.create("SavedData/LimmaVoomRefinedData")
    generateBiologicallyRelevantResults(T1vsT2,T1vsT3,T2vsT3,tumorIDs[i],"LimmaVoomRefinedData")
    
  }
  
}

##### Analyze Results Across Tumor Types ######

###How well do different tumor types bin?
consistentGrowthCounts = vector()
lateGrowthCounts = vector()
earlyGrowthCounts = vector()

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i]))) {
    
    load(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i],"/consistentGrowth.rda"))
    load(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i],"/lateGrowth.rda"))
    load(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i],"/earlyGrowth.rda"))
    
    consistentGrowthCounts[tumorIDs[i]] = dim(consistentGrowth)[1]
    lateGrowthCounts[tumorIDs[i]] = dim(lateGrowth)[1]
    earlyGrowthCounts[tumorIDs[i]] = dim(earlyGrowth)[1]
    
  }
  
}

barplot(consistentGrowthCounts, main = "Consistent Growth Bin Counts")
barplot(lateGrowthCounts, main = "Late Growth Bin Counts")
barplot(earlyGrowthCounts, main = "Early Growth Bin Counts")

### What genes do different tumor types have in common between bins?

# Function to return a frequency table based on a particular bin.
#   parentFolder should be the folder that contains the tumor type folders (e.g. "SavedData/DESeqRefinedData")
#   bin should be one of the three bins that genes are categorized in, such as "consistentGrowth"
generateIntersectionFrequencyTable = function(parentFolder, bin) {
  
  #Initialize the table
  intersectionFrequencyTable = data.table(Gene=character(),Count=integer())
  setkey(intersectionFrequencyTable,Gene)
  
  # Generate the frequency table
  for (i in 1:length(tumorIDs)) {
  
    #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
    if (file.exists(paste0(parentFolder, "/", tumorIDs[i]))) {
  
      #Load in the binned data
      binnedData = get(load(paste0(parentFolder , "/", tumorIDs[i], "/", bin, ".rda")))
  
      #Add any new genes to the frequency table.
      intersectionFrequencyTable = merge(intersectionFrequencyTable,binnedData[,.(Gene)], all=TRUE)
  
      #Increment gene counts based on genes in the current tumor type.
      intersectionFrequencyTable[intersectionFrequencyTable[,Gene] %in% binnedData[,Gene],
                                 Count := Count + 1]
      intersectionFrequencyTable[is.na(Count),
                                 Count := 1]
  
      #Indicate which genes were found in this tumor type.
      intersectionFrequencyTable[intersectionFrequencyTable[,Gene] %in% binnedData[,Gene],
                                 tumorIDs[i] := TRUE]
  
    }
  
  }
  
  # Replace NA's with FALSE.
  for (i in 1:length(tumorIDs)) {
  
    if (tumorIDs[i] %in% colnames(intersectionFrequencyTable)) {
      intersectionFrequencyTable[is.na(intersectionFrequencyTable[,tumorIDs[i],with = FALSE])[,1],tumorIDs[i] := FALSE]
    }
  
  }
  
  # Order the results by most common genes
  setorder(intersectionFrequencyTable,-Count)
  
  # Return the table.
  return(intersectionFrequencyTable)
  
}

#Generate the comparison tables
consistentGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","consistentGrowth")
lateGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","lateGrowth")
earlyGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","earlyGrowth")

#Save the results
dir.create("SavedData/LimmaVoomRefinedData/ComparisonTables")
save(consistentGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/consistentGrowth.rda")
save(lateGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/lateGrowth.rda")
save(earlyGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/earlyGrowth.rda")



##### Misc. #####

#Converting to "geneList" format for clusterProfiler
geneList = setNames(T1vsT3$log2FoldChange,T1vsT2IntersectT1vsT3$Gene)
geneList = sort(geneList, decreasing = TRUE)


### Converting between gene identifiers
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneInfo = as.data.table(getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'description'),
                  filters = "ensembl_gene_id",
                  values = T1vsT2IntersectT1vsT3[,Gene], 
                  mart = ensembl))


### How many patients do we have for each tumor type?
for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))) {
    
    # Load in Data
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    print(paste0(dim(sizeAndExpressionData)[1]," patients in ",tumorIDs[i]))
    
  }
}


###Compare sequencing counts for consistent growth genes in BRCA across other tumor types.

#Load in Data
load("SavedData/LimmaVoomRefinedData/BRCA/consistentGrowth.rda")

#Get a gene to test from BRCA consistentGrowth.
sampleGene = consistentGrowth$Gene[1]

#Compare across tumor types
for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))) {
    
    # Load in Data
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    load(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i],"/consistentGrowth.rda"))
    
    if(sampleGene %in% consistentGrowth$Gene) {
      print(paste0(tumorIDs[i], " displays consistent growth for gene ", sampleGene))
    }
    
    sizes = as.factor(sizeAndExpressionData$stage_event_tnm_categories)
    boxplot(log2(sizeAndExpressionData[[sampleGene]]) ~sizes, main = tumorIDs[i])
    
  }
}


### Compare Prostate Data to BRCA data
# Load in data
load("PRAD_T1vsT3.rda")
load("SavedData/LimmaVoomRawData/BRCA/T1vsT3.rda")

# Merge on Genes
PRADT1vsT3 = as.data.table(R_T1_VS_T3)
setnames(PRADT1vsT3,1,"Gene")
setkey(PRADT1vsT3,Gene)
BRCAT1vsT3 = T1vsT3[padj<0.05]
BRCAIntersectPRADT1vsT3 = merge(PRADT1vsT3,BRCAT1vsT3, suffixes = c("PRAD","BRCA"))
save(BRCAIntersectPRADT1vsT3,file = "BRCA_vs_PRAD_T1vsT3.rda")
dim(BRCAIntersectPRADT1vsT3)[1]