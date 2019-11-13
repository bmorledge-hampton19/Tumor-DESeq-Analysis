#Need to align gene expression to tumor size using patient barcode.
#BUT sometimes multiple samples are taken from one patient.
#Each patient is only associated with at most one tumor size, so those with multiple samples need to be omitted.
#Also, we only care about tumor sizes T1, T2, and T3.
#(Thanks Pete!)
#Then, run DESeq analysis on the data.



##### Set Up the Environment #####

#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

#Include libraries
library(TCGAbiolinks)
library(DESeq2)
library(clusterProfiler)
library(SummarizedExperiment)
library(data.table)
library(org.Hs.eg.db)
library(biomaRt)
library(GO.db)

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
    expressionData = melt(expressionData, id.vars = "bcr_patient_barcode", 
                          variable.name = "Gene", value.name = "Counts")
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



##### Refine Results #####

# Create the necessary Directory for saving data.
dir.create("SavedData/DESeqRefinedData")

# A nice function for creating results from different comparisons
createRefinedResultsTable = function(DESeqResults, genes, size1, size2) {
  
  #Get the results for the indicated comparison
  results = results(DESeqResults, contrast = c("stage_event_tnm_categories",size1,size2), 
                    parallel = TRUE, BPPARAM = SnowParam(7))
  
  #Convert the results to a data.table and add the column for gene ID's.
  DESeqResultsTable = as.data.table(results)[,Gene := genes]
  setcolorder(DESeqResultsTable,c(7,1:6))
  
  #Remove insignificant changes and small changes.
  DESeqResultsTable = DESeqResultsTable[padj < 0.05 & abs(log2FoldChange) > 0.58]
  
  return(DESeqResultsTable)
  
}

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSR.rda"))) {
    
    #Load in DESeq Results
    load(paste0("SavedData/DESeqRawData/TCGA-",tumorIDs[i],"-DSR.rda"))
    load(paste0("SavedData/SizeAndGeneExpression/TCGA-",tumorIDs[i],"-CGE.rda"))
    
    #Create the comparison tables
    T1vsT2 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T2")
    T1vsT3 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T1", "T3")
    T2vsT3 = createRefinedResultsTable(DESeqResults, colnames(sizeAndExpressionData)[-(1:2)], "T2", "T3")
    
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
    dir.create(paste0("SavedData/DESeqRefinedData/",tumorIDs[i]))
    
    #Save the tables
    save(T1vsT2, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2.rda"))
    save(T1vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT3.rda"))
    save(T2vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T2vsT3.rda"))
    save(T1vsT2IntersectT1vsT3, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2_Intersect_T1vsT3.rda"))
    
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

generateDotPlot = function(genes, geneOntologyCategory) {
  geneOntologyResults = enrichGO(genes, 'org.Hs.eg.db',
                                 keyType = "ENSEMBL", ont=geneOntologyCategory, pvalueCutoff = 0.1)
  dotplot(geneOntologyResults)
}

for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRefinedData/",tumorIDs[i]))) {
    
    #Load in data
    load(paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2.rda"))
    load(paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T2vsT3.rda"))
    load(paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT3.rda"))
    
    #Find genes that increase in expression at each increase in tumor size.
    consistentGrowth = merge(T1vsT2[log2FoldChange<0],T2vsT3[log2FoldChange<0],by = "Gene",
                             suffixes = c("T1vsT2","T2vsT3"))
    
    #Sanity Check
    if(length(intersect(consistentGrowth[,Gene],T1vsT3[,Gene])) != length(consistentGrowth)) {
      print("A gene in the set T1vsT3 (Increased expression) is not found in the intersection between
            T1vsT2 and T2vsT3 (Also increased expression) This should be impossible...")
    }
    
    #Find genes that increase in expression from T1 to T3, but do something weird inbetween...
    inconsistentGrowth = setdiff(T1vsT3[log2FoldChange<0],consistentGrowth)
    
    #These genes increase in expression very slowly, or only between T2 and T3.
    veryGradualOrLateGrowth = inconsistentGrowth[!(Gene %in% T1vsT2[log2FoldChange>0])]
    #These genes are expressed less in T2 tumors but more in T3 tumors! (Compared to T1)
    T2IndependentGrowth = inconsistentGrowth[(Gene %in% T1vsT2[log2FoldChange>0])]
  
    #These genes are expressed highly in T2 tumors, but expression is brought back down in T3 tumors.
    earlyGrowth = merge(T1vsT2[log2FoldChange<0],T2vsT3[log2FoldChange>0],
                        suffixes = c("T1vsT2","T2vsT3"))
   
    #Save the results.
    save(consistentGrowth, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/consistentGrowth.rda"))
    save(veryGradualOrLateGrowth, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/verGradualOrLateGrowth.rda"))
    save(T2IndependentGrowth, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T2IndependentGrowth.rda"))
    save(earlyGrowth, file = paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/earlyGrowth.rda"))
    
    #Generate dotplots.
    generateDotPlot(consistentGrowth[,Gene],"BP")
    generateDotPlot(veryGradualOrLateGrowth[,Gene],"BP")
    generateDotPlot(T2IndependentGrowth[,Gene],"BP")
    generateDotPlot(earlyGrowth[,Gene],"BP")
    
    #Generate heatmaps.
    generateHeatMap(consistentGrowth[,Gene],T1vsT3[Gene %in% consistentGrowth[,Gene],log2FoldChange])
    generateHeatMap(veryGradualOrLateGrowth[,Gene],T1vsT3[Gene %in% veryGradualOrLateGrowth[,Gene],log2FoldChange])
    generateHeatMap(T2IndependentGrowth[,Gene])
    generateHeatMap(earlyGrowth[,Gene])
    
  }
  
}

##### Analyze Results Across Tumor Types ######

# Make the base frequency table
intersectionFrequencyTable = data.table(Gene=character(),Count=integer())
setkey(intersectionFrequencyTable,Gene)

# Generate the frequency table
for (i in 1:length(tumorIDs)) {
  
  #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
  if (file.exists(paste0("SavedData/DESeqRefinedData/",tumorIDs[i]))) {
    
    #Load in intersection table
    load(paste0("SavedData/DESeqRefinedData/",tumorIDs[i],"/T1vsT2_Intersect_T1vsT3.rda"))
    
    #Add any new genes to the frequency table.
    intersectionFrequencyTable = merge(intersectionFrequencyTable,T1vsT2IntersectT1vsT3[,.(Gene)], all=TRUE)
    
    #Increment gene counts based on genes in the current tumor type.
    intersectionFrequencyTable[intersectionFrequencyTable[,Gene] %in% T1vsT2IntersectT1vsT3[,Gene],
                               Count := Count + 1]
    intersectionFrequencyTable[is.na(Count),
                               Count := 1]
    
    #Indicate which genes were found in this tumor type.
    intersectionFrequencyTable[intersectionFrequencyTable[,Gene] %in% T1vsT2IntersectT1vsT3[,Gene],
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

#Save the result
save(intersectionFrequencyTable, file = "SavedData/DESeqRefinedData/IntersectionFrequencyTable.rda")

##### Experimenting... #####

genesWithFoldChanges = setNames(T1vsT2IntersectT1vsT3[,T1vsT2.log2FoldChange],T1vsT2IntersectT1vsT3[,Gene])

#hgnc_symbol?
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneInfo = as.data.table(getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'description'),
                  filters = "ensembl_gene_id",
                  values = T1vsT2IntersectT1vsT3[,Gene], 
                  mart = ensembl))

geneInfo = as.data.table(getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'description'),
                  filters = "hgnc_symbol",
                  values = "LACRT", 
                  mart = ensembl))
T1Mean = colMeans(sizeAndExpressionData[stage_event_tnm_categories == "T1", 
                                        geneInfo[,ensembl_gene_id], with = FALSE])
T2Mean = colMeans(sizeAndExpressionData[stage_event_tnm_categories == "T2", 
                                        geneInfo[,ensembl_gene_id], with = FALSE])
T3Mean = colMeans(sizeAndExpressionData[stage_event_tnm_categories == "T3", 
                                        geneInfo[,ensembl_gene_id], with = FALSE])
T1vsT3[Gene == geneInfo[,ensembl_gene_id],log2FoldChange]
log(T1Mean/T3Mean,2)
T1vsT3[Gene == geneInfo[,ensembl_gene_id],baseMean]
colMeans(sizeAndExpressionData[,geneInfo[,ensembl_gene_id], with = FALSE])
#Normalizing for sequence depth?

GOInfo = as.data.table(getBM(attributes=c('ensembl_gene_id',"go_id"),
                              filters = "ensembl_gene_id",
                              values = T1vsT2IntersectT1vsT3[,Gene], 
                              mart = ensembl))
GOTerms = Term(GOTERM)
GOTerms = data.table(GOID = names(GOTerms), GODescription = GOTerms)
GOInfo = merge(GOInfo,GOTerms,by.x = "go_id",by.y = "GOID")

GOInfo = as.data.table(getBM(attributes=c('ensembl_gene_id',"go_id"),
                             filters = c("ensembl_gene_id","go_id"),
                             values = intersectionFrequencyTable[Count >= 4,Gene], 
                             mart = ensembl))
GOInfo = merge(GOInfo,GOTerms,by.x = "go_id",by.y = "GOID")

spindleGenes = GOInfo[grepl("spindle", GOInfo[,GODescription])]
mitoticGenes = GOInfo[grepl("mitotic|mitosis", GOInfo[,GODescription])]

clinical = clinical[grepl("T1|T2|T3", clinical[,stage_event_tnm_categories])]

test = bitr(T1vsT2IntersectT1vsT3[,Gene], fromType = "ENSEMBL",
            toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

test = enrichGO(T1vsT2IntersectT1vsT3[,Gene], 'org.Hs.eg.db',
                keyType = "ENSEMBL", ont="BP", pvalueCutoff=0.05)
barplot(test)
testTable = as.data.table(test)

test2 = enrichGO(T1vsT2IntersectT1vsT3[,Gene], 'org.Hs.eg.db',
                 keyType = "ENSEMBL", ont="BP", pvalueCutoff=0.1)
barplot(test2)
testTable2 = as.data.table(test2)

# Doesn't Work?
test3 = gseGO(T1vsT2IntersectT1vsT3[order(-Gene),Gene], OrgDb = 'org.Hs.eg.db',
              keyType = "ENSEMBL", ont="BP", pvalueCutoff=1)
barplot(test3)
testTable2 = as.data.table(test3)

test4 = enrichGO(T1vsT2IntersectT1vsT3[,Gene], 'org.Hs.eg.db',
                 keyType = "ENSEMBL", ont="MF", pvalueCutoff=0.1)
barplot(test4)
testTable2 = as.data.table(test4)

test5 = setReadable(test4, 'org.Hs.eg.db', 'ENSEMBL')
cnetplot(test5, foldChange = genesWithFoldChanges, circular = TRUE, colorEdge = TRUE)
heatplot(test5, foldChange = genesWithFoldChanges)

geneOntologyResults = enrichGO(veryGradualOrLateGrowth[,Gene], 'org.Hs.eg.db', 
                               keyType = "ENSEMBL", ont="CC", pvalueCutoff = 0.1)
dotplot(geneOntologyResults)
genesWithFoldChanges = setNames(T1vsT2IntersectT1vsT3[,T1vsT2.log2FoldChange],T1vsT2IntersectT1vsT3[,Gene])
heatplot(setReadable(geneOntologyResults, "org.Hs.eg.db", "ENSEMBL"), foldChange = genesWithFoldChanges)

results

# DAVID? Cluster Profiler?  GG Plot?  SEUART?
# Biomart for converting to gene name/description