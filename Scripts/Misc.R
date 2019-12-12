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
PRADT1vsT3 = get(load("SavedData/LimmaVoomRawData/PRAD/T1vsT3.rda"))
BRCAT1vsT3 = get(load("SavedData/LimmaVoomRawData/BRCA/T1vsT3.rda"))

# Merge on Genes
PRADT1vsT3 = PRADT1vsT3[padj<0.05]
BRCAT1vsT3 = BRCAT1vsT3[padj<0.05]
BRCAIntersectPRADT1vsT3 = merge(PRADT1vsT3,BRCAT1vsT3, suffixes = c("PRAD","BRCA"))
save(BRCAIntersectPRADT1vsT3,file = "BRCA_vs_PRAD_T1vsT3.rda")
dim(BRCAIntersectPRADT1vsT3)[1]