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