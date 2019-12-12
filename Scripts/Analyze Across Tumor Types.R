##### Analyze Results Across Tumor Types ######

###How well do different tumor types bin?
displayBinningBarPlots = function(analysis) {
  consistentGrowthCounts = vector()
  lateGrowthCounts = vector()
  earlyGrowthCounts = vector()
  
  for (i in 1:length(tumorIDs)) {
    
    #Ensure that we have data to work with in the first place.  If not, skip this tumor type.
    if (file.exists(paste0("SavedData/LimmaVoomRefinedData/",tumorIDs[i]))) {
      
      load(paste0("SavedData/", analysis, "RefinedData/",tumorIDs[i],"/consistentGrowth.rda"))
      load(paste0("SavedData/", analysis, "RefinedData/",tumorIDs[i],"/lateGrowth.rda"))
      load(paste0("SavedData/", analysis, "RefinedData/",tumorIDs[i],"/earlyGrowth.rda"))
      
      consistentGrowthCounts[tumorIDs[i]] = dim(consistentGrowth)[1]
      lateGrowthCounts[tumorIDs[i]] = dim(lateGrowth)[1]
      earlyGrowthCounts[tumorIDs[i]] = dim(earlyGrowth)[1]
      
    }
    
  }
  
  barplot(consistentGrowthCounts, main = "Consistent Growth Bin Counts", sub = analysis)
  barplot(lateGrowthCounts, main = "Late Growth Bin Counts", sub = analysis)
  barplot(earlyGrowthCounts, main = "Early Growth Bin Counts", sub = analysis)
}

displayBinningBarPlots("DESeq")
displayBinningBarPlots("Limma Voom")



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


##DESeq
#Generate the comparison tables
consistentGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/DESeqRefinedData","consistentGrowth")
lateGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/DESeqRefinedData","lateGrowth")
earlyGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/DESeqRefinedData","earlyGrowth")

#Save the results
dir.create("SavedData/LimmaVoomRefinedData/ComparisonTables")
save(consistentGrowthComparisonTable, file = "SavedData/DESeqRefinedData/ComparisonTables/consistentGrowth.rda")
save(lateGrowthComparisonTable, file = "SavedData/DESeqRefinedData/ComparisonTables/lateGrowth.rda")
save(earlyGrowthComparisonTable, file = "SavedData/DESeqRefinedData/ComparisonTables/earlyGrowth.rda")


##Limma Voom
#Generate the comparison tables
consistentGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","consistentGrowth")
lateGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","lateGrowth")
earlyGrowthComparisonTable = generateIntersectionFrequencyTable("SavedData/LimmaVoomRefinedData","earlyGrowth")

#Save the results
dir.create("SavedData/LimmaVoomRefinedData/ComparisonTables")
save(consistentGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/consistentGrowth.rda")
save(lateGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/lateGrowth.rda")
save(earlyGrowthComparisonTable, file = "SavedData/LimmaVoomRefinedData/ComparisonTables/earlyGrowth.rda")