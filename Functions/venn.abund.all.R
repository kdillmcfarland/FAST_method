#Function to identify taxa shared by S1 and S1 or S2 and M for each subject. Taxon must be at least 0.1% relaive abundance in at least one sample in the comparison for each subject.

venn.abund.all = function(taxa.table, taxa.name){
  #Create empty data frames to hold results
  counts = data.frame() 
  taxa.lists = data.frame()
  
  #Iterate across all subjects
  for(subject.ID in unique(meta$Subject)){
    #Create name IDs for each sample
    S1_ID = paste(subject.ID, "S1", sep="_")
    S2_ID = paste(subject.ID, "S2", sep="_")
    M1_ID = paste(subject.ID, "F1", sep="_")
    M2_ID = paste(subject.ID, "F2", sep="_")
    M3_ID = paste(subject.ID, "F3", sep="_")
    
    ###S1 vs. S2 venns###
    #Subset data to only S1 and S2 samples
    OTU_S12 = taxa.table[rownames(taxa.table) %in% c(S1_ID, S2_ID),]
    #Calculate percent relative abundance
    OTU_S12_perc = OTU_S12/rowSums(OTU_S12)*100
    #Subset data only taxa at least 0.1% relative abundance in S1 and/or S2
    OTU_S12_abund = OTU_S12[, apply(OTU_S12_perc, MARGIN=2, function(x) any(x > 0.1))]
    #Retrieve abundant taxa names (column names) that are present in S1
    S12_S1 = colnames(OTU_S12_abund[S1_ID,apply(OTU_S12_abund[S1_ID,], MARGIN=2, function(x) any(x>0))])
    #Retrieve abundant taxa names (column names) that are present in S2
    S12_S2 = colnames(OTU_S12_abund[S2_ID,apply(OTU_S12_abund[S2_ID,], MARGIN=2, function(x) any(x>0))])
    #Abundant taxa shared by S1 and S2
    S12_both = intersect(S12_S1,S12_S2)
    #Abundant taxa only in S1
    S12_S1only = setdiff(S12_S1,S12_S2)
    #Abundant taxa only in S2
    S12_S2only = setdiff(S12_S2,S12_S1)
    
    ###S2 to M venns###
    #Subset data to only S2 and M samples
    OTU_S2M = taxa.table[rownames(taxa.table) %in% c(S2_ID, M1_ID, M2_ID, M3_ID),]
    #Calculate percent relative abundance
    OTU_S2M_perc = OTU_S2M/rowSums(OTU_S2M)*100
    #Subset data only taxa at least 0.1% relative abundance in S2 and/or M
    OTU_S2M_abund = OTU_S2M[, apply(OTU_S2M_perc, MARGIN=2, function(x) any(x > 0.1))]
    #Retrieve abundant taxa names (column names) that are present in S2
    S2M_S2 = colnames(OTU_S2M_abund[S2_ID,apply(OTU_S2M_abund[S2_ID,], MARGIN=2, function(x) any(x>0))])
    #Retrieve abundant taxa names (column names) that are present in M1
    S2M_M1 = colnames(OTU_S2M_abund[M1_ID,apply(OTU_S2M_abund[M1_ID,], MARGIN=2, function(x) any(x>0))])
    #Retrieve abundant taxa names (column names) that are present in M2
    S2M_M2 = colnames(OTU_S2M_abund[M2_ID,apply(OTU_S2M_abund[M2_ID,], MARGIN=2, function(x) any(x>0))])
    #Retrieve abundant taxa names (column names) that are present in M3
    S2M_M3 = colnames(OTU_S2M_abund[M3_ID,apply(OTU_S2M_abund[M3_ID,], MARGIN=2, function(x) any(x>0))])
    #Combine all M taxa lists and keep only uniques
    S2M_M = unique(c(S2M_M1, S2M_M2, S2M_M3))   
    #Abundant taxa shared by S2 and M
    S2M_both = intersect(S2M_S2,S2M_M)
    #Abundant taxa only in S2
    S2M_S2only = setdiff(S2M_S2,S2M_M)
    #Abundant taxa only in M
    S2M_Monly = setdiff(S2M_M,S2M_S2)
    
    #Combine all taxa lists into a data frame
    taxa.results = data.frame(
      c(S12_both,S12_S1only,S12_S2only,S2M_both,S2M_S2only,S2M_Monly),
      c(rep("S12_both",length(S12_both)),rep("S12_S1only",length(S12_S1only)),rep("S12_S2only",length(S12_S2only)),rep("S2M_both",length(S2M_both)),rep("S2M_S2only",length(S2M_S2only)),rep("S2M_Monly",length(S2M_Monly))),
      c(rep(subject.ID, length(c(S12_both,S12_S1only,S12_S2only,S2M_both,S2M_S2only,S2M_Monly)))))
    
    #Count the number of taxa in each list and combine into data frame  
    count.results = data.frame(
      c(length(S12_both), length(S12_S1only), length(S12_S2only), length(S2M_both), length(S2M_S2only), length(S2M_Monly)),
      c("S12_both","S12_S1only","S12_S2only","S2M_both","S2M_S2only","S2M_Monly"),
      c(rep(subject.ID, 6)))
    
    #Append this subject's results to results tables
    counts = rbind(counts, count.results)
    taxa.lists = rbind(taxa.lists, taxa.results)
  }
  #Create folder for output. Turn off warnings so if folder already exists, nothing happens.
  dir.create(file.path("output_venn.abund.all"), showWarnings = FALSE)
  
  #Remove previous outputs if exist
  output_name1 = paste(taxa.name, "venn.abund.all.list.csv", sep=".")
  output_path1 = paste("output_venn.abund.all/", output_name1, sep="")
  if (file.exists(output_path1)) file.remove(output_path1)

  output_name2 = paste(taxa.name, "venn.abund.all.count.csv", sep=".")
  output_path2 = paste("output_venn.abund.all/", output_name2, sep="")
  if (file.exists(output_path2)) file.remove(output_path2)  
  
    #Write results to .csv files
  write.table(taxa.lists, output_path1, sep=",", append=TRUE, col.names = FALSE)
  write.table(counts, output_path2, sep=",", append=TRUE, col.names = FALSE)
}