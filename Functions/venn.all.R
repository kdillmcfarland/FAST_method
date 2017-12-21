#Function to identify taxa shared by S1 and S1 or S2 and M for each subject. No abundance cutoff.

venn.all = function(taxa.table, taxa.name){
  #Create empty data frames to hold results
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
    #Retrieve taxa names (column names) that are present in S1
    S12_S1 = colnames(OTU_S12[S1_ID,apply(OTU_S12[S1_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Retrieve taxa names (column names) that are present in S2
    S12_S2 = colnames(OTU_S12[S2_ID,apply(OTU_S12[S2_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Taxa shared by S1 and S2
    S12_both = intersect(S12_S1,S12_S2)
    ###S2 to M venns###
    #Subset data to only S2 and M samples    
    OTU_S2M = taxa.table[rownames(taxa.table) %in% c(S2_ID, M1_ID, M2_ID, M3_ID),]
    #Retrieve taxa names (column names) that are present in S2
    S2M_S2 = colnames(OTU_S2M[S2_ID,apply(OTU_S2M[S2_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Retrieve taxa names (column names) that are present in M1
    S2M_M1 = colnames(OTU_S2M[M1_ID,apply(OTU_S2M[M1_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Retrieve taxa names (column names) that are present in M2
    S2M_M2 = colnames(OTU_S2M[M2_ID,apply(OTU_S2M[M2_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Retrieve taxa names (column names) that are present in M3
    S2M_M3 = colnames(OTU_S2M[M3_ID,apply(OTU_S2M[M3_ID,], MARGIN=2,
                                          function(x) any(x>0))])
    #Taxa shared by S2 and M1
    S2M1_both = intersect(S2M_S2,S2M_M1)
    #Taxa shared by S2 and M2
    S2M2_both = intersect(S2M_S2,S2M_M2)
    #Taxa shared by S2 and M3
    S2M3_both = intersect(S2M_S2,S2M_M3)
    
    #Concatenate taxa lists into a data frame      
    taxa.results = data.frame(c(S12_both,S2M1_both,S2M2_both,S2M3_both),
                              c(rep("S12",length(S12_both)),rep("S2M1",length(S2M1_both)),rep("S2M2",length(S2M2_both)),rep("S2M3",length(S2M3_both))),
                              c(rep(subject.ID, length(c(S12_both,S2M1_both,S2M2_both,S2M3_both)))))
    #Add to results table  
    taxa.lists = rbind(taxa.lists, taxa.results)
  }
  #Create folder for output. Turn off warnings so if folder already exists, nothing happens.
  dir.create(file.path("output_venn.all"), showWarnings = FALSE)
  #Remove previous output if exists
  output_name = paste(taxa.name, "venn.all.list.csv", sep=".")
  output_path = paste("output_venn.all/", output_name, sep="")
  
  if (file.exists(output_path)) file.remove(output_path)
  #write results to .csv
  write.table(taxa.lists, output_path, sep=",", append=TRUE, col.names = FALSE)
}
