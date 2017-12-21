#Function to pull shared OTUs out of OTU table and calculate percent relative abundance of those shared OTUs

venn.perc = function(OTU.table){
  #Empty data frame to hold results
  percents = data.frame()
  #Iterate across all subjects
  for(subject.ID in unique(meta$Subject)){
    #Create IDs for each sample for that subject
    S1 = paste(subject.ID, "S1", sep="_")
    S2 = paste(subject.ID, "S2", sep="_")
    #Subset OTU table to sample of inerest  
    OTU.table.S1 = OTU.table[rownames(OTU.table) == S1,]
    OTU.table.S2 = OTU.table[rownames(OTU.table) == S2,]
    #Obtain list of OTUs shared by
    #S1 and S2
    OTU.list.S12 = as.vector(OTU.list[OTU.list$Subject==subject.ID & OTU.list$comparison == "S12",]$OTU)
    #S2 and M1
    OTU.list.S2M1 = as.vector(OTU.list[OTU.list$Subject==subject.ID & OTU.list$comparison == "S2M1",]$OTU)
    #S2 and M2
    OTU.list.S2M2 = as.vector(OTU.list[OTU.list$Subject==subject.ID & OTU.list$comparison == "S2M2",]$OTU)
    #S2 and M3
    OTU.list.S2M3 = as.vector(OTU.list[OTU.list$Subject==subject.ID & OTU.list$comparison == "S2M3",]$OTU)
    
    #Subset S1 data to just OTUs shared with S2
    OTU.S1S2 = OTU.table.S1[,colnames(OTU.table.S1) %in% OTU.list.S12]
    #Subset S2 data to just OTUs shared with each M sample
    OTU.S2M1 = OTU.table.S2[,colnames(OTU.table.S2) %in% OTU.list.S2M1]
    OTU.S2M2 = OTU.table.S2[,colnames(OTU.table.S2) %in% OTU.list.S2M2]
    OTU.S2M3 = OTU.table.S2[,colnames(OTU.table.S2) %in% OTU.list.S2M3]
    
    #Calculate percent relative abundance    
    S1S2.perc = sum(OTU.S1S2)/sum(OTU.table.S1)*100
    S2M1.perc = sum(OTU.S2M1)/sum(OTU.table.S2)*100
    S2M2.perc = sum(OTU.S2M2)/sum(OTU.table.S2)*100
    S2M3.perc = sum(OTU.S2M3)/sum(OTU.table.S2)*100
    
    #Put results in data frame
    perc = data.frame(c(S1S2.perc,S2M1.perc,S2M2.perc,S2M3.perc), c("S1S2","S2M1","S2M2","S2M3"), c("S1S2","S2M","S2M","S2M"), c(rep(subject.ID, 4)))
    #Add to result data
    percents = rbind(percents, perc)
  }
  #Create folder for output. Turn off warnings so if folder already exists, nothing happens.
  dir.create(file.path("output_venn.perc"), showWarnings = FALSE)
  #Remove previous output if exists
  if (file.exists("output_venn.perc/venn.perc.csv")) file.remove("output_venn.perc/venn.perc.csv")
  #Write to .csv
  write.table(percents, "output_venn.perc/venn.perc.csv", sep=",", append=TRUE, col.names = FALSE)
}