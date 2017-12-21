#Function to calculate percentages of taxa recovered and means and standard errors of these percentages

mean.counts = function(counts.list){
  percents = data.frame()
  results = data.frame()
  
  #Iterate through all taxonomic levels
  for(i in 1:length(counts.list)){
    taxa.counts = get(counts.list[i])
    #Iterate through all subjects  
    for(subject.ID in unique(taxa.counts$subject)){
      #Percent of taxa in S1 recovered in S2
      perc.S12 = taxa.counts[taxa.counts$comp == "S12_both" & taxa.counts$subject == subject.ID,]$count / sum(taxa.counts[taxa.counts$comp == "S12_both" & taxa.counts$subject == subject.ID,]$count + taxa.counts[taxa.counts$comp == "S12_S1only" & taxa.counts$subject == subject.ID,]$count)
      #Percent of taxa in S2 recovered in M
      perc.S2M = taxa.counts[taxa.counts$comp == "S2M_both" & taxa.counts$subject == subject.ID,]$count / sum(taxa.counts[taxa.counts$comp == "S2M_both" & taxa.counts$subject == subject.ID,]$count + taxa.counts[taxa.counts$comp == "S2M_S2only" & taxa.counts$subject == subject.ID,]$count)
      #Concatenate results into data frame with labels
      percents.iter = data.frame(c(perc.S12, perc.S2M), c("S12","S2M"), c(subject.ID, subject.ID), c(rep(counts.list[i],2)))
      #Add to results data frame    
      percents = rbind(percents, percents.iter)
    }}
  #Give column names to percents table
  colnames(percents) = c("perc","comp","subject","taxa") 
  
  for(i in 1:length(counts.list)){
    taxa.counts.name = counts.list[i]
    #Calculate means for each taxonomic level and comparison type
    mean.S12 = mean(percents[percents$taxa==taxa.counts.name & percents$comp == "S12",]$perc)
    mean.S2M = mean(percents[percents$taxa==taxa.counts.name & percents$comp == "S2M",]$perc)
    
    #Calculate standard errors for each taxonomic level and comparison type
    se.S12 = sd(percents[percents$taxa==taxa.counts.name & percents$comp == "S12",]$perc) / sqrt(length(percents[percents$taxa==taxa.counts.name & percents$comp == "S12",]$perc))
    se.S2M = sd(percents[percents$taxa==taxa.counts.name & percents$comp == "S2M",]$perc) / sqrt(length(percents[percents$taxa==taxa.counts.name & percents$comp == "S12",]$perc))
    #Put this iters results in a data frame
    results.iter = data.frame(c(mean.S12, mean.S2M), c(se.S12, se.S2M), c("S12","S2M"), c(rep(taxa.counts.name,2)))
    #Add to results data frame    
    results = rbind(results, results.iter)
  }
  #Create folder for output. Turn off warnings so if folder already exists, nothing happens.
  dir.create(file.path("output_mean.counts"), showWarnings = FALSE)
  #Remove previous output if exists
  if (file.exists("output_mean.counts/mean.counts.csv")) file.remove("output_mean.counts/mean.counts.csv")
  #Write results to .csv files
  write.table(results, "output_mean.counts/mean.counts.csv", sep=",", append=TRUE, col.names = FALSE)
}