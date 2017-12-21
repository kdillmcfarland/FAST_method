#To obtain lists of beta-diversity values of mouse samples compared to samples from their matched subject donor (DONOR) or compared to any other subject in the dataset (OTHER).

#Create distance matrix from `Subject` ID comparisons so that you can pull out donor-mouse for one subject (ID distance = 0)

donor.other = function(subject, groups, beta){
  #Empty results table
  match.unmatch.results = data.frame()
  
  #Create vector of subject IDs
  ID.vec = as.vector(subject)
  #Calculate distance matrix comparing IDs. A value of 0 means the two samples have the same subject ID (DONOR) and any nonzero values means they are from different subjects (OTHER)
  ID.dist=as.matrix(dist(ID.vec, method="euclidean"))
  
  #Create matrix of group comparisons so that you can pull out specifc samples, specifically S1 or S2 vs. F
  #Create vector of groups
  Link.vec = as.vector(groups)
  #Using a loop, create a matrix of link ID comparisons
  Link.mat = matrix(numeric(0), length(ID.vec), length(ID.vec))
  for(j in 1:dim(Link.mat)[1] ){
    for(k in 1:dim(Link.mat)[2]){
      Link.mat[j, k] = paste(Link.vec[j], Link.vec[k], sep = "")
    }}

  #Pull out pairs of interest
  ###Mouse fecal vs. human sample: Link.mat == "FS1" or "FS2"
  ###From the same subject: ID.dist == 0
  ###From different subject: ID.dist =/= 0
  
  #Loop across beta-diversity metrics
  for(metric in 1:length(beta)){
    #Obtain distance matrix
    dist.name = paste(beta[metric], ".dist", sep="")
    distances = get(dist.name)
    
    #DONOR
    match.FS1 = distances[ID.dist == 0 & Link.mat == "FS1"]
    match.FS2 = distances[ID.dist == 0 & Link.mat == "FS2"]
      #Combine into 1 vector
      match = c(match.FS1, match.FS2)
    
    #OTHER
    unmatch.FS1 = distances[ID.dist != 0 & Link.mat == "FS1"]
    unmatch.FS2 = distances[ID.dist != 0 & Link.mat == "FS2"]
      #Combine into 1 vector
      unmatch = c(unmatch.FS1, unmatch.FS2)
    
    #Add to results table
    results.table = data.frame(c(match, unmatch),
                               c(rep("match", length(match)), rep("unmatch", length(unmatch))),
                               c(rep(beta[metric], length(match)+length(unmatch))))
    match.unmatch.results = rbind(match.unmatch.results, results.table)
  } 
  #Create folder for output. Turn off warnings so if folder already exists, nothing happens.
  dir.create(file.path("output_donor.other"), showWarnings = FALSE)
  #Remove previous output if exists
  if (file.exists("output_donor.other/donor.other.csv")) file.remove("output_donor.other/donor.other.csv")
  #Write data to file
  write.table(match.unmatch.results, "output_donor.other/donor.other.csv", sep=",", append=TRUE, col.names = FALSE)
}