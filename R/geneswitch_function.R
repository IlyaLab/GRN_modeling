#function for turning gene mutations on/off

#combine the table 7 rows into columns, so one column is equal to one patient 
#patient_profile_df <- s7_df %>%
#dplyr::group_by(labId) %>%
#dplyr::summarise(symbol = paste(symbol, collapse = ","))

#FOR NOW IM MAKING A TEST DATAFRAME
#patient_profile_test_df <- patient_profile_df[1:5, ]

#create the dictonary
genes <- Dict$new(  
  FLT3 = "time.series[1,col] <- 1" ,
  CEBPA = "time.series[2,col] <- 0" ,
  RUNX1 = "time.series[3,col] <- 1" ,			
  RAD21 = "time.series[4,col] <- 0" ,
  #AKT 
  NRAS = "time.series[6,col] <- 1",
  BRAF = "time.series[7,col] <- 1",
  #AMPK,
  NPM1 = "time.series[9,col] <- 0",
  MEK = "time.series[10,col] <- 0",
  FOXO = "time.series[11,col] <- 1",
  FBXW7 = "time.series[12,col] <- 1",  
  #ERK,
  MYC = "time.series[14,col] <- 1",
  IDH2 = "time.series[15,col] <- 1", 
  IDH1 = "time.series[16,col] <- 1",
  #OXO2,
  CDKN2A = "time.series[18,col] <- 0", 
  TET2 = "time.series[19,col] <- 0",
  MDM2 = "time.series[20,col] <- 1", 
  WT1 = "time.series[21,col] <- 1",			
  TP53 = "time.series[22,col] <- 0",				
  BCL2 = "time.series[23,col] <- 1",
  PTPN11 = "time.series[24,col] <- 1",	
  .class = "character",
  .overwrite = TRUE)

#genes["FLT3"]= time.series[j,col] <- 1

#read the mutation profile
geneswitch <- function(mutation_profile) {
  
  #seperate the values by the commas
  mutation_profile <- unlist(strsplit(mutation_profile,","))
  
  #remove repeated values
  mutation_profile <- unique(mutation_profile)
  
  #save as a list
  mutation_list <- as.list(mutation_profile)
  
  #so now it should look like this (example)
  #List values: 
    #"FLT3"
    #"DNMT3A"
    #"NPM1"
  
  #now we loop through the list 
  for (p in mutation_list) {
    cat(genes[p], sep="\n")
  }}


#how the function will be used -> can execute
#read a row from table
#for (i in 1:nrow(patient_profile_test_df)){
#patient_id <- patient_profile_test_df[i, "labId"]          #stores the labID as a patient_id variable for each row
#mutation_profile <- as.vector(patient_profile_test_df$symbol[i]) #stores the mutations as a mutation_profile variable for each row

#cat(geneswitch(mutation_profile))}
