### Created by Tazein Shah 
### Last modified on 7/25/2023
### This code takes the FLT3-NPM1-IDH2 BNM from the Palma paper and personalizes it per patient profile (from s7_df) 
### Some of the mutations are inaccurate/missing

#creating the FLT3-NPM1-IDH2 BNM 
#load BoolNet
library(BoolNet)
library(dplyr)

#creating the file for the  FLT3-NPM1-IDH
fil <- tempfile(pattern = "testNet")
sink(fil)
cat("targets, factors\n")
cat("FLT3, FLT3\n")
cat("CEBPA, !FLT3\n")
cat("RUNX1, FLT3\n")
cat("RAD21, !RUNX1\n")
cat("PTPN11, FLT3\n")
cat("AKT, AKT\n")
cat("NRAS, PTPN11\n")
cat("BRAF, NRAS\n")
cat("AMPK, AMPK\n") 
cat("NPM1, NPM1\n")
cat("MEK, BRAF\n")
cat("FOXO, !AKT & AMPK\n")
cat("FBXW7, NPM1\n")
cat("ERK, MEK\n")
cat("MYC, !(RUNX1 | CEBPA | FBXW7) & ERK\n")
cat("IDH2, IDH2 \n")
cat("IDH1, FOXO\n")
cat("OXO2 , IDH1 & IDH2\n")
cat("CDKN2A, !MYC & NPM1\n")
cat("TET2, OXO2 & AMPK\n")
cat("MDM2, !CDKN2A & TP53\n")
cat("WT1, TET2\n")
cat("TP53, !MDM2\n")
cat("BCL2, !TP53 & ERK\n")
cat("Proliferation, !(FOXO | WT1 | CDKN2A | TP53) & MYC\n")
cat("Differentiation, RUNX1 & CEBPA\n")
cat("Apoptosis, !BCL2 & TP53\n")

sink()

#this loads the network and saves it as flt3_npm1_idh_net
flt3_idh_net <- loadNetwork(fil)
print(flt3_idh_net)

#now we move on to table 7
#combine the table 7 rows into columns, so one column is equal to one patient 
patient_profile_df <- s7_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#FOR NOW IM MAKING A TEST DATAFRAME
patient_profile_test_df <- patient_profile_df[1:5, ]

#In subgroup diagram but cannot find mutation 
#AKT 
#AMPK
#MEK
#ERK
#OXO2

fixIndices <- c()
values <- c()
score_list <- list() 

#read a row from table
for (i in 1:nrow(patient_profile_test_df)){
  patient_id <- patient_profile_test_df[i, "labId"]          #stores the labID as a patient_id variable for each row
  mutation_profile <- patient_profile_test_df[i, "symbol"]   #stores the mutations as a mutation_profile variable for each row
  
  if (any(grepl("FLT3", mutation_profile))) {                #if the mutation is present, it moves on to the loop
    fixIndices <- c(fixIndices, "FLT3")                    #adds mutation to fixIndices
    values <- c(values, 1)}                                #adds the mutation value to values
  
  if (any(grepl("CEBPA", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "CEBPA")               
    values <- c(values, 0)}		
  
  if (any(grepl("RUNX1", mutation_profile))) {       			#low classification consistency         
    fixIndices <- c(fixIndices, "RUNX1")               
    values <- c(values, 1)}
  
  if (any(grepl("RAD21", mutation_profile))) {                #low classification consistency
    fixIndices <- c(fixIndices, "RAD21")               
    values <- c(values, 0)}
  
  if (any(grepl("PTPN11", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "PTPN11")               
    values <- c(values, 1)}
  
  if (any(grepl("NRAS", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "NRAS")               
    values <- c(values, 1)}
  
  if (any(grepl("NPM1", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "NPM1")               
    values <- c(values, 0)}
  
  if (any(grepl("MEK", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "MEK")               
    values <- c(values, 0)}
  
  if (any(grepl("MYC", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "MYC")               
    values <- c(values, 1)}
  
  if (any(grepl("IDH2", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "IDH2")               
    values <- c(values, 1)}
  
  if (any(grepl("IDH1", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "IDH1")               
    values <- c(values, 1)}
  
  if (any(grepl("CDKN2A", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "CDKN2A")               
    values <- c(values, 0)}
  
  if (any(grepl("TET2", mutation_profile))) {                 
    fixIndices <- c(fixIndices, "TET2")               
    values <- c(values, 0)}
  
  if (any(grepl("WT1", mutation_profile))) {                 #low classification consistency (change to 1)
    fixIndices <- c(fixIndices, "WT1")               
    values <- c(values, 0)}
  
  if (any(grepl("TP53", mutation_profile))) {                 #low classification consistency
    fixIndices <- c(fixIndices, "TP53")               
    values <- c(values, 0)}
  
  if (any(grepl("BRAF", mutation_profile))) {                 #from COSMIC
    fixIndices <- c(fixIndices, "BRAF")               
    values <- c(values, 1)}
  
  if (any(grepl("FOXO", mutation_profile))) {                 #from COSMIC
    fixIndices <- c(fixIndices, "FOXO")               
    values <- c(values, 1)}
  
  if (any(grepl("FBXW7", mutation_profile))) {                #from COSMIC
    fixIndices <- c(fixIndices, "FBXW7")               
    values <- c(values, 1)}
  
  if (any(grepl("MDM2", mutation_profile))) {                 #from COSMIC
    fixIndices <- c(fixIndices, "MDM2")               
    values <- c(values, 1)}
  
  if (any(grepl("BCL2", mutation_profile))) {                 #from COSMIC
    fixIndices <- c(fixIndices, "BCL2")               
    values <- c(values, 1)}
  
  #save this profile as it's own BNM figure
  patient_network <- fixGenes(flt3_idh_net, fixIndices= fixIndices, values= values)
  assign(paste0("network_", patient_id), patient_network)
  score_list[[i]] <- patient_network
  
  #get attractors
  attr <- getAttractors(patient_network, method = c("sat.exhaustive"))
  assign(paste0("attr_", patient_id), attr)
  score_list[[i]] <- attr
  
  print(paste0("Finish attrs for ", patient_id$labId))
  
  #get matrices
  matrix <- as.data.frame(plotAttractors(attr)$'4')
  assign(paste0("matrix_", patient_id), matrix)
  score_list[[i]] <- matrix
  
  print(paste0("Finish ", patient_id$labId))
  
  #reset for next loop 
  fixIndices <- c()
  values <- c()
}


