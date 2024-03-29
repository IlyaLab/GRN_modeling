### Created by Tazein Shah 
### Last modified on 7/25/2023
### This code takes the FLT3-NPM1-IDH2 BNM from the Palma paper and personalizes it per patient profile (from s7_df)
### Uses a function called scorePersonal() to knock in/out genes (main difference from idh2_personalized_bnm.R) 
### Some of the mutations are inaccurate/missing


#load needed packages
library(BoolNet)
library(reshape2)
library(dplyr)
library(zoo)
library(Dict)
library(ggplot2)
library(ggrepel)
#library(writexl)
#library(readxl)

#importing the excel sheets 
#s7_df <- read_excel("(insert_path_to_excel_sheet)/s7_table.xlsx")
#s5_df <- read_excel("(insert_path_to_excel_sheet)/s5_table.xlsx")

#combine the table 7 rows into columns, so one row is equal to one patient 
patient_profile_df <- s7_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#FOR NOW IM MAKING A TEST DATAFRAME
#chooses the first 200 samples
  patient_profile_test_df <- patient_profile_df[1:10, ]
#chooses a random 200 samples of the dataframe
  #patient_profile_test_df <- patient_profile_df[sample(nrow(patient_profile_df), 200), ]

#create the dictonary
genes_values <- Dict$new(  
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

n.obs <- 20500 # number of observation aka steps 
n.half <- (n.obs/2) #half the steps (for plot making)
p <- 24 # number of variables

#how the mutation profile will be formatted
mutation_profile <- as.vector(patient_profile_test_df$symbol[2]) 

#this is a threshold to flip the node 
chance.to.flip <- 0.01

############# scorePersonal() function #############

scorePersonal <- function(mutation_profile) { 
  
  #separate the values by the commas
  mutation_profile <- (unlist(strsplit(mutation_profile,",")))
  
  #remove repeated values
  mutation_profile <- unique(mutation_profile)
  
  #save as a list
  #mutation_list <- as.list(mutation_profile)
  
  #so now it should look like this (example)
  #List values: 
  #"FLT3"
  #"DNMT3A"
  #"NPM1"
  
  print(paste0(mutation_profile))
  print("preparing to create time series")
  
  # 24 time-series of noise for 24 variables
  zeta <- matrix(runif(n.obs * p, 0, 1), ncol = p)
  
  zeta <- t(zeta) # adjust zeta to a wide-format
  
  time.series <- matrix(rep(0, p * n.obs), nrow =p , ncol =  n.obs)
  
  time.series[,1] <- runif(p,0, 1) # initial value of the three variables were randomly chosen
  
  time.series[,1] <-ifelse(time.series[,1] >= .5, TRUE,FALSE) # if the randomly generated number is greater than .5, then we convert the value to TRUE, otherwise to FALSE
  
  for (col in 2:n.obs){
    # simulate each time step based on the previous time step and Boolean functions
    time.series[1,col] <-time.series[1,col-1] 									            #FLT3= FLT3
    time.series[2,col] <-!time.series[1,col-1] 									            #CEBPA = !FLT3
    time.series[3,col] <-time.series[1,col-1] 									            #RUNX1 = FLT3
    time.series[4,col] <-!time.series[3,col-1] 									            #RAD21 = !RUNX1
    time.series[5,col] <-time.series[5,col-1] 									            #AKT = AKT
    time.series[6,col] <-time.series[24,col-1] 									            #NRAS = PTPN11
    time.series[7,col] <-time.series[6,col-1] 									            #BRAF = NRAS
    time.series[8,col] <-time.series[8,col-1] 									            #AMPK = AMPK
    time.series[9,col] <-time.series[9,col-1] 								            	#NPM1 = NPM1
    time.series[10,col] <-time.series[7,col-1] 							            		#MEK = BRAF
    time.series[11,col] <-!time.series[5,col-1] & time.series[8,col-1] 						#FOXO = !AKT & AMPK
    time.series[12,col] <-time.series[9,col-1] 									            #FBXW7 = NPM1
    time.series[13,col] <-time.series[10,col-1] 									          #ERK = MEK
    time.series[14,col] <-!(time.series[3,col-1] | time.series[2,col-1] | time.series[12,col-1]) & time.series[13,col-1] 	#MYC = !(RUNX1 | CEBPA | FBXW7) & ERK
    time.series[15,col] <-time.series[15,col-1] 									          #IDH2 = IDH2
    time.series[16,col] <-time.series[11,col-1] 									          #IDH1 = FOXO
    time.series[17,col] <-time.series[16,col-1] & time.series[15,col-1]							#OXO2 = IDH1 & IDH2
    time.series[18,col] <-!time.series[14,col-1]  & time.series[9,col-1]						#CDKN2A = !MYC & NPM1
    time.series[19,col] <-time.series[17,col-1] & time.series[8,col-1]							#TET2 = OXO2 & AMPK
    time.series[20,col] <- !time.series[18,col-1] & time.series[22,col-1] 					#MDM2 = !CDKN2A & TP53
    time.series[21,col] <-time.series[19,col-1] 									          #WT1 = TET2
    time.series[22,col] <-!time.series[20,col-1] 									          #TP53 = !MDM2
    time.series[23,col] <-!time.series[22,col-1] & time.series[13,col-1]						#BCL2 = !TP53 & ERK
    time.series[24,col] <-time.series[1,col-1] 								            	#PTPN11 = FLT3
    
    #adding noise 
    time.series[1,col] <- ifelse(zeta[1,col] < chance.to.flip, 
                                 !time.series[1,col],
                                 time.series[1,col] )
    time.series[2,col] <- ifelse(zeta[2,col] < chance.to.flip, 
                                 !time.series[2,col],
                                 time.series[2,col] )
    time.series[3,col] <- ifelse(zeta[3,col] < chance.to.flip, 
                                 !time.series[3,col],
                                 time.series[3,col])
    time.series[4,col] <- ifelse(zeta[4,col] < chance.to.flip, 
                                 !time.series[4,col],
                                 time.series[,col] )
    time.series[5,col] <- ifelse(zeta[5,col] < chance.to.flip, 
                                 !time.series[5,col],
                                 time.series[5,col] )
    time.series[6,col] <- ifelse(zeta[6,col] < chance.to.flip, 
                                 !time.series[6,col],
                                 time.series[6,col])
    time.series[7,col] <- ifelse(zeta[7,col] < chance.to.flip, 
                                 !time.series[7,col],
                                 time.series[7,col] )
    time.series[8,col] <- ifelse(zeta[8,col] < chance.to.flip, 
                                 !time.series[8,col],
                                 time.series[8,col] )
    time.series[9,col] <- ifelse(zeta[9,col] < chance.to.flip, 
                                 !time.series[9,col],
                                 time.series[9,col])
    time.series[10,col] <- ifelse(zeta[10,col] < chance.to.flip, 
                                  !time.series[10,col],
                                  time.series[10,col])
    time.series[11,col] <- ifelse(zeta[11,col] < chance.to.flip, 
                                  !time.series[11,col],
                                  time.series[11,col])
    time.series[12,col] <- ifelse(zeta[12,col] < chance.to.flip, 
                                  !time.series[12,col],
                                  time.series[12,col])
    time.series[13,col] <- ifelse(zeta[13,col] < chance.to.flip, 
                                  !time.series[13,col],
                                  time.series[13,col])
    time.series[14,col] <- ifelse(zeta[14,col] < chance.to.flip, 
                                  !time.series[14,col],
                                  time.series[14,col])
    time.series[15,col] <- ifelse(zeta[15,col] < chance.to.flip, 
                                  !time.series[15,col],
                                  time.series[15,col])
    time.series[16,col] <- ifelse(zeta[16,col] < chance.to.flip, 
                                  !time.series[16,col],
                                  time.series[16,col])
    time.series[17,col] <- ifelse(zeta[17,col] < chance.to.flip, 
                                  !time.series[17,col],
                                  time.series[17,col])
    time.series[18,col] <- ifelse(zeta[18,col] < chance.to.flip, 
                                  !time.series[18,col],
                                  time.series[18,col])
    time.series[19,col] <- ifelse(zeta[19,col] < chance.to.flip, 
                                  !time.series[19,col],
                                  time.series[19,col])
    time.series[20,col] <- ifelse(zeta[20,col] < chance.to.flip, 
                                  !time.series[20,col],
                                  time.series[20,col])
    time.series[21,col] <- ifelse(zeta[21,col] < chance.to.flip, 
                                  !time.series[21,col],
                                  time.series[21,col])
    time.series[22,col] <- ifelse(zeta[22,col] < chance.to.flip, 
                                  !time.series[22,col],
                                  time.series[22,col])
    time.series[23,col] <- ifelse(zeta[23,col] < chance.to.flip, 
                                  !time.series[23,col],
                                  time.series[23,col])
    time.series[24,col] <- ifelse(zeta[24,col] < chance.to.flip, 
                                  !time.series[24,col],
                                  time.series[24,col])
    
    #fixing the genes that are mutated
    #all info for gene switches are in the dictionary 'genes'
    if (any(grepl("FLT3", mutation_profile))) {time.series[1,col] <- 1}  
    if (any(grepl("CEBPA", mutation_profile))) {time.series[2,col] <- 0}			  
    if (any(grepl("RUNX1", mutation_profile))) {time.series[3,col] <- 1}	  
    if (any(grepl("RAD21", mutation_profile))) {time.series[4,col] <- 0}	  
    if (any(grepl("NRAS", mutation_profile))) {time.series[6,col] <- 1}	  
    if (any(grepl("BRAF", mutation_profile))) {time.series[7,col] <- 1}
    if (any(grepl("NPM1", mutation_profile))) {time.series[9,col] <- 0}	  
    if (any(grepl("MEK", mutation_profile))) {time.series[10,col] <- 0}	  
    if (any(grepl("FOXO", mutation_profile))) {time.series[11,col] <- 1}	  
    if (any(grepl("FBXW7", mutation_profile))) {time.series[12,col] <- 0}		  
    if (any(grepl("MYC", mutation_profile))) {time.series[14,col] <- 1}	  
    if (any(grepl("IDH2", mutation_profile))) {time.series[15,col] <- 1}	  
    if (any(grepl("IDH1", mutation_profile))) {time.series[16,col] <- 1}	  
    if (any(grepl("CDKN2A", mutation_profile))) {time.series[18,col] <- 0}	  
    if (any(grepl("TET2", mutation_profile))) {time.series[19,col] <- 0}	  
    if (any(grepl("MDM2", mutation_profile))) {time.series[20,col] <- 1}	  
    if (any(grepl("WT1", mutation_profile))) {time.series[21,col] <- 1}	  
    if (any(grepl("TP53", mutation_profile))) {time.series[22,col] <- 0}	  
    if (any(grepl("BCL2", mutation_profile))) {time.series[23,col] <- 1}	  
    if (any(grepl("PTPN11", mutation_profile))) {time.series[24,col] <- 1}
    
  }
  
  print("completed time series")
  
  # adjust the time-series to long-format
  time.series <- t(time.series)
  time.series <- data.frame(time.series)
  names(time.series) <- c("FLT3", "CEBPA", "RUNX1", "RAD21", "AKT", "NRAS", "BRAF", "AMPK", "NPM1", 
                          "MEK", "FOXO", "FBXW7", "ERK", "MYC", "IDH2", "IDH1", "OXO2", "CDKN2A", "TET2",
                          "MDM2", "WT1", "TP53", "BCL2", "PTPN11")
  
  #getting the network score and saving it as a vector
  time.series$Proliferation <- with(time.series, ((-1 * FOXO) + (-1 * WT1) + (-1 * CDKN2A) + (-1 * TP53) + MYC))
  time.series$Differentiation <- with(time.series, (CEBPA + RUNX1))
  time.series$Apoptosis <- with(time.series, ((-1 * BCL2) + TP53))
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  print("finished getting the running means")
  
  #final score based on window value of running mean 
  upper_bound <- 10000 
  lower_bound <- 20000    
  
  final_score <- mean(time.series$running_mean[upper_bound:lower_bound], na.rm=TRUE)
  print(paste0("the final score is ", final_score))
  
  return(final_score)    
}

############# using the scorePersonal()  #############
scores <- c()
id <- c()
PB_Blast <- c()
BM_Blast <- c()

for (i in 1:nrow(patient_profile_test_df)){
  patient_id <- as.vector(patient_profile_test_df$labId[i])
  mutation_profile <- as.vector(patient_profile_test_df$symbol[i])
  
  PB_score <- as.numeric(s5_df$'%.Blasts.in.PB'[i])
  BM_score <- as.numeric(s5_df$'%.Blasts.in.BM'[i])
  
  #get the result score
  result <- scorePersonal(mutation_profile)
  
  #save all info into vectors
  scores <- c(scores, result)
  id <- c(id, patient_id)
  PB_Blast <- c(PB_Blast, PB_score)
  BM_Blast <- c(BM_Blast, BM_score)
  
  PB_patient_df <- na.omit(data.frame(id,scores,PB_Blast))
  BM_patient_df <- na.omit(data.frame(id,scores,BM_Blast))
}

#saving as excel sheet
#write_xlsx(
  #PB_patient_df,"path to location/PB_patient.xlsx",
  #col_names = TRUE,
  #format_headers = TRUE)

#making the scatterplots
PB_scatterplot <- ggplot(data = PB_patient_df, aes(x = PB_Blast, y = scores)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Patient PB Blast % vs Score ", x = "% PB Blast", y = "Score") +
  theme_minimal()
PB_scatterplot 

BM_scatterplot <- ggplot(data = BM_patient_df, aes(x = BM_Blast, y = scores)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Patient BM Blast % vs Score ", x = "% BM Blast", y = "Score") +
  theme_minimal()
BM_scatterplot 

#doing a p-test 
PB_pearson <- cor(PB_patient_df$PB_Blast, PB_patient_df$score)
print(PB_pearson)

BM_pearson <- cor(BM_patient_df$BM_Blast, BM_patient_df$score)
print(BM_pearson)