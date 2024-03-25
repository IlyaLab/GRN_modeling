### Created by Tazein Shah
### Last modified on 10/30/2023
### This code runs the FLT3-NPM1-DNMT3A subgroup from the Palma paper, and uses scorePersonal() function to knock in/out genes based on the patient mutation profile

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
  patient_profile_test_df <- patient_profile_df[1:100, ]
#chooses a random 200 samples of the dataframe
  #patient_profile_test_df <- patient_profile_df[sample(nrow(patient_profile_df), 200), ]

#create the dictonary
genes_values <- Dict$new(  
  FLT3 = "time.series[1,col] <- 1" ,
  AKT = "time.series[2,col] <- 1}",	
  CEBPA = "time.series[3,col] <- 0",
  DNMT3A = "time.series[4,col] <- 0",
  GSK3B = "time.series[5,col] <- 0",
  NPM1 = "time.series[6,col] <- 0",
  ARF = "time.series[7,col] <- 0",	   
  HOXA9 = "time.series[8,col] <- 1",
  FBXW7 = "time.series[9,col] <- 0", 
  ERK = "time.series[10,col] <- 1",
  CDKN2A = "time.series[11,col] <- 0",
  STAT5A = "time.series[12,col] <- 1",
  SOX4 = "time.series[13,col] <- 1",
  CCND1 = "time.series[14,col] <- 1",
  MEIS1 = "time.series[15,col] <- 1",
  MYC = "time.series[16,col] <- 1",
  ETV6 = "time.series[17,col] <- 0",
  TP53 = "time.series[18,col] <- 0",
  BCL2 = "time.series[19,col] <- 1",
  .class = "character",
  .overwrite = TRUE)

#genes["FLT3"]= time.series[j,col] <- 1

n.obs <- 65000 # number of observation aka steps 
n.half <- (n.obs/2) #half the steps (for plot making)
p <- 19 # number of variables

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
    time.series[1,col] <-time.series[1,col-1] 									            		#FLT3 = FLT3
    time.series[2,col] <-time.series[1,col-1] 									            		#AKT= FLT3
    time.series[3,col] <-!time.series[1,col-1] 									            		#CEBPA=!FLT3
    time.series[4,col] <-time.series[4,col-1] 									            		#DNMT3A= DNMT3A
    time.series[5,col] <-!time.series[2,col-1] 									           			#GSK3B= !AKT
    time.series[6,col] <-time.series[6,col-1] 									           		 	#NPM1= NPM1
    time.series[7,col] <-time.series[6,col-1] 									            		#ARF= NPM1
    time.series[8,col] <-!time.series[6,col-1] 									           		 	#HOXA9= !NPM1
    time.series[9,col] <-time.series[6,col-1] 								            			#FBXW7= NPM1
    time.series[10,col] <-time.series[1,col-1] 							            				#ERK= FLT3
    time.series[11,col] <-time.series[6,col-1] 														      #CDKN2A= NPM1
    time.series[12,col] <-time.series[1,col-1] 									            		#STAT5A= FLT3
    time.series[13,col] <-!time.series[3,col-1] 									        	  	#SOX4= !CEBPA
    time.series[14,col] <-!(time.series[4,col-1] | time.series[5,col-1])                 		        #CCND1= !(DNMT3A | GSK3B)
    time.series[15,col] <-!(time.series[4,col-1] | ! time.series[8,col-1])   		        		        #MEIS1=!(DNMT3A & !HOXA9)
    time.series[16,col] <-!(time.series[5,col-1] & time.series[9,col-1]) & time.series[10,col-1]    #MYC= !(GSK3B & FBXW7) & ERK
    time.series[17,col] <-!time.series[10,col-1]													      #ETV6=!ERK
    time.series[18,col] <-!time.series[7,col-1] 													      #TP53= ARF
    time.series[19,col] <-time.series[10,col-1] & !time.series[18,col-1]				#BCL2= ERK & !TP53
    
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
    
    #fixing the genes that are mutated
    #all info for gene switches are in the dictionary 'genes'
    if (any(grepl("FLT3", mutation_profile))) {time.series[1,col] <- 1}  
    if (any(grepl("AKT", mutation_profile))) {time.series[2,col] <- 1}	   #look up info		  
    if (any(grepl("CEBPA", mutation_profile))) {time.series[3,col] <- 0}	  
    if (any(grepl("DNMT3A", mutation_profile))) {time.series[4,col] <- 0}	  
    if (any(grepl("GSK3B", mutation_profile))) {time.series[5,col] <- 0}   #look up info
    if (any(grepl("NPM1", mutation_profile))) {time.series[6,col] <- 0}    #look up info
    if (any(grepl("ARF", mutation_profile))) {time.series[7,col] <- 0}	   #look up info
    if (any(grepl("HOXA9", mutation_profile))) {time.series[8,col] <- 1}   #look up info
    if (any(grepl("FBXW7", mutation_profile))) {time.series[9,col] <- 0}   #look up info
    if (any(grepl("ERK", mutation_profile))) {time.series[10,col] <- 1}    #look up info
    if (any(grepl("CDKN2A", mutation_profile))) {time.series[11,col] <- 0} #look up info
    if (any(grepl("STAT5A", mutation_profile))) {time.series[12,col] <- 1}  #look up info
    if (any(grepl("SOX4", mutation_profile))) {time.series[13,col] <- 1}	  #look up info
    if (any(grepl("CCND1", mutation_profile))) {time.series[14,col] <- 1}	  #look up info
    if (any(grepl("MEIS1", mutation_profile))) {time.series[15,col] <- 1}   #look up info
    if (any(grepl("MYC", mutation_profile))) {time.series[16,col] <- 1}	
    if (any(grepl("ETV6", mutation_profile))) {time.series[17,col] <- 0}
    if (any(grepl("TP53", mutation_profile))) {time.series[18,col] <- 0}	  #look up info  
    if (any(grepl("BCL2", mutation_profile))) {time.series[19,col] <- 1}	
    
  }
  
  print("completed time series")
  
  # adjust the time-series to long-format
  time.series <- t(time.series)
  time.series <- data.frame(time.series)
  names(time.series) <- c("FLT3", "AKT", "CEBPA", "DNMT3A", "GSK3B", "NPM1", "ARF", "HOXA9", "FBXW7", "ERK", 
                          "CDKN2A", "STAT5A", "SOX4", "CCND1", "MEIS1", "MYC", "ETV6",
                          "TP53", "BCL2")
  
  #getting the network score and saving it as a vector
  time.series$Proliferation <- with(time.series, (MYC + CCND1 + SOX4 + MEIS1 + STAT5A))
  time.series$Differentiation <- with(time.series, (CEBPA + ETV6 + (-1 *MEIS1)))
  time.series$Apoptosis <- with(time.series, ((-1 * BCL2) + TP53))
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  print("finished getting the running means")
  
  #final score based on window value of running mean 
  upper_bound <- 60000  
  lower_bound <- 64000    
  
  final_score <- mean(time.series$running_mean[upper_bound:lower_bound], na.rm=TRUE)
  print(paste0("the final score is ", final_score))
  
  return(final_score)    
}

############# for finding the convergence #############

#rolling median of the means 
#rolling_median <- rollmedian(time.series$running_mean, k = window)

#diagrams
#hist(time.series[window:n.obs,19],
     #xlim = c(-8,2),
     #ylim = c(0, 20000),
     #main = "Histogram of running mean",
     #xlab = "running mean"
#) # plot the distribution of the running mean 

#hist(rolling_median,
  #   xlim = c(-5,-2),
  #   ylim = c(0, 20000),
  #   main = "Histogram of rolling median",
  #   xlab = "rolling median"
#) # plot the distribution of the rolling median 

#plot(seq(1,n.obs), time.series$running_mean[upper_bound:lower_bound], xlab = "time points", ylab = "running mean")

#getting the relative error points to check if the score converges
#s1 <- mean(time.series$Score[500:1000])  #for 1000
#s2 <- mean(time.series$Score[500:2000])  #for 2000
#s3 <- mean(time.series$Score[500:4000])  #for 4000
#s4 <- mean(time.series$Score[500:8000])  #for 8000
#s5 <- mean(time.series$Score[500:16000]) #for 16000
#s6 <- mean(time.series$Score[500:32000]) #for 32000
#s7 <- mean(time.series$Score[500:64000]) #for 64000

#getting the relative error points
#point_1 <- ((s2-s1)/s1)
#point_2 <- ((s3-s2)/s2)
#point_3 <- ((s4-s3)/s3)
#point_4 <- ((s5-s4)/s4)
#point_5 <- ((s6-s5)/s5)
#point_6 <- ((s7-s6)/s6)

#saving into vectors
#relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6)
#time <- c(1000, 2000, 4000, 8000, 16000,32000)

#create a new dataframe for the relative_error
#mean_time_points <- data.frame(time, relative_error)

#mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = time, y = relative_error)) +
 # geom_point() +
 # labs(title = "Relative Error Scatterplot (FLT3-NPM1-DNMT3A)", x = "time interval", y = "relative error") +
 # theme_minimal()
#mean_time_scatterplot + coord_cartesian(xlim = c(800,65000), ylim = c(-2, 1))

#final score based on window value of running mean 
#upper_bound <- 60000 
#lower_bound <- 64000    
#the last 500 values (without NA) are considered in the mean 

#final_score <- mean(time.series$running_mean[upper_bound:lower_bound], na.rm=TRUE)
#print(final_score)

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
  labs(title = "Patient PB Blast % vs Score (FLT3-NPM1-DNMT3A)", x = "% PB Blast", y = "Score") +
  theme_minimal()
PB_scatterplot 

BM_scatterplot <- ggplot(data = BM_patient_df, aes(x = BM_Blast, y = scores)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Patient BM Blast % vs Score (FLT3-NPM1-DNMT3A)", x = "% BM Blast", y = "Score") +
  theme_minimal()
BM_scatterplot 

#doing a p-test 
PB_pearson <- cor(PB_patient_df$PB_Blast, PB_patient_df$score)
print(PB_pearson)

BM_pearson <- cor(BM_patient_df$BM_Blast, BM_patient_df$score)
print(BM_pearson)