### Created by Tazein Shah (uses FLT3_NPM1_DNMT3A_model1.xlsx from Guangrong Qin)
### Last modified on 11/13/2023
### Using the FLT3_NPM1_DNMT3A_model1 from Guangrong, this code simulates the time series, calculates apop, prolif, diff, and network scores. Then it gets the final score + scatteplots


#load needed packages
library(BoolNet)
library(reshape2)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggrepel)
library(writexl)
library(readxl)
library(dplyr)

#importing the excel sheets 
s7_df <- read_excel("/path/raw_data/s7_table.xlsx")
s5_df <- read_excel("/path/raw_data/s5_table.xlsx")

patient_profile_df <- s7_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#patient_profile_test_df <- patient_profile_df[1:5, ] #first 5 samples for the test
patient_profile_test_df <- patient_profile_df #doing the whole dataframe

#how the mutation profile will be formatted
mutation_profile <- as.vector(patient_profile_test_df$symbol[2])

n.obs <- 80000 # number of observation aka steps 
p <- 25 # number of variables

#this is a threshold to flip the node 
chance.to.flip <- 0.01

############# scorePersonal() function #############

scorePersonal <- function(mutation_profile) { 
  
  #separate the values by the commas
  mutation_profile <- (unlist(strsplit(mutation_profile,",")))
  
  #remove repeated values
  mutation_profile <- unique(mutation_profile)
  
  #so now it should look like this (example)
  #List values: 
  #"FLT3"
  #"DNMT3A"
  #"NPM1"
  
  print(paste0(mutation_profile))
  print("preparing to create time series")
  
  # 22 time-series of noise for 22 variables
  zeta <- matrix(runif(n.obs * p, 0, 1), ncol = p)
  
  zeta <- t(zeta) # adjust zeta to a wide-format
  
  time.series <- matrix(rep(0, p * n.obs), nrow =p , ncol =  n.obs)
  
  time.series[,1] <- runif(p,0, 1) # initial value of the three variables were randomly chosen
  
  time.series[,1] <-ifelse(time.series[,1] >= .5, TRUE,FALSE) # if the randomly generated number is greater than .5, then we convert the value to TRUE, otherwise to FALSE
  
  for (col in 2:n.obs){
    # simulate each time step based on the previous time step and Boolean functions
    time.series[1,col] = !(time.series[2,col-1] | time.series[8,col-1]) & (time.series[18,col-1] | time.series[2,col-1] | time.series[11,col-1]) #AKT1, !(BRAF|GSK3B) & (BCL2|BRAF|MAPK1)
    time.series[2,col] = time.series[10,col-1] #BRAF, MAP2K2
    time.series[3,col] = !(time.series[13,col-1]) & time.series[17,col-1] #CDKN2A, !(NPM1) & (TP53)
    time.series[4,col] = !(time.series[12,col-1] | time.series[19,col-1]) #CEBPA, !(MYC|SOX4)
    time.series[5,col] = !(time.series[20,col-1] | time.series[3,col-1] | time.series[9,col-1] | time.series[21,col-1]) #DNMT3A, !(CCND1|CDKN2A|HOXA9|MEIS1)
    time.series[6,col] = !(time.series[12,col-1]) #FBXW7, !MYC
    time.series[7,col] = !(time.series[4,col-1]) & (time.series[1,col-1] | time.series[12,col-1] | time.series[15,col-1] | time.series[16,col-1]) #FLT3, !(CEBPA) & (AKT1|MYC|PTPN11|STAT5A)
    time.series[8,col] = !(time.series[12,col-1] | time.series[20,col-1]) & (time.series[4,col-1] | time.series[17,col-1])	#GSK3B, !(MYC|CCND1) & (CEBPA|TP53)
    time.series[9,col] = time.series[21,col-1] #HOXA9, MEIS1
    time.series[10,col] = time.series[11,col] #MAP2K2, MAPK1
    time.series[11,col] = !(time.series[2,col-1] | time.series[4,col-1] | time.series[22,col-1] |time.series[8,col-1]) & (time.series[18,col-1] | time.series[12,col-1] | time.series[16,col-1]) #MAPK1, !(BRAF|CEBPA|ETV6|GSK3B) & (BCL2|MYC|STAT5A)
    time.series[12,col] = !(time.series[3,col-1]) & (time.series[20,col-1] | time.series[5,col-1]) #MYC, !(CDKN2A) & (CCND1|DNMT3A)
    time.series[13,col] = !(time.series[9,col-1]) & (time.series[3,col-1] | time.series[6,col-1]) #NPM1, !(HOXA9) & (CDKN2A|FBXW7)
    time.series[14,col] = time.series[2,col-1] #NRAS, BRAF
    time.series[15,col] = time.series[14,col-1] #PTPN11, NRAS
    time.series[16,col] = !(time.series[5,col-1]) & time.series[5,col-1] #STAT5A, !(DNMT3A) & (DNMT3A)	
    time.series[17,col] = !(time.series[18,col-1]) #TP53, !(BCL2) 
    time.series[18,col] = time.series[18,col-1] #BCL2, BCL2
    time.series[19,col] = time.series[19,col-1] #SOX4, SOX4
    time.series[20,col] = time.series[20,col-1] #CCND1, CCND1
    time.series[21,col] = time.series[21,col-1]  #MEIS1, MEIS1
    time.series[22,col] = time.series[22,col-1] #ETV6, ETV6
    time.series[23,col] =  !(time.series[3,col-1]) & (time.series[16,col-1] | time.series[12,col-1]| time.series[21,col-1] | time.series[20,col-1]) #Proliferation
    time.series[24,col] =  !(time.series[19,col-1]) & (time.series[4,col-1] |  time.series[22,col-1]) #Differentiation
    time.series[25,col] =  !(time.series[18,col-1]) & time.series[17,col-1] #Apoptosis

    
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
                                  time.series[20,col] )
    time.series[21,col] <- ifelse(zeta[21,col] < chance.to.flip, 
                                  !time.series[21,col],
                                  time.series[21,col] )
    time.series[22,col] <- ifelse(zeta[22,col] < chance.to.flip, 
                                  !time.series[22,col],
                                  time.series[22,col])
    time.series[23,col] <- ifelse(zeta[23,col] < chance.to.flip, 
                                  !time.series[23,col],
                                  time.series[23,col])
    time.series[24,col] <- ifelse(zeta[24,col] < chance.to.flip, 
                                  !time.series[24,col],
                                  time.series[24,col])
    time.series[25,col] <- ifelse(zeta[25,col] < chance.to.flip, 
                                  !time.series[25,col],
                                  time.series[25,col])
      }
  
  # adjust the time-series to long-format
  time.series <- t(time.series)
  time.series <- data.frame(time.series)
  names(time.series) <- c("AKT1","BRAF", "CDKN2A", "CEBPA", "DNMT3A", "FBXW7", "FLT3", "GSK3B", "HOXA9", "MAP2K2", "MAPK1",
                          "MYC", "NPM1", "NRAS", "PTPN11", "STAT5A", "TP53", "BCL2", "SOX4", "CCND1", 
                          "MEIS1", "ETV6", "Proliferation", "Differentiation", "Apoptosis")
  
  
  #fixing the genes that are mutated
  #knock out whole columns based on mutation profile
  for (i in mutation_profile){
    switch(i, 
           "AKT1"={time.series$AKT1 <- 1},
           "BRAF"={time.series$BRAF <- 1},  
           "CDKN2A"={time.series$CDKN2A <- 0},  
           "CEBPA"={time.series$CEBPA <- 0},  
           "DNMT3A"={time.series$DNMT3A <- 0},  
           "FBXW7"={time.series$FBXW7 <- 0},  
           "FLT3"={time.series$FLT3 <- 1},  
           "GSK3B"={time.series$GSK3B <- 1}, 
           "HOXA9"={time.series$HOXA9 <- 1},  
           "MAP2K2"={time.series$MAP2K2 <- 1},
           "MAPK1"={time.series$MAPK1 <- 1},  
           "MYC"={time.series$MYC <- 1},  
           "NPM1"={time.series$NPM1 <- 1},  
           "NRAS"={time.series$NRAS <- 1},  
           "PTPN11"={time.series$PTPN11 <- 1},  
           "STAT5A"={time.series$STAT5A <- 1},  
           "TP53"={time.series$TP53 <- 0},  
           "BCL2"={time.series$BCL2 <- 1},  
           "SOX4"={time.series$SOX4 <- 1},  
           "CCND1"={time.series$CCND1 <- 1},  
           "MEIS1"={time.series$MEIS1 <- 1},  
           "ETV6"={time.series$ETV6 <- 0},     
    )}
  
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  upper_bound <- n.obs  
  lower_bound <- upper_bound-5000 #I added more values so I can get a more accurate score    
  
  #saving into the matrix 
  prolif <- mean(time.series$Proliferation[lower_bound:upper_bound], na.rm=TRUE)
  diff <- mean(time.series$Differentiation[lower_bound:upper_bound], na.rm=TRUE)
  apop <- mean(time.series$Apoptosis[lower_bound:upper_bound], na.rm=TRUE)
  final_score <- mean(time.series$running_mean[lower_bound:upper_bound], na.rm=TRUE)
  
  print("scores have been calculated")
  
  score_vector <- c(prolif,diff,apop,final_score)
  
  return(score_vector)    
}

############# using the scorePersonal()  #############
scores = matrix(, nrow = nrow(patient_profile_test_df), ncol = 8)
counter <- 0

for (i in 1:nrow(patient_profile_test_df)){
  mutation_profile <- as.vector(patient_profile_test_df$symbol[i])
  
  scores[i,1] <- as.vector(patient_profile_test_df$labId[i]) 
  scores[i,2] <- as.vector(patient_profile_test_df$symbol[i])
  scores[i,3] <- as.numeric(s5_df$'%.Blasts.in.PB'[i])
  scores[i,4] <- as.numeric(s5_df$'%.Blasts.in.BM'[i])
  
  #get the result score
  result <- scorePersonal(mutation_profile)
  
  scores[i,5] <- result[1] #proliferation
  scores[i,6] <- result[2] #differentiation
  scores[i,7] <- result[3] #apoptosis
  scores[i,8] <- result[4] #final_score
  
  
  counter <- counter + 1
  
  print(paste0("patient ", counter))
}

#saving the scores matrix as a dataframe and adding colnames 
scores_df <- as.data.frame(scores)
#specify column names
colnames(scores_df) <- c('LabId', 'Symbol', 'PB_Blast', 'BM_Blast', 
                         'Proliferation', 'Differentiation', 'Apoptosis',
                         'Final_Score')

scores_df <- scores_df %>% mutate_at(c('Final_Score', 'Apoptosis','Differentiation','Proliferation','BM_Blast','PB_Blast'), as.numeric)

#saving as excel sheet
#write_xlsx(
#scores_df,"/path/raw_data/dnmt3a_model1.xlsx",
#col_names = TRUE,
#format_headers = TRUE)

############# making the scatterplots  #############

#load the file from excel sheet
dnmt3a_model1 <- read_excel("/path/raw_data/dnmt3a_model1.xlsx")
dnmt3a_model1 <- na.omit(dnmt3a_model1)

### proliferation ###

prolif_PB_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = PB_Blast, y = Proliferation)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "PB Blast vs Proliferation (DNMT3A Model 1)", x = "% PB Blast", y = "Proliferation Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
prolif_PB_scatterplot 

prolif_BM_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = BM_Blast, y = Proliferation)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "BM Blast vs Proliferation (DNMT3A Model 1)", x = "% BM Blast", y = "Proliferation Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
prolif_BM_scatterplot 

print((cor(dnmt3a_model1$PB_Blast, dnmt3a_model1$Proliferation)))
print((cor(dnmt3a_model1$BM_Blast, dnmt3a_model1$Proliferation)))

### differentiation ###

diff_PB_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = PB_Blast, y = Differentiation)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "PB Blast vs Differentiation (DNMT3A Model 1)", x = "% PB Blast", y = "Differentiation Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
diff_PB_scatterplot 

diff_BM_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = BM_Blast, y = Differentiation)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "BM Blast vs Differentiation (DNMT3A Model 1)", x = "% BM Blast", y = "Differentiation Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
diff_BM_scatterplot 

print((cor(dnmt3a_model1$PB_Blast, dnmt3a_model1$Differentiation)))
print((cor(dnmt3a_model1$BM_Blast, dnmt3a_model1$Differentiation)))

### apoptosis ###

apop_PB_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = PB_Blast, y = Apoptosis)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "PB Blast vs Apoptosis (DNMT3A Model 1)", x = "% PB Blast", y = "Apoptosis Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
apop_PB_scatterplot 

apop_BM_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = BM_Blast, y = Apoptosis)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "BM Blast vs Apoptosis (DNMT3A Model 1)", x = "% BM Blast", y = "Apoptosis Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
apop_BM_scatterplot 

print((cor(dnmt3a_model1$PB_Blast, dnmt3a_model1$Apoptosis)))
print((cor(dnmt3a_model1$BM_Blast, dnmt3a_model1$Apoptosis)))

### final score ###

final_PB_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = PB_Blast, y = Final_Score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "PB Blast vs Final Score (DNMT3A Model 1)", x = "% PB Blast", y = "Final Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
final_PB_scatterplot 

final_BM_scatterplot <- ggplot(data = dnmt3a_model1, aes(x = BM_Blast, y = Final_Score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "BM Blast vs Final Score (DNMT3A Model 1)", x = "% BM Blast", y = "Final Score") +
  #scale_x_discrete(breaks = seq(0, 100, by = 20)) + 
  #scale_y_discrete(breaks = seq(0, 100, by = 20)) + 
  theme_minimal()
final_BM_scatterplot 

print((cor(dnmt3a_model1$PB_Blast, dnmt3a_model1$Final_Score)))
print((cor(dnmt3a_model1$BM_Blast, dnmt3a_model1$Final_Score)))
