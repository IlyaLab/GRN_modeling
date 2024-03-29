### Created by Tazein Shah
### Last modified on 10/4/2023
### This code runs the whole AML network from fig 3 in the Palma paper
### Uses a function called scorePersonal() to knock in/out genes based on patient mutation profiles

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

#genes["FLT3"]= time.series[j,col] <- 1

n.obs <- 35000 # number of observation aka steps 
n.half <- (n.obs/2) #half the steps (for plot making)
p <- 78 # number of variables

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
    time.series[1,col]<- time.series[1,col-1]																							#KITLG = KITLG
    time.series[2,col]<- time.series[2,col-1]																							#FLT3LG = FLT3LG
    time.series[3,col]<- time.series[6,col-1] & time.series[11,col-1] & time.series[10,col-1]											#GRB2 = BCR_ABL | FLT3 | KIT
    time.series[4,col]<- time.series[4,col-1]																							#KDM5A = KDM5A
    time.series[5,col]<- !time.series[4,col-1] & time.series[5,col-1] 										#PTEN = !KDM5A & PTEN
    time.series[6,col]<- !time.series[5,col-1]																						#BCR_ABL = !PTEN
    time.series[7,col]<- time.series[7,col-1]																							#CBLB = CBLB
    time.series[8,col]<- time.series[3,col-1]																							#CBL = GRB2
    time.series[9,col]<- time.series[3,col-1]																							#SOS1 = GRB2
    time.series[10,col]<- !(time.series[8,col-1] | time.series[7,col-1]) & (time.series[16,col-1] | time.series[1,col-1])				#KIT = !(CBL | CBLB) & (AML1_ETO | KITLG)
    time.series[11,col]<- !(time.series[8,col-1] | time.series[7,col-1]) &	time.series[2,col-1]										#FLT3 = !(CBL | CBLB ) & FLT3LG
    time.series[12,col]<- time.series[12,col-1]																						#CDK1 = CDK1
    time.series[13,col]<- time.series[13,col-1]																						#STAG2 = STAG2 
    time.series[14,col]<- time.series[14,col-1]																						#SH2B3 = SH2B3
    time.series[15,col]<- !time.series[11,col-1]																					#PTPN6 =!FLT3
    time.series[16,col]<- time.series[16,col-1]																						#AML1_ETO = AML1_ETO
    time.series[17,col]<- time.series[17,col-1]																						#PML_RARalpha = PML_RARalpha
    time.series[18,col]<- time.series[11,col-1] & time.series[10,col-1] 									#PTPN11 = FLT3 | KIT
    time.series[19,col]<- !time.series[12,col-1]																					#CUX1 =!CDK1
    time.series[20,col]<- time.series[13,col-1]																						#RAD21 = STAG2
    time.series[21,col]<- !time.series[14,col-1] & time.series[29,col-1]									#BCL2L1 = !SH2B3 & STAT5A
    time.series[22,col]<- !(time.series[16,col-1] | time.series[15,col-1] | time.series[14,col-1]) & (time.series[6,col-1] | time.series[10,col-1]) 	#JAK2 = !(AML1_ETO | PTPN6 | SH2B3) & (BCR_ABL | KIT)
    time.series[23,col]<- time.series[17,col-1]																						#CCNA1 = PML_RARalpha
    time.series[24,col]<- !(time.series[16,col-1] | time.series[11,col-1]) | (time.series[17,col-1]) 									#CEBPA =!(AML1_ETO | FLT3 | PML_RARalpha)
    time.series[25,col]<- !time.series[19,col-1]																					#PIK3IP1 =!CUX1
    time.series[26,col]<- time.series[9,col-1] & time.series[18,col-1]										#NRAS = SOS1 | PTPN11
    time.series[27,col]<- !time.series[20,col-1] & time.series[17,col-1] 									#GATA2 = !RAD21 & PML_RARalpha
    time.series[28,col]<- !time.series[21,col-1] 																					#BAX = !BCL2L1
    time.series[29,col]<- time.series[22,col-1] & time.series[6,col-1]										#STAT5A = JAK2 | BCR_ABL 
    time.series[30,col]<- !time.series[24,col-1] 																					#SOX4 = !CEBPA
    time.series[31,col]<- !(time.series[25,col-1]) & (time.series[11,col-1] | time.series[10,col-1] | time.series[26,col-1] | time.series[32,col-1])	#PI3K = !PIK3IP1 & (FLT3 | KIT | NRAS | SPI1)
    time.series[32,col]<- !time.series[27,col-1] & time.series[24,col-1]									#SPI1 = !GATA2 & CEBPA
    time.series[33,col]<- time.series[29,col-1] 																					#PIM = STAT5A
    time.series[34,col]<- time.series[30,col-1] & time.series[11,col-1]										#CTNNB1 =  SOX4 | FLT3 
    time.series[35,col]<- time.series[31,col-1] 																					#AKT = PI3K 
    time.series[36,col]<- !time.series[32,col-1] & time.series[27,col-1]									#GATA1 = !SPI1 & GATA2
    time.series[37,col]<- time.series[37,col-1] 																					#NUP98_Fusion = NUP98_Fusion
    time.series[38,col]<- time.series[38,col-1] 																					#MLL_Fusion = MLL_Fusion
    time.series[39,col]<- !time.series[35,col-1]																					#FOXO = !AKT
    time.series[40,col]<- time.series[35,col-1] 																					#EP300 = AKT
    time.series[41,col]<- !time.series[35,col-1] & time.series[26,col-1]									#BRAF = !AKT & NRAS
    time.series[42,col]<- time.series[36,col-1]  																					#CBFB = GATA1
    time.series[43,col]<- time.series[37,col-1]																						#CDK6 = NUP98_Fusion
    time.series[44,col]<- time.series[38,col-1]																						#MECOM = MLL_Fusion
    time.series[45,col]<- time.series[38,col-1]																						#CBFbeta_MYH11 = MLL_Fusion
    time.series[46,col]<- time.series[46,col-1]																						#NPM1 = NPM1
    time.series[47,col]<- time.series[39,col-1]  																					#IDH1 = FOXO
    time.series[48,col]<- time.series[48,col-1]																						#IDH2 = IDH2
    time.series[49,col]<- time.series[41,col-1] 																					#ERK1/2 = BRAF 
    time.series[50,col]<- !(time.series[38,col-1] | time.series[45,col-1] | time.series[20,col-1]) & (time.series[11,col-1] | time.series[43,col-1] | time.series[42,col-1] | time.series[50,col-1])	#RUNX1 = !(MLL_Fusion | CBFbeta_MYH11 | RAD21) & (FLT3 | CDK6 | CBFB | RUNX1)
    time.series[51,col]<- time.series[45,col-1]																						#DOT1L = CBFbeta_MYH11
    time.series[52,col]<- time.series[46,col-1]																						#FBXW7 = NPM1
    time.series[53,col]<- time.series[53,col-1]																						#AMPK = AMPK 
    time.series[54,col]<- time.series[54,col-1]																						#ASXL2 = ASXL2
    time.series[55,col]<- time.series[47,col-1] & time.series[48,col-1]										#OXO2 = IDH1 | IDH2
    time.series[56,col]<- !time.series[49,col-1] 																					#ETV6 = !ERK1/2
    time.series[57,col]<- time.series[49,col-1] 																					#AP1 = ERK1/2 
    time.series[58,col]<- time.series[58,col-1]																						#DNMT3A = DNMT3A
    time.series[59,col]<- !time.series[35,col-1]																					#AEZH2 = !AKT
    time.series[60,col]<- time.series[60,col-1]																						#ASXL1 = ASXL1
    time.series[61,col]<- !(time.series[52,col-1] | time.series[50,col-1]) & (time.series[49,col-1] | time.series[51,col-1])			#MYC = !(FBXW7 | RUNX1) & (ERK1/2 | DOT1L)
    time.series[62,col]<- !time.series[53,col-1] & !time.series[35,col-1] 								#MTOR = !(AMPK | AKT)
    time.series[63,col]<- time.series[55,col-1] | time.series[53,col-1] | time.series[54,col-1]											#TET2 = OXO2 | AMPK | ASXL2
    time.series[64,col]<- time.series[64,col-1] 																					#PHF6 = PHF6
    time.series[65,col]<- !time.series[58,col-1] 																					#CCND1 = !DNMT3A 
    time.series[66,col]<- !(time.series[59,col-1] | time.series[58,col-1] | time.series[60,col-1]) & time.series[51,col-1]				#HOXA9 = !(AEZH2 | DNMT3A | ASXL1) & DOT1L 
    time.series[67,col]<- !(time.series[58,col-1] | time.series[61,col-1] | time.series[16,col-1]) & (time.series[60,col-1]	| time.series[46,col-1])	#CDKN2A = !( DNMT3A | MYC | AML1_ETO) & (ASXL1 | NPM1)
    time.series[68,col]<- time.series[68,col-1] 																					#SRSF2 = SRSF2
    time.series[69,col]<- time.series[63,col-1]        																		#WT1 = TET2
    time.series[70,col]<- time.series[70,col-1]																						#BCOR = BCOR
    time.series[71,col]<- !time.series[64,col-1]																					#UBTF = !PHF6
    time.series[72,col]<- !time.series[58,col-1] & (time.series[66,col-1] | time.series[51,col-1])										#MEIS1 =!DNMT3A & (HOXA9 | DOT1L)
    time.series[73,col]<- !time.series[67,col-1] & (time.series[68,col-1] | time.series[77,col-1] | time.series[35,col-1])	 			#MDM2 = !CDKN2A & (SRSF2 | TP53 | AKT)
    time.series[74,col]<- time.series[74,col-1]																						#U2AF1 = U2AF1
    time.series[75,col]<- time.series[75,col-1]																						#XPO1 = XPO1
    time.series[76,col]<- time.series[76,col-1]																						#CREBBP = CREBBP
    time.series[77,col]<- !(time.series[73,col-1] | time.series[45,col-1] | time.series[75,col-1]) & (time.series[74,col-1]) & time.series[76,col-1]	#TP53 = !(MDM2 | CBFbeta_MYH11 | XPO1 | U2AF1) & (CREBBP)
    time.series[78,col]<- !time.series[77,col-1] & (time.series[49,col-1] | time.series[51,col-1])										#BCL2 = !TP53 & (ERK1/2 | DOT1L)    
    
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
                                  time.series[,col] )								 
    time.series[24,col] <- ifelse(zeta[24,col] < chance.to.flip, 
                                  !time.series[24,col],
                                  time.series[24,col] )
    time.series[25,col] <- ifelse(zeta[25,col] < chance.to.flip, 
                                  !time.series[25,col],
                                  time.series[25,col])
    time.series[26,col] <- ifelse(zeta[26,col] < chance.to.flip, 
                                  !time.series[26,col],
                                  time.series[26,col] )
    time.series[27,col] <- ifelse(zeta[27,col] < chance.to.flip, 
                                  !time.series[27,col],
                                  time.series[27,col] )
    time.series[28,col] <- ifelse(zeta[28,col] < chance.to.flip, 
                                  !time.series[28,col],
                                  time.series[28,col])
    time.series[29,col] <- ifelse(zeta[29,col] < chance.to.flip, 
                                  !time.series[29,col],
                                  time.series[29,col])
    time.series[30,col] <- ifelse(zeta[30,col] < chance.to.flip, 
                                  !time.series[30,col],
                                  time.series[30,col])							  
    time.series[31,col] <- ifelse(zeta[31,col] < chance.to.flip, 
                                  !time.series[31,col],
                                  time.series[31,col])
    time.series[32,col] <- ifelse(zeta[32,col] < chance.to.flip, 
                                  !time.series[32,col],
                                  time.series[32,col])
    time.series[33,col] <- ifelse(zeta[33,col] < chance.to.flip, 
                                  !time.series[33,col],
                                  time.series[33,col])								  
    time.series[34,col] <- ifelse(zeta[34,col] < chance.to.flip, 
                                  !time.series[34,col],
                                  time.series[34,col])								  
    time.series[35,col] <- ifelse(zeta[35,col] < chance.to.flip, 
                                  !time.series[35,col],
                                  time.series[35,col])								  
    time.series[36,col] <- ifelse(zeta[36,col] < chance.to.flip, 
                                  !time.series[36,col],
                                  time.series[36,col])								  
    time.series[37,col] <- ifelse(zeta[37,col] < chance.to.flip, 
                                  !time.series[37,col],
                                  time.series[37,col])
    time.series[38,col] <- ifelse(zeta[38,col] < chance.to.flip, 
                                  !time.series[38,col],
                                  time.series[38,col])
    time.series[39,col] <- ifelse(zeta[39,col] < chance.to.flip, 
                                  !time.series[39,col],
                                  time.series[39,col] )
    time.series[40,col] <- ifelse(zeta[40,col] < chance.to.flip, 
                                  !time.series[40,col],
                                  time.series[40,col] )
    time.series[41,col] <- ifelse(zeta[41,col] < chance.to.flip, 
                                  !time.series[41,col],
                                  time.series[41,col])
    time.series[42,col] <- ifelse(zeta[42,col] < chance.to.flip, 
                                  !time.series[42,col],
                                  time.series[42,col] )									 
    time.series[43,col] <- ifelse(zeta[43,col] < chance.to.flip, 
                                  !time.series[43,col],
                                  time.series[43,col] )							 
    time.series[44,col] <- ifelse(zeta[44,col] < chance.to.flip, 
                                  !time.series[44,col],
                                  time.series[44,col])
    time.series[45,col] <- ifelse(zeta[45,col] < chance.to.flip, 
                                  !time.series[45,col],
                                  time.series[45,col] )
    time.series[46,col] <- ifelse(zeta[46,col] < chance.to.flip, 
                                  !time.series[46,col],
                                  time.series[46,col] )
    time.series[47,col] <- ifelse(zeta[47,col] < chance.to.flip, 
                                  !time.series[47,col],
                                  time.series[47,col])
    time.series[48,col] <- ifelse(zeta[48,col] < chance.to.flip, 
                                  !time.series[48,col],
                                  time.series[48,col])
    time.series[49,col] <- ifelse(zeta[49,col] < chance.to.flip, 
                                  !time.series[49,col],
                                  time.series[49,col])
    time.series[50,col] <- ifelse(zeta[50,col] < chance.to.flip, 
                                  !time.series[50,col],
                                  time.series[50,col])
    time.series[51,col] <- ifelse(zeta[51,col] < chance.to.flip, 
                                  !time.series[51,col],
                                  time.series[51,col])
    time.series[52,col] <- ifelse(zeta[52,col] < chance.to.flip, 
                                  !time.series[52,col],
                                  time.series[52,col])
    time.series[53,col] <- ifelse(zeta[53,col] < chance.to.flip, 
                                  !time.series[53,col],
                                  time.series[53,col])
    time.series[54,col] <- ifelse(zeta[54,col] < chance.to.flip, 
                                  !time.series[54,col],
                                  time.series[54,col])
    time.series[55,col] <- ifelse(zeta[55,col] < chance.to.flip, 
                                  !time.series[55,col],
                                  time.series[55,col])
    time.series[56,col] <- ifelse(zeta[56,col] < chance.to.flip, 
                                  !time.series[56,col],
                                  time.series[56,col])
    time.series[57,col] <- ifelse(zeta[57,col] < chance.to.flip, 
                                  !time.series[57,col],
                                  time.series[57,col])
    time.series[58,col] <- ifelse(zeta[58,col] < chance.to.flip, 
                                  !time.series[58,col],
                                  time.series[58,col] )
    time.series[59,col] <- ifelse(zeta[59,col] < chance.to.flip, 
                                  !time.series[59,col],
                                  time.series[59,col] )
    time.series[60,col] <- ifelse(zeta[60,col] < chance.to.flip, 
                                  !time.series[60,col],
                                  time.series[60,col])
    time.series[61,col] <- ifelse(zeta[61,col] < chance.to.flip, 
                                  !time.series[61,col],
                                  time.series[,col] )
    time.series[62,col] <- ifelse(zeta[62,col] < chance.to.flip, 
                                  !time.series[62,col],
                                  time.series[62,col] )
    time.series[63,col] <- ifelse(zeta[63,col] < chance.to.flip, 
                                  !time.series[63,col],
                                  time.series[63,col])
    time.series[64,col] <- ifelse(zeta[64,col] < chance.to.flip, 
                                  !time.series[64,col],
                                  time.series[64,col] )
    time.series[65,col] <- ifelse(zeta[65,col] < chance.to.flip, 
                                  !time.series[65,col],
                                  time.series[65,col] )
    time.series[66,col] <- ifelse(zeta[66,col] < chance.to.flip, 
                                  !time.series[66,col],
                                  time.series[66,col])
    time.series[67,col] <- ifelse(zeta[67,col] < chance.to.flip, 
                                  !time.series[67,col],
                                  time.series[67,col])
    time.series[68,col] <- ifelse(zeta[68,col] < chance.to.flip, 
                                  !time.series[68,col],
                                  time.series[68,col])
    time.series[69,col] <- ifelse(zeta[69,col] < chance.to.flip, 
                                  !time.series[69,col],
                                  time.series[69,col])
    time.series[70,col] <- ifelse(zeta[70,col] < chance.to.flip, 
                                  !time.series[70,col],
                                  time.series[70,col])
    time.series[71,col] <- ifelse(zeta[71,col] < chance.to.flip, 
                                  !time.series[71,col],
                                  time.series[71,col])
    time.series[72,col] <- ifelse(zeta[72,col] < chance.to.flip, 
                                  !time.series[72,col],
                                  time.series[72,col])
    time.series[73,col] <- ifelse(zeta[73,col] < chance.to.flip, 
                                  !time.series[73,col],
                                  time.series[73,col])
    time.series[74,col] <- ifelse(zeta[74,col] < chance.to.flip, 
                                  !time.series[74,col],
                                  time.series[74,col])
    time.series[75,col] <- ifelse(zeta[75,col] < chance.to.flip, 
                                  !time.series[75,col],
                                  time.series[75,col])
    time.series[76,col] <- ifelse(zeta[76,col] < chance.to.flip, 
                                  !time.series[76,col],
                                  time.series[76,col])
    time.series[77,col] <- ifelse(zeta[77,col] < chance.to.flip, 
                                  !time.series[77,col],
                                  time.series[77,col] )
    time.series[78,col] <- ifelse(zeta[78,col] < chance.to.flip, 
                                  !time.series[78,col],
                                  time.series[78,col] )
    
    #fixing the genes that are mutated
    #all info for gene switches are in the dictionary 'genes'
    if (any(grepl("KITLG", mutation_profile))) {time.series[1,col] <- 1}  
    if (any(grepl("FLT3LG", mutation_profile))) {time.series[2,col] <- 1}	   		  
    if (any(grepl("GRB2", mutation_profile))) {time.series[3,col] <- 1}	  
    if (any(grepl("KDM5A", mutation_profile))) {time.series[4,col] <- 1}	  
    if (any(grepl("PTEN", mutation_profile))) {time.series[5,col] <- 0}  
    if (any(grepl("BCR_ABL", mutation_profile))) {time.series[6,col] <- 1}    
    if (any(grepl("CBLB", mutation_profile))) {time.series[7,col] <- 0}   
    if (any(grepl("CBL", mutation_profile))) {time.series[8,col] <- 0}  
    if (any(grepl("SOS1", mutation_profile))) {time.series[9,col] <- 1}
    if (any(grepl("KIT", mutation_profile))) {time.series[10,col] <- 1}		
    if (any(grepl("FLT3", mutation_profile))) {time.series[11,col] <- 1}	  
    if (any(grepl("CDK1", mutation_profile))) {time.series[12,col] <- 0} 
    if (any(grepl("STAG2", mutation_profile))) {time.series[13,col] <- 0}    
    if (any(grepl("SH2B3", mutation_profile))) {time.series[14,col] <- 0} 
    if (any(grepl("PTPN6", mutation_profile))) {time.series[15,col] <- 0}
    if (any(grepl("AML1_ETO", mutation_profile))) {time.series[16,col] <- 1}	   		  
    if (any(grepl("PML_RARalpha", mutation_profile))) {time.series[17,col] <- 1}	  #check the lettering!
    if (any(grepl("PTPN11", mutation_profile))) {time.series[18,col] <- 1}	  
    if (any(grepl("CUX1", mutation_profile))) {time.series[19,col] <- 1}   
    if (any(grepl("RAD21", mutation_profile))) {time.series[20,col] <- 1} 
    if (any(grepl("BCL2L1", mutation_profile))) {time.series[21,col] <- 1}   
    if (any(grepl("JAK2", mutation_profile))) {time.series[22,col] <- 1}  
    if (any(grepl("CCNA1", mutation_profile))) {time.series[23,col] <- 1}	 
    if (any(grepl("CEBPA", mutation_profile))) {time.series[24,col] <- 0}	  
    if (any(grepl("PIK3IP1", mutation_profile))) {time.series[25,col] <- 0}	  
    if (any(grepl("NRAS", mutation_profile))) {time.series[26,col] <- 1}   
    if (any(grepl("GATA2", mutation_profile))) {time.series[27,col] <- 1}
    if (any(grepl("BAX", mutation_profile))) {time.series[28,col] <- 0}   
    if (any(grepl("STAT5A", mutation_profile))) {time.series[29,col] <- 1}  
    if (any(grepl("SOX4", mutation_profile))) {time.series[30,col] <- 1}		
    if (any(grepl("PI3K", mutation_profile))) {time.series[31,col] <- 1}
    if (any(grepl("SPI1", mutation_profile))) {time.series[32,col] <- 0}	  
    if (any(grepl("PIM", mutation_profile))) {time.series[33,col] <- 1}  
    if (any(grepl("CTNNB1", mutation_profile))) {time.series[34,col] <- 1}    
    if (any(grepl("AKT", mutation_profile))) {time.series[35,col] <- 1}   
    if (any(grepl("GATA1", mutation_profile))) {time.series[36,col] <- 1}  
    if (any(grepl("NUP98_Fusion", mutation_profile))) {time.series[37,col] <- 0}	   	#check the lettering!	  
    if (any(grepl("MLL_Fusion", mutation_profile))) {time.series[38,col] <- 1}	  #check the lettering!
    if (any(grepl("FOXO", mutation_profile))) {time.series[39,col] <- 0}	  
    if (any(grepl("EP300", mutation_profile))) {time.series[40,col] <- 0}  	
    if (any(grepl("BRAF", mutation_profile))) {time.series[41,col] <- 1}    
    if (any(grepl("CBFB", mutation_profile))) {time.series[42,col] <- 0}   
    if (any(grepl("CDK6", mutation_profile))) {time.series[43,col] <- 1}  
    if (any(grepl("MECOM", mutation_profile))) {time.series[44,col] <- 1}	   		  
    if (any(grepl("CBFbeta_MYH11", mutation_profile))) {time.series[45,col] <- 0}	 #check the lettering! 
    if (any(grepl("NPM1", mutation_profile))) {time.series[46,col] <- 1}
    if (any(grepl("IDH1", mutation_profile))) {time.series[47,col] <- 1}   
    if (any(grepl("IDH2", mutation_profile))) {time.series[48,col] <- 1}    
    if (any(grepl("ERK1/2", mutation_profile))) {time.series[49,col] <- 1}   #check the lettering!	
    if (any(grepl("RUNX1", mutation_profile))) {time.series[50,col] <- 0} 	
    if (any(grepl("DOT1L", mutation_profile))) {time.series[51,col] <- 1}	
    if (any(grepl("FBXW7", mutation_profile))) {time.series[52,col] <- 0}	  
    if (any(grepl("AMPK", mutation_profile))) {time.series[53,col] <- 1}	#check the value  
    if (any(grepl("ASXL2", mutation_profile))) {time.series[54,col] <- 0}   
    if (any(grepl("OXO2", mutation_profile))) {time.series[55,col] <- 1}    #check the value
    if (any(grepl("ETV6", mutation_profile))) {time.series[56,col] <- 0}   
    if (any(grepl("AP1", mutation_profile))) {time.series[57,col] <- 1}  
    if (any(grepl("DNMT3A", mutation_profile))) {time.series[58,col] <- 0}	
    if (any(grepl("AEZH2", mutation_profile))) {time.series[59,col] <- 0}
    if (any(grepl("ASXL1", mutation_profile))) {time.series[60,col] <- 0}	
    if (any(grepl("MYC", mutation_profile))) {time.series[61,col] <- 1}   
    if (any(grepl("MTOR", mutation_profile))) {time.series[62,col] <- 1}    #check the value
    if (any(grepl("TET2", mutation_profile))) {time.series[63,col] <- 0}   
    if (any(grepl("PHF6", mutation_profile))) {time.series[64,col] <- 0}  
    if (any(grepl("CCND1", mutation_profile))) {time.series[65,col] <- 1}	
    if (any(grepl("HOXA9", mutation_profile))) {time.series[66,col] <- 1}	  
    if (any(grepl("CDKN2A", mutation_profile))) {time.series[67,col] <- 0}	  
    if (any(grepl("SRSF2", mutation_profile))) {time.series[68,col] <- 1}   #check the value
    if (any(grepl("WT1", mutation_profile))) {time.series[69,col] <- 0}  	
    if (any(grepl("BCOR", mutation_profile))) {time.series[70,col] <- 0}	
    if (any(grepl("UBTF", mutation_profile))) {time.series[71,col] <- 1}  
    if (any(grepl("MEIS1", mutation_profile))) {time.series[72,col] <- 1}	   		  
    if (any(grepl("MDM2", mutation_profile))) {time.series[73,col] <- 1}	 #check the value  
    if (any(grepl("U2AF1", mutation_profile))) {time.series[74,col] <- 1}	  
    if (any(grepl("XPO1", mutation_profile))) {time.series[75,col] <- 1}   
    if (any(grepl("CREBBP", mutation_profile))) {time.series[76,col] <- 0}    
    if (any(grepl("TP53", mutation_profile))) {time.series[77,col] <- 0}   
    if (any(grepl("BCL2", mutation_profile))) {time.series[78,col] <- 1}   
    
  }
  
  print("completed time series")
  
  # adjust the time-series to long-format
  time.series <- t(time.series)
  time.series <- data.frame(time.series)
  names(time.series) <- c("KITLG", "FLT3LG", "GRB2", 
                          "KDM5A", "PTEN", "BCR_ABL", 
                          "CBLB", "CBL", "SOS1", 
                          "KIT", "FLT3", "CDK1", 
                          "STAG2", "SH2B3", "PTPN6", 
                          "AML1_ETO", "PML_RARalpha","PTPN11",						  
                          "CUX1", "RAD21", "BCL2L1",						  
                          "JAK2", "CCNA1", "CEBPA", 
                          "PIK3IP1", "NRAS", "GATA2", 						  
                          "BAX", "STAT5A", "SOX4", 
                          "PI3K", "SPI1", "PIM", 
                          "CTNNB1", "AKT","GATA1", 						  
                          "NUP98_Fusion", "MLL_Fusion", "FOXO", 
                          "EP300", "BRAF", "CBFB", 						  
                          "CDK6", "MECOM", "CBFbeta_MYH11", 						  
                          "NPM1", "IDH1", "IDH2", 						  
                          "ERK1/2", "RUNX1", "DOT1L", 						  
                          "FBXW7", "AMPK","ASXL2",						  
                          "OXO2", "ETV6", "AP1", 
                          "DNMT3A", "AEZH2", "ASXL1", 
                          "MYC", "MTOR", "TET2", 						  
                          "PHF6", "CCND1", "HOXA9", 
                          "CDKN2A", "SRSF2", "WT1", 
                          "BCOR", "UBTF","MEIS1", 						  
                          "MDM2", "U2AF1", "XPO1", 
                          "CREBBP", "TP53", "BCL2")
  
  #getting the network score and saving it as a vector
  time.series$Proliferation <- with(time.series, (-1 *CDKN2A) + (-1 * WT1) + (-1 * BCOR) + STAT5A + CTNNB1 + MYC + AP1 + CCNA1 + MEIS1 + CCND1 + UBTF)
  time.series$Differentiation <- with(time.series, ((-1 * SOX4) + ( -1 * CBFbeta_MYH11) + SPI1 + RUNX1 + CEBPA + ETV6 + EP300))
  time.series$Apoptosis <- with(time.series, ((-1 * BCL2) + (-1 * PIM) + FOXO + TP53 + BAX))
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  print("finished getting the running means")
  
  #final score based on window value of running mean 
  upper_bound <- n.obs  
  lower_bound <- n.obs-1500 
  
  final_score <- mean(time.series$running_mean[upper_bound:lower_bound], na.rm=TRUE)
  print(paste0("the final score is ", final_score))
  
  return(final_score)    
}

############# for finding the convergence #############

#rolling median of the means 
#rolling_median <- rollmedian(time.series$running_mean, k = window)

#diagrams
hist(time.series[window:n.obs, p],
xlim = c(-1,1),
ylim = c(0, n.obs),
main = "Histogram of running mean",
xlab = "running mean"
) # plot the distribution of the running mean 

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
counter <- 0

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
  
  counter <- counter + 1
  
  print(paste0("point ", counter))
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
  labs(title = "Patient PB Blast % vs Score (Fig 3)", x = "% PB Blast", y = "Score") +
  theme_minimal()
PB_scatterplot 

BM_scatterplot <- ggplot(data = BM_patient_df, aes(x = BM_Blast, y = scores)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Patient BM Blast % vs Score (Fig 3)", x = "% BM Blast", y = "Score") +
  theme_minimal()
BM_scatterplot 

#doing a p-test 
PB_pearson <- cor(PB_patient_df$PB_Blast, PB_patient_df$score)
print(PB_pearson)

BM_pearson <- cor(BM_patient_df$BM_Blast, BM_patient_df$score)
print(BM_pearson)