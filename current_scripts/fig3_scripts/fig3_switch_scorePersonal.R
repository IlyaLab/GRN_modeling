### Created by Tazein Shah
### Last modified on 10/4/2023
### This code runs the whole AML network from fig 3 in the Palma paper
### Uses a function called scorePersonal() to knock in/out genes based on patient mutation profiles
### HOWEVER, this uses a switch() versus for loop in original scorePersonal() script (fig3_scorePersonal.R)

#load needed packages
library(BoolNet)
library(reshape2)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggrepel)
library(writexl)
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

n.obs <- 50000 # number of observation aka steps 
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
    
  }
  
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
                          "ERK", "RUNX1", "DOT1L", 						  
                          "FBXW7", "AMPK","ASXL2",						  
                          "OXO2", "ETV6", "AP1", 
                          "DNMT3A", "AEZH2", "ASXL1", 
                          "MYC", "MTOR", "TET2", 						  
                          "PHF6", "CCND1", "HOXA9", 
                          "CDKN2A", "SRSF2", "WT1", 
                          "BCOR", "UBTF","MEIS1", 						  
                          "MDM2", "U2AF1", "XPO1", 
                          "CREBBP", "TP53", "BCL2")
  
  #fixing the genes that are mutated
  #knock out whole columns based on mutation profile
  for (i in mutation_profile){
    switch(i, 
           "KITLG"={time.series$KITLG <- 1},  
           "FLT3LG"={time.series$FLT3LG <- 1},	   		  
           "GRB2"={time.series$GRB2 <- 1},  
           "KDM5A"={time.series$KDM5A <- 1},	  
           "PTEN"={time.series$PTEN <- 0},   
           "BCR_ABL"={time.series$BCR_ABL <- 1},    
           "CBLB"={time.series$CBLB <- 0},   
           "CBL"={time.series$CBL <- 0},  
           "SOS1"={time.series$SOS1 <- 1},	   		  
           "KIT"={time.series$KIT <- 1},		
           "FLT3"={time.series$FLT3 <- 1},	  
           "CDK1"={time.series$CDK1 <- 0}, 
           "STAG2"={time.series$STAG2 <- 0},    
           "SH2B3"={time.series$SH2B3 <- 0}, 
           "PTPN6"={time.series$PTPN6 <- 0},  
           "AML1_ETO"={time.series$AML1_ETO <- 1},	   		  
           "PML_RARalpha"={time.series$PML_RARalpha <- 1},	  #check the lettering!
           "PTPN11"={time.series$PTPN11 <- 1},  
           "CUX1"={time.series$CUX1 <- 1},   
           "RAD21"={time.series$RAD21 <- 1},  	
           "BCL2L1"={time.series$BCL2L1 <- 1},   
           "JAK2"={time.series$JAK2 <- 1},  
           "CCNA1"={time.series$CCNA1 <- 1},	   		  
           "CEBPA"={time.series$CEBPA <- 0},	  
           "PIK3IP1"={time.series$PIK3IP1 <- 0},	  
           "NRAS"={time.series$NRAS <- 1},   
           "GATA2"={time.series$GATA2 <- 1},    
           "BAX"={time.series$BAX <- 0},   
           "STAT5A"={time.series$STAT5A <- 1}, 
           "SOX4"={time.series$SOX4 <- 1},		
           "PI3K"={time.series$PI3K <- 1},	  
           "SPI1"={time.series$SPI1 <- 0},	  
           "PIM"={time.series$PIM <- 1},   
           "CTNNB1"={time.series$CTNNB1 <- 1},    
           "AKT"={time.series$AKT <- 1},   
           "GATA1"={time.series$GATA1 <- 1},  
           "NUP98_Fusion"={time.series$NUP98_Fusion <- 0},	   	#check the lettering!	  
           "MLL_Fusion"={time.series$MLL_Fusion <- 1},	  #check the lettering!
           "FOXO"={time.series$FOXO <- 0},	  
           "EP300"={time.series$EP300 <- 0},  	
           "BRAF"={time.series$BRAF <- 1},    
           "CBFB"={time.series$CBFB <- 0},   
           "CDK6"={time.series$CDK6 <- 1},  
           "MECOM"={time.series$MECOM <- 1},	   		  
           "CBFbeta_MYH11"={time.series$CBFbeta_MYH11 <- 0},	 #check the lettering! 
           "NPM1"={time.series$NPM1 <- 1},	  
           "IDH1"={time.series$IDH1 <- 1},   
           "IDH2"={time.series$IDH2 <- 1},    
           "ERK1/2"={time.series$ERK <- 1},   #check the lettering!	
           "RUNX1"={time.series$RUNX1 <- 0}, 	
           "DOT1L"={time.series$DOT1L <- 1},	   		  
           "FBXW7"={time.series$FBXW7 <- 0},	  
           "AMPK"={time.series$AMPK <- 1},	#check the value  
           "ASXL2"={time.series$ASXL2 <- 0},   
           "OXO2"={time.series$OXO2 <- 1},    #check the value
           "ETV6"={time.series$ETV6 <- 0},   
           "AP1"={time.series$AP1 <- 1},  
           "DNMT3A"={time.series$DNMT3A <- 0},	
           "AEZH2"={time.series$AEZH2 <- 0},	  
           "ASXL1"={time.series$ASXL1 <- 0},	
           "MYC"={time.series$MYC <- 1},   
           "MTOR"={time.series$MTOR <- 1},    #check the value
           "TET2"={time.series$TET2 <- 0},   
           "PHF6"={time.series$PHF6 <- 0},  
           "CCND1"={time.series$CCND1 <- 1},	   		  
           "HOXA9"={time.series$HOXA9 <- 1},	  
           "CDKN2A"={time.series$CDKN2A <- 0},	  
           "SRSF2"={time.series$SRSF2 <- 1},   #check the value
           "WT1"={time.series$WT1 <- 0},  	
           "BCOR"={time.series$BCOR <- 0},	
           "UBTF"={time.series$UBTF <- 1},  
           "MEIS1"={time.series$MEIS1 <- 1},	   		  
           "MDM2"={time.series$MDM2 <- 1},	 #check the value  
           "U2AF1"={time.series$U2AF1 <- 1},	  
           "XPO1"={time.series$XPO1 <- 1},   
           "CREBBP"={time.series$CREBBP <- 0},    
           "TP53"={time.series$TP53 <- 0},   
           "BCL2"={time.series$BCL2 <- 1},   
           {
             
           }
    ) 
  } 
  
  #getting the network score and saving it as a vector
  time.series$Proliferation <- with(time.series, ((-1 * CDKN2A) + (-1 * WT1) + (-1 * BCOR) + STAT5A + CTNNB1 + MYC + AP1 + CCNA1 + MEIS1 + CCND1 + UBTF))
  time.series$Differentiation <- with(time.series, ((-1 * SOX4) + ( -1 * CBFbeta_MYH11) + SPI1 + RUNX1 + CEBPA + ETV6 + EP300))
  time.series$Apoptosis <- with(time.series, ((-1 * BCL2) + (-1 * PIM) + FOXO + TP53 + BAX))
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  print("finished getting the running means")
  
  #final score based on window value of running mean 
  upper_bound <- 50000  
  lower_bound <- upper_bound-1000   
  
  final_score <- mean(time.series$running_mean[lower_bound:upper_bound], na.rm=TRUE)
  print(paste0("the final score is ", final_score))
  
  return(final_score)    
}

############# for finding the convergence #############

#rolling median of the means 
#rolling_median <- rollmedian(time.series$running_mean, k = window)

#diagrams
#hist(time.series[window:n.obs,19],
#xlim = c(0,1),
#ylim = c(0, n.obs),
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
#time <- c(1000, 2000, 4000, 8000, 16000, 32000)

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
write_xlsx(
PB_patient_df,"path/PB_switch_patient.xlsx",
col_names = TRUE,
format_headers = TRUE)

write_xlsx(
  BM_patient_df,"path/BM_switch_patient.xlsx",
  col_names = TRUE,
  format_headers = TRUE)

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