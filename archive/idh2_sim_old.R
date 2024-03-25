### Created by Tazein Shah
### Last modified: 7/19/2023
### This code was used to model the FLT3-NPM1-IDH subgroup from the Palma paper, however, because of a feedback loop, there were cyclic attractors. Decided to not use this method.

#creating the file for the FLT3-NPM1-IDH subgroup
#load BoolNet
library(BoolNet)

#creating the file for the  FLT3-NPM1-IDH
fil <- tempfile(pattern = "testNet")
sink(fil)
cat("targets, factors\n")
cat("FLT3, FLT3\n")
cat("CEBPA, !FLT3\n")
cat("RUNX1, FLT3\n")
cat("RAD21, !RUNX1\n")
cat("PTPN11, FLT3\n")
#look into AKT
cat("AKT, AKT\n")
cat("NRAS, PTPN11\n")
cat("BRAF, NRAS\n")
#look into AMPK
cat("AMPK, AMPK\n") 
cat("NPM1, NPM1\n")
cat("MEK, BRAF\n")
cat("FOXO, !AKT & AMPK)\n")
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

#simulates the flt3_idh_net to get the mutation profiles
#when IDH2 is mutated it turns into 1
#the feedback loop affects NPM1 which means when it's mutated it's cyclic

#WT
figure_a <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(0,1,0))
#NPM1
figure_b <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(0,0,0))
#FLT3
figure_c <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(1,1,0))
#IDH2
figure_d <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(0,1,1))
#FLT3-NPM1
figure_e <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(1,0,0))
#FLT3-IDH2
figure_f <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(1,1,1))
#NPM1-IDH2
figure_g <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(0,0,1))
#FLT3-NPM1-IDH2
figure_h <- fixGenes(flt3_idh_net, fixIndices=c("FLT3","NPM1","IDH2"), values=c(1,0,1))

#getting the attractors of the new mutated profiles 
attr_a <- getAttractors(figure_a)
attr_b <- getAttractors(figure_b)
attr_c <- getAttractors(figure_c)
attr_d <- getAttractors(figure_d)
attr_e <- getAttractors(figure_e)
attr_f <- getAttractors(figure_f)
attr_g <- getAttractors(figure_g)
attr_h <- getAttractors(figure_h)

#load BoolNet
library(BoolNet)

#plot attractors to better see each value 
matrix_a <- as.data.frame(plotAttractors(attr_a)$'1')
matrix_b <- as.data.frame(plotAttractors(attr_b)$'4')
matrix_c <- as.data.frame(plotAttractors(attr_c)$'1')
matrix_d <- as.data.frame(plotAttractors(attr_d)$'1')
matrix_e <- as.data.frame(plotAttractors(attr_e)$'4')
matrix_f <- as.data.frame(plotAttractors(attr_f)$'1')
matrix_g <- as.data.frame(plotAttractors(attr_g)$'4')
matrix_h <- as.data.frame(plotAttractors(attr_h)$'4')

#from the loop, use if you want to save proliferation, differentiation, and apoptosis
# Save the apoptosis as a separate object
#score_name <- paste0("apoptosis_", letters[i])
#score_list[[score_name]] <- apoptosis
#assign(score_name, score_list[[score_name]])

# Save the proliferation as a separate object
#score_name <- paste0("proliferation_", letters[i])
#score_list[[score_name]] <- proliferation
#assign(score_name, score_list[[score_name]])

# Save the differentiation as a separate object
#score_name <- paste0("differentiation_", letters[i])
#score_list[[score_name]] <- differentiation
#assign(score_name, score_list[[score_name]])

score_list <- list()  # List to store the network phenotype scores

#trying to make the loop!
for (i in 1:6) {
  # Select the appropriate color_matrix based on the loop index
  matrix <- switch(i,
                   matrix_b,
                   matrix_c,
                   matrix_e,
                   matrix_f,
                   matrix_g,
                   matrix_h,
  )
  
  #attractor 1
  if (all(c("Attr1.1") %in% colnames(matrix))) {
    # Perform the calculations
    proliferation <- (-1 * (matrix["FOXO", "Attr1.1"])) + (-1 * (matrix["WT1", "Attr1.1"])) +
      (-1 * (matrix["CDKN2A", "Attr1.1"])) + (-1 * (matrix["TP53", "Attr1.1"])) +
      ((matrix["MYC", "Attr1.1"]))
    print(proliferation)
    
    differentiation <- (matrix["CEBPA", "Attr1.1"]) +
      (matrix["RUNX1", "Attr1.1"])
    print(differentiation)
    
    apoptosis <- ((-1 * matrix["BCL2", "Attr1.1"])) +
      ((matrix["TP53", "Attr1.1"]))
    print(apoptosis)
    
    attr_1_network_score <- proliferation - (apoptosis + differentiation)
  }
  
  #attractor 2
  if (all(c("Attr2.1") %in% colnames(matrix))) {
    # Perform the calculations
    proliferation <- (-1 * (matrix["FOXO", "Attr2.1"])) + (-1 * (matrix["WT1", "Attr2.1"])) +
      (-1 * (matrix["CDKN2A", "Attr2.1"])) + (-1 * (matrix["TP53", "Attr2.1"])) +
      ((matrix["MYC", "Attr2.1"]))
    print(proliferation)
    
    differentiation <- (matrix["CEBPA", "Attr2.1"]) +
      (matrix["RUNX1", "Attr2.1"])
    print(differentiation)
    
    apoptosis <- ((-1 * matrix["BCL2", "Attr2.1"])) +
      ((matrix["TP53", "Attr2.1"]))
    print(apoptosis)
    
    attr_2_network_score <- proliferation - (apoptosis + differentiation)
  }
  
  #attractor 3
  if (all(c("Attr3.1") %in% colnames(matrix))) {
    # Perform the calculations
    proliferation <- (-1 * (matrix["FOXO", "Attr3.1"])) + (-1 * (matrix["WT1", "Attr3.1"])) +
      (-1 * (matrix["CDKN2A", "Attr3.1"])) + (-1 * (matrix["TP53", "Attr3.1"])) +
      ((matrix["MYC", "Attr3.1"]))
    print(proliferation)
    
    differentiation <- (matrix["CEBPA", "Attr3.1"]) +
      (matrix["RUNX1", "Attr3.1"])
    print(differentiation)
    
    apoptosis <- ((-1 * matrix["BCL2", "Attr3.1"])) +
      ((matrix["TP53", "Attr3.1"]))
    print(apoptosis)
    
    attr_3_network_score <- proliferation - (apoptosis + differentiation)
  }
  
  #attractor 4
  if (all(c("Attr4.1") %in% colnames(matrix))){
    # Perform the calculations
    proliferation <- (-1 * (matrix["FOXO", "Attr4.1"])) + (-1 * (matrix["WT1", "Attr4.1"])) +
      (-1 * (matrix["CDKN2A", "Attr4.1"])) + (-1 * (matrix["TP53", "Attr4.1"])) +
      ((matrix["MYC", "Attr4.1"]))
    print(proliferation)
    
    differentiation <- (matrix["CEBPA", "Attr4.1"]) +
      (matrix["RUNX1", "Attr4.1"])
    print(differentiation)
    
    apoptosis <- ((-1 * matrix["BCL2", "Attr4.1"])) +
      ((matrix["TP53", "Attr4.1"]))
    print(apoptosis)
    
    attr_4_network_score <- proliferation - (apoptosis + differentiation)
  }
  network_phenotype_score <- mean(attr_1_network_score,attr_2_network_score,attr_3_network_score,attr_4_network_score)
  
  # Save the network phenotype score as a separate object
  score_name <- paste0("network_score_", letters[i+1])
  score_list[[score_name]] <- network_phenotype_score
  assign(score_name, score_list[[score_name]])
}


#for the cyclic attractors
for (i in 1:2) {
  # Select the appropriate color_matrix based on the loop index
  matrix <- switch(i,
                   matrix_a,
                   matrix_d)
  
  #attractor 1 cyclic
  if (all(c("Attr1.2") %in% colnames(matrix))) {
    # Perform the calculations
    #proliferation
    proliferation_1 <- (-1 * (matrix["FOXO", "Attr1.1"])) + (-1 * (matrix["WT1", "Attr1.1"])) +
      (-1 * (matrix["CDKN2A", "Attr1.1"])) + (-1 * (matrix["TP53", "Attr1.1"])) +
      ((matrix["MYC", "Attr1.1"]))
    
    proliferation_2 <- (-1 * (matrix["FOXO", "Attr1.2"])) + (-1 * (matrix["WT1", "Attr1.2"])) +
      (-1 * (matrix["CDKN2A", "Attr1.2"])) + (-1 * (matrix["TP53", "Attr1.2"])) +
      ((matrix["MYC", "Attr1.2"]))
    
    proliferation_3 <- (-1 * (matrix["FOXO", "Attr1.3"])) + (-1 * (matrix["WT1", "Attr1.3"])) +
      (-1 * (matrix["CDKN2A", "Attr1.3"])) + (-1 * (matrix["TP53", "Attr1.3"])) +
      ((matrix["MYC", "Attr1.3"]))
    
    proliferation_4 <- (-1 * (matrix["FOXO", "Attr1.4"])) + (-1 * (matrix["WT1", "Attr1.4"])) +
      (-1 * (matrix["CDKN2A", "Attr1.4"])) + (-1 * (matrix["TP53", "Attr1.4"])) +
      ((matrix["MYC", "Attr1.4"]))
    
    proliferation <- mean(proliferation_1,proliferation_2,proliferation_3,proliferation_4)
    #print(proliferation)
    
    #differentiation
    differentiation_1 <- (matrix["CEBPA", "Attr1.1"]) + (matrix["RUNX1", "Attr1.1"])
    differentiation_2 <- (matrix["CEBPA", "Attr1.2"]) + (matrix["RUNX1", "Attr1.2"])
    differentiation_3 <- (matrix["CEBPA", "Attr1.3"]) + (matrix["RUNX1", "Attr1.3"])
    differentiation_4 <- (matrix["CEBPA", "Attr1.4"]) + (matrix["RUNX1", "Attr1.4"])
    
    differentiation <- mean(differentiation_1,differentiation_2,differentiation_3,differentiation_4)
    #print(differentiation)
    
    #apoptosis
    apoptosis_1 <- ((-1 * matrix["BCL2", "Attr1.1"])) + ((matrix["TP53", "Attr1.1"]))
    apoptosis_2 <- ((-1 * matrix["BCL2", "Attr1.2"])) + ((matrix["TP53", "Attr1.2"]))
    apoptosis_3 <- ((-1 * matrix["BCL2", "Attr1.3"])) + ((matrix["TP53", "Attr1.3"]))
    apoptosis_4 <- ((-1 * matrix["BCL2", "Attr1.4"])) + ((matrix["TP53", "Attr1.4"]))
    
    apoptosis <- mean(apoptosis_1,apoptosis_2,apoptosis_3,apoptosis_4)
    #print(apoptosis)
    
    #network phenotype score
    attr_1_cyclic_network_score <- proliferation - (apoptosis + differentiation)
  }
  
  #attractor 2 cyclic
  if (all(c("Attr2.2") %in% colnames(matrix))) {
    # Perform the calculations
    #proliferation
    proliferation_1 <- (-1 * (matrix["FOXO", "Attr2.1"])) + (-1 * (matrix["WT1", "Attr2.1"])) +
      (-1 * (matrix["CDKN2A", "Attr2.1"])) + (-1 * (matrix["TP53", "Attr2.1"])) +
      ((matrix["MYC", "Attr2.1"]))
    
    proliferation_2 <- (-1 * (matrix["FOXO", "Attr2.2"])) + (-1 * (matrix["WT1", "Attr2.2"])) +
      (-1 * (matrix["CDKN2A", "Attr2.2"])) + (-1 * (matrix["TP53", "Attr2.2"])) +
      ((matrix["MYC", "Attr2.2"]))
    
    proliferation_3 <- (-1 * (matrix["FOXO", "Attr2.3"])) + (-1 * (matrix["WT1", "Attr2.3"])) +
      (-1 * (matrix["CDKN2A", "Attr2.3"])) + (-1 * (matrix["TP53", "Attr2.3"])) +
      ((matrix["MYC", "Attr2.3"]))
    
    proliferation_4 <- (-1 * (matrix["FOXO", "Attr2.4"])) + (-1 * (matrix["WT1", "Attr2.4"])) +
      (-1 * (matrix["CDKN2A", "Attr2.4"])) + (-1 * (matrix["TP53", "Attr2.4"])) +
      ((matrix["MYC", "Attr2.4"]))
    
    proliferation <- mean(proliferation_1,proliferation_2,proliferation_3,proliferation_4)
    #print(proliferation)
    
    #differentiation
    differentiation_1 <- (matrix["CEBPA", "Attr2.1"]) + (matrix["RUNX1", "Attr2.1"])
    differentiation_2 <- (matrix["CEBPA", "Attr2.2"]) + (matrix["RUNX1", "Attr2.2"])
    differentiation_3 <- (matrix["CEBPA", "Attr2.3"]) + (matrix["RUNX1", "Attr2.3"])
    differentiation_4 <- (matrix["CEBPA", "Attr2.4"]) + (matrix["RUNX1", "Attr2.4"])
    
    differentiation <- mean(differentiation_1,differentiation_2,differentiation_3,differentiation_4)
    #print(differentiation)
    
    #apoptosis
    apoptosis_1 <- ((-1 * matrix["BCL2", "Attr2.1"])) + ((matrix["TP53", "Attr2.1"]))
    apoptosis_2 <- ((-1 * matrix["BCL2", "Attr2.2"])) + ((matrix["TP53", "Attr2.2"]))
    apoptosis_3 <- ((-1 * matrix["BCL2", "Attr2.3"])) + ((matrix["TP53", "Attr2.3"]))
    apoptosis_4 <- ((-1 * matrix["BCL2", "Attr2.4"])) + ((matrix["TP53", "Attr2.4"]))
    
    apoptosis <- mean(apoptosis_1,apoptosis_2,apoptosis_3,apoptosis_4)
    #print(apoptosis)
    
    #network phenotype score
    attr_2_cyclic_network_score <- proliferation - (apoptosis + differentiation)
  }
  
  #attractor 3 cyclic
  if (all(c("Attr3.2") %in% colnames(matrix))) {
    # Perform the calculations
    #proliferation
    proliferation_1 <- (-1 * (matrix["FOXO", "Attr3.1"])) + (-1 * (matrix["WT1", "Attr3.1"])) +
      (-1 * (matrix["CDKN2A", "Attr3.1"])) + (-1 * (matrix["TP53", "Attr3.1"])) +
      ((matrix["MYC", "Attr3.1"]))
    
    proliferation_2 <- (-1 * (matrix["FOXO", "Attr3.2"])) + (-1 * (matrix["WT1", "Attr3.2"])) +
      (-1 * (matrix["CDKN2A", "Attr3.2"])) + (-1 * (matrix["TP53", "Attr3.2"])) +
      ((matrix["MYC", "Attr3.2"]))
    
    proliferation_3 <- (-1 * (matrix["FOXO", "Attr3.3"])) + (-1 * (matrix["WT1", "Attr3.3"])) +
      (-1 * (matrix["CDKN2A", "Attr3.3"])) + (-1 * (matrix["TP53", "Attr3.3"])) +
      ((matrix["MYC", "Attr3.3"]))
    
    proliferation_4 <- (-1 * (matrix["FOXO", "Attr3.4"])) + (-1 * (matrix["WT1", "Attr3.4"])) +
      (-1 * (matrix["CDKN2A", "Attr3.4"])) + (-1 * (matrix["TP53", "Attr3.4"])) +
      ((matrix["MYC", "Attr3.4"]))
    
    proliferation <- mean(proliferation_1,proliferation_2,proliferation_3,proliferation_4)
    #print(proliferation)
    
    #differentiation
    differentiation_1 <- (matrix["CEBPA", "Attr3.1"]) + (matrix["RUNX1", "Attr3.1"])
    differentiation_2 <- (matrix["CEBPA", "Attr3.2"]) + (matrix["RUNX1", "Attr3.2"])
    differentiation_3 <- (matrix["CEBPA", "Attr3.3"]) + (matrix["RUNX1", "Attr3.3"])
    differentiation_4 <- (matrix["CEBPA", "Attr3.4"]) + (matrix["RUNX1", "Attr3.4"])
    
    differentiation <- mean(differentiation_1,differentiation_2,differentiation_3,differentiation_4)
    #print(differentiation)
    
    #apoptosis
    apoptosis_1 <- ((-1 * matrix["BCL2", "Attr3.1"])) + ((matrix["TP53", "Attr3.1"]))
    apoptosis_2 <- ((-1 * matrix["BCL2", "Attr3.2"])) + ((matrix["TP53", "Attr3.2"]))
    apoptosis_3 <- ((-1 * matrix["BCL2", "Attr3.3"])) + ((matrix["TP53", "Attr3.3"]))
    apoptosis_4 <- ((-1 * matrix["BCL2", "Attr3.4"])) + ((matrix["TP53", "Attr3.4"]))
    
    apoptosis <- mean(apoptosis_1,apoptosis_2,apoptosis_3,apoptosis_4)
    #print(apoptosis)
    
    #network phenotype score
    attr_3_cyclic_network_score <- proliferation - (apoptosis + differentiation)
  }  
  
  # attractor 4 cyclic
  if (all(c("Attr4.2") %in% colnames(matrix))) {
    #perform the calculations
    #proliferation
    proliferation_1 <- (-1 * (matrix["FOXO", "Attr4.1"])) + (-1 * (matrix["WT1", "Attr4.1"])) +
      (-1 * (matrix["CDKN2A", "Attr4.1"])) + (-1 * (matrix["TP53", "Attr4.1"])) +
      ((matrix["MYC", "Attr4.1"]))
    
    proliferation_2 <- (-1 * (matrix["FOXO", "Attr4.2"])) + (-1 * (matrix["WT1", "Attr4.2"])) +
      (-1 * (matrix["CDKN2A", "Attr4.2"])) + (-1 * (matrix["TP53", "Attr4.2"])) +
      ((matrix["MYC", "Attr4.2"]))
    
    proliferation_3 <- (-1 * (matrix["FOXO", "Attr4.3"])) + (-1 * (matrix["WT1", "Attr4.3"])) +
      (-1 * (matrix["CDKN2A", "Attr4.3"])) + (-1 * (matrix["TP53", "Attr4.3"])) +
      ((matrix["MYC", "Attr4.3"]))
    
    proliferation_4 <- (-1 * (matrix["FOXO", "Attr1.4"])) + (-1 * (matrix["WT1", "Attr1.4"])) +
      (-1 * (matrix["CDKN2A", "Attr1.4"])) + (-1 * (matrix["TP53", "Attr1.4"])) +
      ((matrix["MYC", "Attr1.4"]))
    
    proliferation <- mean(proliferation_1,proliferation_2,proliferation_3,proliferation_4)
    #print(proliferation)
    
    #differentiation
    differentiation_1 <- (matrix["CEBPA", "Attr4.1"]) + (matrix["RUNX1", "Attr4.1"])
    differentiation_2 <- (matrix["CEBPA", "Attr4.2"]) + (matrix["RUNX1", "Attr4.2"])
    differentiation_3 <- (matrix["CEBPA", "Attr4.3"]) + (matrix["RUNX1", "Attr4.3"])
    differentiation_4 <- (matrix["CEBPA", "Attr1.4"]) + (matrix["RUNX1", "Attr1.4"])
    
    differentiation <- mean(differentiation_1,differentiation_2,differentiation_3,differentiation_4)
    #print(differentiation)
    
    #apoptosis
    apoptosis_1 <- ((-1 * matrix["BCL2", "Attr4.1"])) + ((matrix["TP53", "Attr4.1"]))
    apoptosis_2 <- ((-1 * matrix["BCL2", "Attr4.2"])) + ((matrix["TP53", "Attr4.2"]))
    apoptosis_3 <- ((-1 * matrix["BCL2", "Attr4.3"])) + ((matrix["TP53", "Attr4.3"]))
    apoptosis_4 <- ((-1 * matrix["BCL2", "Attr1.4"])) + ((matrix["TP53", "Attr1.4"]))
    
    apoptosis <- mean(apoptosis_1,apoptosis_2,apoptosis_3,apoptosis_4)
    #print(apoptosis)
    
    #network phenotype score
    attr_4_cyclic_network_score <- proliferation - (apoptosis + differentiation)
  }
  network_phenotype_score <- mean(attr_1_cyclic_network_score,attr_2_cyclic_network_score,
                                  attr_3_cyclic_network_score,attr_4_cyclic_network_score)
  
  #save the network phenotype score as a separate object
  score_name <- paste0("network_score_cyclic_", letters[i])
  score_list[[score_name]] <- network_phenotype_score
  assign(score_name, score_list[[score_name]])}


#now we are making the scatterplots comparing the newtwork scores to clinical data
#load necessary packages for cleaning up the tables, making scatterplots
library(dplyr)
library(stringr)
require(ggplot2)
library(tidyr)
library(ggsci)
library(ggrepel) 

#had to change the column name of FLT3-ITD
#colnames(s5_df)[4] ="FLT3_ITD"

#renaming s5_df column LabId to labId
colnames(s5_df)[1] ="labId"

#getting what's in common between both s5 and s7 labId's
lab_ID_common <- intersect(s5_df$labId,s7_df$labId)

#making new df with only the labId that overlap
s5_common_df <- subset(s5_df, labId %in% c(lab_ID_common))
s7_common_df <- subset(s7_df, labId %in% c(lab_ID_common))

#checking to see if they do share the same values
#sum(str_detect(s7_common_df$labId, '^09-00705$')) > 0
#sum(str_detect(s5_common_df$labId, '^16-00702$')) > 0

#make a merged df of the s5_common_df and s5_common_df
common_merge_df <- s7_common_df %>% full_join(s5_common_df)

#remove duplicate rows
common_df <- distinct(common_merge_df)

#replacing the NA values with the mean
#getting the values into their own vector 
common_PB_Blast <- common_df$`%.Blasts.in.PB`
common_BM_Blast <- common_df$`%.Blasts.in.BM`
#Replace NA with mean 
PB_Blast_replacement <- ifelse(is.na(common_PB_Blast),(mean(common_df$`%.Blasts.in.PB`, na.rm = TRUE)),common_PB_Blast)
BM_Blast_replacement <- ifelse(is.na(common_BM_Blast),(mean(common_df$`%.Blasts.in.BM`, na.rm = TRUE)),common_BM_Blast)

#lets put back the new variable into the dataframe
common_wo_NA_df <- common_df
common_wo_NA_df$`%.Blasts.in.PB` <- PB_Blast_replacement
common_wo_NA_df$`%.Blasts.in.BM` <- BM_Blast_replacement
common_wo_NA_df

#sort by labID name
common_wo_NA_df <- common_wo_NA_df[order(common_wo_NA_df[,'labId']), ]

#combines the symbols rows into symbol column
symbol_common <- common_wo_NA_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#with the combined symbol column, separate the symbols into separate columns
symbol_common <- symbol_common %>%
  separate_rows(symbol, convert = TRUE) %>%
  group_by(labId) %>%
  mutate(col_num = row_number()) %>%
  spread(col_num, symbol, sep = "_") %>%
  ungroup()

#Making all values that are NOT FLT3, NPM1, or IDH2 NA
symbol_common <- symbol_common %>%
  mutate_at(vars(col_num_1:col_num_70), ~ifelse(. %in% c("NPM1", "IDH2","FLT3"), ., NA))

#Recombining all the columns, and the NA values "fall out" 
symbol_common <- symbol_common %>%
  unite(combined_column, col_num_1:col_num_70, sep = ",", na.rm = TRUE)

#grouping the new symbols based on diagram in Fig 5
symbol_common <- symbol_common %>%
  mutate(group_value = 
           ifelse(combined_column == "NPM1", 2, 
                  ifelse(combined_column == "FLT3", 3,
                         ifelse(combined_column == "IDH2", 4, 
                                ifelse(combined_column == "FLT3,NPM1", 5,
                                       ifelse(combined_column == "FLT3,IDH2", 6,
                                              ifelse(combined_column == "IDH2,NPM1", 7,
                                                     ifelse(combined_column == "FLT3,IDH2,NPM1", 8,
                                                            ifelse(combined_column == "", 1,NA)))))))))

#printing the IDs of every group
group_1_id <- symbol_common %>%
  filter(group_value == 1) %>%
  pull(labId)
print(group_1_id)

group_2_id <- symbol_common %>%
  filter(group_value == 2) %>%
  pull(labId)
print(group_2_id)

group_3_id <- symbol_common %>%
  filter(group_value == 3) %>%
  pull(labId)
print(group_3_id)

group_4_id <- symbol_common %>%
  filter(group_value == 4) %>%
  pull(labId)
print(group_4_id)

group_5_id <- symbol_common %>%
  filter(group_value == 5) %>%
  pull(labId)
print(group_5_id)

group_6_id <- symbol_common %>%
  filter(group_value == 6) %>%
  pull(labId)
print(group_6_id)

group_7_id <- symbol_common %>%
  filter(group_value == 7) %>%
  pull(labId)
print(group_7_id)

group_8_id <- symbol_common %>%
  filter(group_value == 8) %>%
  pull(labId)
print(group_8_id)

#getting the averages for the group IDs
WT_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_1_id)$`%.Blasts.in.PB`, na.rm = TRUE)
WT_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_1_id)$`%.Blasts.in.BM`, na.rm = TRUE)

NPM1_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_2_id)$`%.Blasts.in.PB`, na.rm = TRUE)
NPM1_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_2_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_3_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_3_id)$`%.Blasts.in.BM`, na.rm = TRUE)

IDH2_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_4_id)$`%.Blasts.in.PB`, na.rm = TRUE)
IDH2_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_4_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_NPM1_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_5_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_NPM1_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_5_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_IDH2_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_6_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_IDH2_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_6_id)$`%.Blasts.in.BM`, na.rm = TRUE)

NPM1_IDH2_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_7_id)$`%.Blasts.in.PB`, na.rm = TRUE)
NPM1_IDH2_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_7_id)$`%.Blasts.in.BM`, na.rm = TRUE)

all_mutation_PB <- mean(subset(common_wo_NA_df, labId %in% group_8_id)$`%.Blasts.in.PB`, na.rm = TRUE)
all_mutation_BM <- mean(subset(common_wo_NA_df, labId %in% group_8_id)$`%.Blasts.in.BM`, na.rm = TRUE)

#the network_score_ comes from the other R Script
#making dataframes for graphs
#PB Blast vs Network Score
PB_blast_vs_network_score_IDH2 <- data.frame( 
  PB_BLAST = c(get("WT_mean_PB"), get("NPM1_mean_PB"), get("FLT3_mean_PB"), get("IDH2_mean_PB"), get("FLT3_NPM1_mean_PB"), get("FLT3_IDH2_mean_PB"), get("NPM1_IDH2_mean_PB"),get("all_mutation_PB")),
  Network_Score = c(get("network_score_cyclic_a"), get("network_score_b"), get("network_score_c"), get("network_score_cyclic_b"), get("network_score_d"), get("network_score_e"), get("network_score_f"),get("network_score_g")),
  row.names = c('WT', 'NPM1', 'FLT3', 'IDH2', 'FLT3_NPM1', 'FLT3_IDH2', 'NPM1_IDH2', 'FLT3_NPM1_IDH2'))

print(PB_blast_vs_network_score_IDH2)

#the network_score_ comes from the other R Script
#making dataframes for graphs
#BM Blast vs Network Score
BM_blast_vs_network_score_IDH2 <- data.frame( 
  BM_BLAST = c(get("WT_mean_BM"), get("NPM1_mean_BM"), get("FLT3_mean_BM"), get("IDH2_mean_BM"), get("FLT3_NPM1_mean_BM"), get("FLT3_IDH2_mean_BM"), get("NPM1_IDH2_mean_BM"),get("all_mutation_PB")),
  Network_Score = c(get("network_score_cyclic_a"), get("network_score_b"), get("network_score_c"), get("network_score_cyclic_b"), get("network_score_d"), get("network_score_e"), get("network_score_f"),get("network_score_g")),
  row.names = c('WT', 'NPM1', 'FLT3', 'IDH2', 'FLT3_NPM1', 'FLT3_IDH2', 'NPM1_IDH2', 'FLT3_NPM1_IDH2'))

print(BM_blast_vs_network_score_IDH2)

#plotting the PB Blast vs Network Score scatter plot
PB_network_scatterplot_IDH2 <- ggplot(data = PB_blast_vs_network_score_IDH2, aes(x = PB_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(PB_blast_vs_network_score_IDH2)), vjust = -1) +
  labs(title = "% PB Blast vs Network Phenotypes", x = "% PB Blast", y = "Network Phenotype Score") +
  theme_minimal()
PB_network_scatterplot_IDH2 + coord_cartesian(xlim = c(20,85), ylim = c(-6, 0))

BM_network_scatterplot_IDH2 <- ggplot(data = BM_blast_vs_network_score_IDH2, aes(x = BM_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(BM_blast_vs_network_score_IDH2)), vjust = -1) +
  labs(title = "% BM Blast vs Network Phenotypes", x = "% BM Blast", y = "Network Phenotype Score") +
  theme_minimal()
BM_network_scatterplot_IDH2 + coord_cartesian(xlim = c(20,85), ylim = c(-6, 0))

#doing a p test 
PB_diagram <- PB_blast_vs_network_score_IDH2$PB_BLAST
Network_diagram <- PB_blast_vs_network_score_IDH2$Network_Score
BM_diagram <- BM_blast_vs_network_score_IDH2$BM_BLAST

t.test(PB_diagram, Network_diagram, paired = FALSE)
t.test(BM_diagram, Network_diagram, paired = FALSE)

#doing a pearson correlation
PB_network_pearson_IDH2 <- cor(PB_blast_vs_network_score_IDH2$PB_BLAST, PB_blast_vs_network_score_IDH2$Network_Score)
print(PB_network_pearson_IDH2)

BM_network_pearson_IDH2 <- cor(BM_blast_vs_network_score_IDH2$BM_BLAST, BM_blast_vs_network_score_IDH2$Network_Score)
print(BM_network_pearson_IDH2)

