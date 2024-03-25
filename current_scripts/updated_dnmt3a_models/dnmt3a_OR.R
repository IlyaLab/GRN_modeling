### Created by Tazein Shah
### Last modified on 12/12/2023
### In the paper the equation for MYC is MYC = !(GSK3B & FBXW7) & ERK, which is an exception to the rule of "or" seperating inhibitors from inhibitors/ activators from activators. 
### This code tests to see if there are any major differences in the network when this MYC rule is changed to match the standard rules

#load BoolNet
library(BoolNet)

#creating the file for the FLT3-NPM1-DNMT3A subgroup
fil <- tempfile(pattern = "testNet")
sink(fil)
cat("targets, factors\n")
cat("FLT3, FLT3\n")
cat("AKT, FLT3\n")
cat("CEBPA,!FLT3\n")
cat("DNMT3A, DNMT3A\n")
cat("GSK3B, !AKT\n")
cat("NPM1, NPM1\n")
cat("ARF, NPM1\n")
cat("HOXA9, !NPM1\n")
cat("FBXW7, NPM1\n")
cat("ERK, FLT3\n")
cat("CDKN2A, NPM1\n")
cat("STAT5A, FLT3\n")
cat("SOX4, !CEBPA\n")
cat("CCND1, !(DNMT3A | GSK3B)\n")
cat("MEIS1,!(DNMT3A | !HOXA9)\n")
cat("MYC, !(GSK3B | FBXW7) & ERK\n")
cat("ETV6,!ERK\n")
cat("TP53, ARF\n")
cat("BCL2, ERK & !TP53\n")
cat("Apoptosis, TP53 & !BCL2\n")
cat("Differentiation, (CEBPA | ETV6) & !MEIS1 \n")
cat("Proliferation, (MYC| CCND1| SOX4 | MEIS1| STAT5A)\n")

sink()

#this loads the network and saves it as aml_test_net
aml_test_net <- loadNetwork(fil)
print(aml_test_net)

#simulates the aml_test_net to get figure 5 examples
figure_5a <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(0,1,1))
figure_5b <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(0,0,1))
figure_5c <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(1,1,1))
figure_5d <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(0,1,0))
figure_5e <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(1,0,1))
figure_5f <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(1,1,0))
figure_5g <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(0,0,0))
figure_5h <- fixGenes(aml_test_net, fixIndices=c("FLT3","NPM1","DNMT3A"), values=c(1,0,0))

#getting the attractors of the new mutated profiles 
attr_a <- getAttractors(figure_5a)
attr_b <- getAttractors(figure_5b)
attr_c <- getAttractors(figure_5c)
attr_d <- getAttractors(figure_5d)
attr_e <- getAttractors(figure_5e)
attr_f <- getAttractors(figure_5f)
attr_g <- getAttractors(figure_5g)
attr_h <- getAttractors(figure_5h)

#plot attractors to better see each value 
matrix_a <- print(plotAttractors(attr_a)$'1')
matrix_b <- print(plotAttractors(attr_b)$'1')
matrix_c <- print(plotAttractors(attr_c)$'1')
matrix_d <- print(plotAttractors(attr_d)$'1')
matrix_e <- print(plotAttractors(attr_e)$'1')
matrix_f <- print(plotAttractors(attr_f)$'1')
matrix_g <- print(plotAttractors(attr_g)$'1')
matrix_h <- print(plotAttractors(attr_h)$'1')

#making the 0 into -1 for the color score
#color_matrix_a[color_matrix_a==0] <- -1
#color_matrix_b[color_matrix_b==0] <- -1
#color_matrix_c[color_matrix_c==0] <- -1
#color_matrix_d[color_matrix_d==0] <- -1
#color_matrix_e[color_matrix_e==0] <- -1
#color_matrix_f[color_matrix_f==0] <- -1
#color_matrix_g[color_matrix_g==0] <- -1
#color_matrix_h[color_matrix_h==0] <- -1

# Define the loop to create tables
table_list <- list()  # List to store the tables
score_list <- list()  # List to store the network phenotype scores

for (i in 1:8) {
  # Select the appropriate matrix based on the loop index
  matrix <- switch(i,
                   matrix_a,
                   matrix_b,
                   matrix_c,
                   matrix_d,
                   matrix_e,
                   matrix_f,
                   matrix_g,
                   matrix_h)
  
  # Perform the calculations
  proliferation <- ((matrix["STAT5A", "Attr1.1"])) +
    (matrix["SOX4", "Attr1.1"]) +
    (matrix["CCND1", "Attr1.1"]) +
    (matrix["MEIS1", "Attr1.1"]) +
    (matrix["MYC", "Attr1.1"])
  print(proliferation)
  
  differentiation <- (matrix["CEBPA", "Attr1.1"]) +
    (( -1 * matrix["MEIS1", "Attr1.1"])) +
    (matrix["ETV6", "Attr1.1"])
  print(differentiation)
  
  apoptosis <- ((-1 * matrix["BCL2", "Attr1.1"])) +
    ((matrix["TP53", "Attr1.1"]))
  print(apoptosis)
  
  network_phenotype_score <- proliferation - (apoptosis + differentiation)
  
  # Save the apoptosis as a separate object
  score_name <- paste0("apoptosis_", letters[i])
  score_list[[score_name]] <- apoptosis
  assign(score_name, score_list[[score_name]])
  
  # Save the proliferation as a separate object
  score_name <- paste0("proliferation_", letters[i])
  score_list[[score_name]] <- proliferation
  assign(score_name, score_list[[score_name]])
  
  # Save the differentiation as a separate object
  score_name <- paste0("differentiation_", letters[i])
  score_list[[score_name]] <- differentiation
  assign(score_name, score_list[[score_name]])
  
  # Save the network phenotype score as a separate object
  score_name <- paste0("network_score_", letters[i])
  score_list[[score_name]] <- network_phenotype_score
  assign(score_name, score_list[[score_name]])
}

#Making the mega table
Figure = c("figure_a","figure_b","figure_c",
           "figure_d","figure_e","figure_f",
           "figure_g","figure_h")

FLT3 = c(matrix_a["FLT3", "Attr1.1"],
         matrix_b["FLT3", "Attr1.1"],
         matrix_c["FLT3", "Attr1.1"],
         matrix_d["FLT3", "Attr1.1"],
         matrix_e["FLT3", "Attr1.1"],
         matrix_f["FLT3", "Attr1.1"],
         matrix_g["FLT3", "Attr1.1"],
         matrix_h["FLT3", "Attr1.1"])

DNMT3A = c(matrix_a["DNMT3A", "Attr1.1"],
           matrix_b["DNMT3A", "Attr1.1"],
           matrix_c["DNMT3A", "Attr1.1"],
           matrix_d["DNMT3A", "Attr1.1"],
           matrix_e["DNMT3A", "Attr1.1"],
           matrix_f["DNMT3A", "Attr1.1"],
           matrix_g["DNMT3A", "Attr1.1"],
           matrix_h["DNMT3A", "Attr1.1"])

NPM1 = c(matrix_a["NPM1", "Attr1.1"],
         matrix_b["NPM1", "Attr1.1"],
         matrix_c["NPM1", "Attr1.1"],
         matrix_d["NPM1", "Attr1.1"],
         matrix_e["NPM1", "Attr1.1"],
         matrix_f["NPM1", "Attr1.1"],
         matrix_g["NPM1", "Attr1.1"],
         matrix_h["NPM1", "Attr1.1"])

AKT = c(matrix_a["AKT", "Attr1.1"],
        matrix_b["AKT", "Attr1.1"],
        matrix_c["AKT", "Attr1.1"],
        matrix_d["AKT", "Attr1.1"],
        matrix_e["AKT", "Attr1.1"],
        matrix_f["AKT", "Attr1.1"],
        matrix_g["AKT", "Attr1.1"],
        matrix_h["AKT", "Attr1.1"])

CEBPA = c(matrix_a["CEBPA", "Attr1.1"],
          matrix_b["CEBPA", "Attr1.1"],
          matrix_c["CEBPA", "Attr1.1"],
          matrix_d["CEBPA", "Attr1.1"],
          matrix_e["CEBPA", "Attr1.1"],
          matrix_f["CEBPA", "Attr1.1"],
          matrix_g["CEBPA", "Attr1.1"],
          matrix_h["CEBPA", "Attr1.1"])

GSK3B = c(matrix_a["GSK3B", "Attr1.1"],
          matrix_b["GSK3B", "Attr1.1"],
          matrix_c["GSK3B", "Attr1.1"],
          matrix_d["GSK3B", "Attr1.1"],
          matrix_e["GSK3B", "Attr1.1"],
          matrix_f["GSK3B", "Attr1.1"],
          matrix_g["GSK3B", "Attr1.1"],
          matrix_h["GSK3B", "Attr1.1"])

ARF = c(matrix_a["ARF", "Attr1.1"],
        matrix_b["ARF", "Attr1.1"],
        matrix_c["ARF", "Attr1.1"],
        matrix_d["ARF", "Attr1.1"],
        matrix_e["ARF", "Attr1.1"],
        matrix_f["ARF", "Attr1.1"],
        matrix_g["ARF", "Attr1.1"],
        matrix_h["ARF", "Attr1.1"])

HOXA9 = c(matrix_a["HOXA9", "Attr1.1"],
          matrix_b["HOXA9", "Attr1.1"],
          matrix_c["HOXA9", "Attr1.1"],
          matrix_d["HOXA9", "Attr1.1"],
          matrix_e["HOXA9", "Attr1.1"],
          matrix_f["HOXA9", "Attr1.1"],
          matrix_g["HOXA9", "Attr1.1"],
          matrix_h["HOXA9", "Attr1.1"])

FBXW7 = c(matrix_a["FBXW7", "Attr1.1"],
          matrix_b["FBXW7", "Attr1.1"],
          matrix_c["FBXW7", "Attr1.1"],
          matrix_d["FBXW7", "Attr1.1"],
          matrix_e["FBXW7", "Attr1.1"],
          matrix_f["FBXW7", "Attr1.1"],
          matrix_g["FBXW7", "Attr1.1"],
          matrix_h["FBXW7", "Attr1.1"])

ERK = c(matrix_a["ERK", "Attr1.1"],
        matrix_b["ERK", "Attr1.1"],
        matrix_c["ERK", "Attr1.1"],
        matrix_d["ERK", "Attr1.1"],
        matrix_e["ERK", "Attr1.1"],
        matrix_f["ERK", "Attr1.1"],
        matrix_g["ERK", "Attr1.1"],
        matrix_h["ERK", "Attr1.1"])

CDKN2A = c(matrix_a["CDKN2A", "Attr1.1"],
           matrix_b["CDKN2A", "Attr1.1"],
           matrix_c["CDKN2A", "Attr1.1"],
           matrix_d["CDKN2A", "Attr1.1"],
           matrix_e["CDKN2A", "Attr1.1"],
           matrix_f["CDKN2A", "Attr1.1"],
           matrix_g["CDKN2A", "Attr1.1"],
           matrix_h["CDKN2A", "Attr1.1"])

STAT5A = c(matrix_a["STAT5A", "Attr1.1"],
           matrix_b["STAT5A", "Attr1.1"],
           matrix_c["STAT5A", "Attr1.1"],
           matrix_d["STAT5A", "Attr1.1"],
           matrix_e["STAT5A", "Attr1.1"],
           matrix_f["STAT5A", "Attr1.1"],
           matrix_g["STAT5A", "Attr1.1"],
           matrix_h["STAT5A", "Attr1.1"])

SOX4 = c(matrix_a["SOX4", "Attr1.1"],
         matrix_b["SOX4", "Attr1.1"],
         matrix_c["SOX4", "Attr1.1"],
         matrix_d["SOX4", "Attr1.1"],
         matrix_e["SOX4", "Attr1.1"],
         matrix_f["SOX4", "Attr1.1"],
         matrix_g["SOX4", "Attr1.1"],
         matrix_h["SOX4", "Attr1.1"])

CCND1 = c(matrix_a["CCND1", "Attr1.1"],
          matrix_b["CCND1", "Attr1.1"],
          matrix_c["CCND1", "Attr1.1"],
          matrix_d["CCND1", "Attr1.1"],
          matrix_e["CCND1", "Attr1.1"],
          matrix_f["CCND1", "Attr1.1"],
          matrix_g["CCND1", "Attr1.1"],
          matrix_h["CCND1", "Attr1.1"])

MEIS1 = c(matrix_a["MEIS1", "Attr1.1"],
          matrix_b["MEIS1", "Attr1.1"],
          matrix_c["MEIS1", "Attr1.1"],
          matrix_d["MEIS1", "Attr1.1"],
          matrix_e["MEIS1", "Attr1.1"],
          matrix_f["MEIS1", "Attr1.1"],
          matrix_g["MEIS1", "Attr1.1"],
          matrix_h["MEIS1", "Attr1.1"])

MYC = c(matrix_a["MYC", "Attr1.1"],
        matrix_b["MYC", "Attr1.1"],
        matrix_c["MYC", "Attr1.1"],
        matrix_d["MYC", "Attr1.1"],
        matrix_e["MYC", "Attr1.1"],
        matrix_f["MYC", "Attr1.1"],
        matrix_g["MYC", "Attr1.1"],
        matrix_h["MYC", "Attr1.1"])

ETV6 = c(matrix_a["ETV6", "Attr1.1"],
         matrix_b["ETV6", "Attr1.1"],
         matrix_c["ETV6", "Attr1.1"],
         matrix_d["ETV6", "Attr1.1"],
         matrix_e["ETV6", "Attr1.1"],
         matrix_f["ETV6", "Attr1.1"],
         matrix_g["ETV6", "Attr1.1"],
         matrix_h["ETV6", "Attr1.1"])

TP53 = c(matrix_a["TP53", "Attr1.1"],
         matrix_b["TP53", "Attr1.1"],
         matrix_c["TP53", "Attr1.1"],
         matrix_d["TP53", "Attr1.1"],
         matrix_e["TP53", "Attr1.1"],
         matrix_f["TP53", "Attr1.1"],
         matrix_g["TP53", "Attr1.1"],
         matrix_h["TP53", "Attr1.1"])

BCL2 = c(matrix_a["BCL2", "Attr1.1"],
         matrix_b["BCL2", "Attr1.1"],
         matrix_c["BCL2", "Attr1.1"],
         matrix_d["BCL2", "Attr1.1"],
         matrix_e["BCL2", "Attr1.1"],
         matrix_f["BCL2", "Attr1.1"],
         matrix_g["BCL2", "Attr1.1"],
         matrix_h["BCL2", "Attr1.1"])


Apoptosis = c((apoptosis_a), (apoptosis_b), (apoptosis_c), 
              (apoptosis_d), (apoptosis_e), (apoptosis_f),
              (apoptosis_g),  (apoptosis_h))

Differentiation <-c((differentiation_a),  (differentiation_b),  (differentiation_c),
                    (differentiation_d),  (differentiation_e),  (differentiation_f),
                    (differentiation_g),  (differentiation_h))

Proliferation <-c( (proliferation_a),  (proliferation_b),  (proliferation_c),
                   (proliferation_d),  (proliferation_e),  (proliferation_f),
                   (proliferation_g),  (proliferation_h))    

Network_Score <-c( (network_score_a),  (network_score_b),  (network_score_c),
                   (network_score_d),  (network_score_e),  (network_score_f),
                   (network_score_g),  (network_score_h))    


network_values_fixed_df<- data.frame(Figure, FLT3, DNMT3A, NPM1, AKT, CEBPA, GSK3B, ARF, HOXA9,
                                     FBXW7,ERK,CDKN2A,STAT5A,SOX4,CCND1,MEIS1,MYC, ETV6,TP53,BCL2,
                                     Apoptosis, Differentiation, Proliferation, Network_Score)

#saves the dataframe as a txt file 
#write.table(network_values_fixed_df, "path to file",sep=","row.names=FALSE)

#now we are making the scatterplots comparing the newtwork scores to clinical data
#load necessary packages for cleaning up the tables, making scatterplots
library(dplyr)
library(stringr)
require(ggplot2)
library(tidyr)
library(ggsci)
library(ggrepel) 

#had to change the column name of FLT3-ITD
colnames(s5_df)[4] ="FLT3_ITD"

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

#Making all values that are NOT DNMT3A,FLT3, or NPM1 NA
symbol_common <- symbol_common %>%
  mutate_at(vars(col_num_1:col_num_70), ~ifelse(. %in% c("NPM1", "DNMT3A","FLT3"), ., NA))

#Recombining all the columns, and the NA values "fall out" 
symbol_common <- symbol_common %>%
  unite(combined_column, col_num_1:col_num_70, sep = ",", na.rm = TRUE)

#grouping the new symbols based on diagram in Fig 5
symbol_common <- symbol_common %>%
  mutate(group_value = 
           ifelse(combined_column == "NPM1", 2, 
                  ifelse(combined_column == "FLT3", 3,
                         ifelse(combined_column == "DNMT3A", 4, 
                                ifelse(combined_column == "FLT3,NPM1", 5,
                                       ifelse(combined_column == "FLT3,DNMT3A", 6,
                                              ifelse(combined_column == "DNMT3A,NPM1", 7,
                                                     ifelse(combined_column == "FLT3,DNMT3A,NPM1", 8,
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

common_wo_NA_df$labID <- as.numeric(common_wo_NA_df$labId)
common_wo_NA_df$`%.Blasts.in.BM` <- as.numeric(common_wo_NA_df$`%.Blasts.in.BM`)
common_wo_NA_df$`%.Blasts.in.PB` <- as.numeric(common_wo_NA_df$`%.Blasts.in.PB`)

#getting the averages for the group IDs
WT_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_1_id)$`%.Blasts.in.PB`, na.rm = TRUE)
WT_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_1_id)$`%.Blasts.in.BM`, na.rm = TRUE)

NPM1_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_2_id)$`%.Blasts.in.PB`, na.rm = TRUE)
NPM1_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_2_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_3_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_3_id)$`%.Blasts.in.BM`, na.rm = TRUE)

DNMT3A_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_4_id)$`%.Blasts.in.PB`, na.rm = TRUE)
DNMT3A_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_4_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_NPM1_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_5_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_NPM1_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_5_id)$`%.Blasts.in.BM`, na.rm = TRUE)

FLT3_DNMT3A_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_6_id)$`%.Blasts.in.PB`, na.rm = TRUE)
FLT3_DNMT3A_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_6_id)$`%.Blasts.in.BM`, na.rm = TRUE)

NPM1_DNMT3A_mean_PB <- mean(subset(common_wo_NA_df, labId %in% group_7_id)$`%.Blasts.in.PB`, na.rm = TRUE)
NPM1_DNMT3A_mean_BM <- mean(subset(common_wo_NA_df, labId %in% group_7_id)$`%.Blasts.in.BM`, na.rm = TRUE)

all_mutation_PB <- mean(subset(common_wo_NA_df, labId %in% group_8_id)$`%.Blasts.in.PB`, na.rm = TRUE)
all_mutation_BM <- mean(subset(common_wo_NA_df, labId %in% group_8_id)$`%.Blasts.in.BM`, na.rm = TRUE)

#the network_score_ comes from the other R Script
#making dataframes for graphs
#PB Blast vs Network Score
PB_blast_vs_network_score <- data.frame( 
  PB_BLAST = c(get("WT_mean_PB"), get("NPM1_mean_PB"), get("FLT3_mean_PB"), get("DNMT3A_mean_PB"), get("FLT3_NPM1_mean_PB"), get("FLT3_DNMT3A_mean_PB"), get("NPM1_DNMT3A_mean_PB"),get("all_mutation_PB")),
  Network_Score = c(get("network_score_a"), get("network_score_b"), get("network_score_c"), get("network_score_d"), get("network_score_e"), get("network_score_f"), get("network_score_g"),get("network_score_h")),
  row.names = c('WT', 'NPM1', 'FLT3', 'DNMT3A', 'FLT3_NPM1', 'FLT3_DNMT3A', 'NPM1_DNMT3A', 'FLT3_NPM1_DNMT3A'))

print(PB_blast_vs_network_score)

#making dataframes for graphs
#BM Blast vs Network Score
BM_blast_vs_network_score <- data.frame( 
  BM_BLAST = c(get("WT_mean_BM"), get("NPM1_mean_BM"), get("FLT3_mean_BM"), get("DNMT3A_mean_BM"), get("FLT3_NPM1_mean_BM"), get("FLT3_DNMT3A_mean_BM"), get("NPM1_DNMT3A_mean_PB"),get("all_mutation_PB")),
  Network_Score = c(get("network_score_a"), get("network_score_b"), get("network_score_c"), get("network_score_d"), get("network_score_e"), get("network_score_f"), get("network_score_g"),get("network_score_h")),
  row.names = c('WT', 'NPM1', 'FLT3', 'DNMT3A', 'FLT3_NPM1', 'FLT3_DNMT3A', 'NPM1_DNMT3A', 'FLT3_NPM1_DNMT3A'))

print(BM_blast_vs_network_score)

#plotting the PB Blast vs Network Score scatter plot
PB_network_scatterplot <- ggplot(data = PB_blast_vs_network_score, aes(x = PB_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(PB_blast_vs_network_score)), vjust = -1) +
  labs(title = "% PB Blast vs Network Phenotypes (DNMT3A OR)", x = "% PB Blast", y = "Network Phenotype Score") +
  theme_minimal()
PB_network_scatterplot + coord_cartesian(xlim = c(20,85), ylim = c(-4, 10))

BM_network_scatterplot <- ggplot(data = BM_blast_vs_network_score, aes(x = BM_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(BM_blast_vs_network_score)), vjust = -1) +
  labs(title = "% BM Blast vs Network Phenotypes (DNMT3A OR)", x = "% BM Blast", y = "Network Phenotype Score") +
  theme_minimal()
BM_network_scatterplot + coord_cartesian(xlim = c(20,85), ylim = c(-4, 10))

#doing a p test 
PB_diagram <- PB_blast_vs_network_score$PB_BLAST
Network_diagram <- PB_blast_vs_network_score$Network_Score
BM_diagram <- BM_blast_vs_network_score$BM_BLAST

t.test(PB_diagram, Network_diagram, paired = FALSE)
t.test(BM_diagram, Network_diagram, paired = FALSE)

#doing a pearson correlation
PB_network_pearson <- cor(PB_blast_vs_network_score$PB_BLAST, PB_blast_vs_network_score$Network_Score)
print(PB_network_pearson)

BM_network_pearson <- cor(BM_blast_vs_network_score$BM_BLAST, BM_blast_vs_network_score$Network_Score)
print(BM_network_pearson)
