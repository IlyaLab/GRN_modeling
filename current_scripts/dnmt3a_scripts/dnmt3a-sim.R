### Created by Tazein Shah
### Last modified on 2/20/2024
### This code runs the FLT3-NPM1-DNMT3A subgroup from the Palma paper, and runs the 8 groups from the fig 5, it then compares the scores from these eight profiles with the PB and BM Blast %

#load packages
library(BoolNet)
library(dplyr)
library(stringr)
require(ggplot2)
library(tidyr)
library(ggsci)
library(ggrepel) 
library(readxl)

#import the s5 and s7 dataframes
s7_combined_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s7_data_combined.xlsx", col_names = c('labId','symbol'))
s5_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s5_table.xlsx")

#changing the colnames of s5_df and making them numeric
colnames(s5_df)[2] = 'BM_Blast'
colnames(s5_df)[3] = 'PB_Blast'

s5_df$BM_Blast <- as.numeric(s5_df$BM_Blast)
s5_df$PB_Blast <- as.numeric(s5_df$PB_Blast)

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
cat("MEIS1,!(DNMT3A & !HOXA9)\n")
cat("MYC, !(GSK3B & FBXW7) & ERK\n")
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


#####getting the proliferation, apoptosis, and differentiation##### 

#define the lists to store the values
table_list <- list()  #list to store the tables
score_list <- list()  #list to store the network phenotype scores

for (i in 1:8) {
  matrix <- switch(i,
                   matrix_a,
                   matrix_b,
                   matrix_c,
                   matrix_d,
                   matrix_e,
                   matrix_f,
                   matrix_g,
                   matrix_h)
  
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
  
  #save the apoptosis as a separate object
  score_name <- paste0("apoptosis_", letters[i])
  score_list[[score_name]] <- apoptosis
  assign(score_name, score_list[[score_name]])
  
  #save the proliferation as a separate object
  score_name <- paste0("proliferation_", letters[i])
  score_list[[score_name]] <- proliferation
  assign(score_name, score_list[[score_name]])
  
  #save the differentiation as a separate object
  score_name <- paste0("differentiation_", letters[i])
  score_list[[score_name]] <- differentiation
  assign(score_name, score_list[[score_name]])
  
  #save the network phenotype score as a separate object
  score_name <- paste0("network_score_", letters[i])
  score_list[[score_name]] <- network_phenotype_score
  assign(score_name, score_list[[score_name]])
}


#######making the mega table of values#######

genes <- c("FLT3", "DNMT3A", "NPM1", "AKT", "CEBPA", "GSK3B", "ARF", "HOXA9", "FBXW7", "ERK",
           "CDKN2A", "STAT5A", "SOX4", "CCND1", "MEIS1", "MYC", "ETV6", "TP53", "BCL2")

#initialize empty lists to store gene values
gene_values <- list()

for (gene in genes) {
  #extract values for each gene from matrices
  gene_values[[gene]] <- c(matrix_a[gene, "Attr1.1"], matrix_b[gene, "Attr1.1"],
                           matrix_c[gene, "Attr1.1"], matrix_d[gene, "Attr1.1"],
                           matrix_e[gene, "Attr1.1"], matrix_f[gene, "Attr1.1"],
                           matrix_g[gene, "Attr1.1"], matrix_h[gene, "Attr1.1"])
}

network_values_fixed_df <- data.frame(Figure = c("figure_a","figure_b","figure_c",
                                                 "figure_d","figure_e","figure_f",
                                                 "figure_g","figure_h"),
                                      do.call(cbind, gene_values),
                                      Apoptosis = c(apoptosis_a, apoptosis_b, apoptosis_c,
                                                    apoptosis_d, apoptosis_e, apoptosis_f,
                                                    apoptosis_g, apoptosis_h),
                                      Differentiation = c(differentiation_a, differentiation_b, differentiation_c,
                                                          differentiation_d, differentiation_e, differentiation_f,
                                                          differentiation_g, differentiation_h),
                                      Proliferation = c(proliferation_a, proliferation_b, proliferation_c,
                                                        proliferation_d, proliferation_e, proliferation_f,
                                                        proliferation_g, proliferation_h),
                                      Network_Score = c(network_score_a, network_score_b, network_score_c,
                                                        network_score_d, network_score_e, network_score_f,
                                                        network_score_g, network_score_h))

#saves the dataframe as a txt file 
#write.table(network_values_fixed_df, "path to file",sep=","row.names=FALSE)

#######making the dataframes for the plots####### 
#combines the symbols rows into symbol column
s7_combined_df <- s7_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#merging together the s7_combined_df and s5_df
common <- intersect(s5_df$labId,s7_combined_df$labId)

s5_merged_df <- subset(s5_df, labId %in% c(common))
s7_merged_df <- subset(s7_combined_df, labId %in% c(common))

merged_df <- s7_merged_df %>% full_join(s5_merged_df)

#replacing the NA values with the mean
merged_df <- merged_df %>%
  mutate(PB_Blast = ifelse(is.na(PB_Blast), mean(PB_Blast, na.rm = TRUE), PB_Blast),
         BM_Blast = ifelse(is.na(BM_Blast), mean(BM_Blast, na.rm = TRUE), BM_Blast))

#with the combined symbol column, separate the symbols into separate columns
merged_df <- merged_df %>%
  separate_rows(symbol, convert = TRUE) %>%
  group_by(labId) %>%
  mutate(col_num = row_number()) %>%
  spread(col_num, symbol, sep = "_") %>%
  ungroup()

#making all values that are NOT DNMT3A,FLT3, or NPM1 NA
merged_df <- merged_df %>%
  mutate_at(vars(starts_with("col")), ~ifelse(. %in% c("NPM1", "DNMT3A", "FLT3"), ., NA))

#recombining all the columns, and the NA values "fall out" 
merged_df <- merged_df %>%
  unite(combined_column, starts_with("col"), sep = ",", na.rm = TRUE)

#removing duplicate values in merged_df combined column
for (i in seq_along(merged_df$combined_column)) {
  merged_df$combined_column[i] <- toString(unique(unlist(strsplit(merged_df$combined_column[i], ","))))
}

#grouping the new symbols based on diagram in Fig 5
merged_df <- merged_df %>%
  mutate(group_value = 
           ifelse(combined_column == "NPM1", 2, 
                  ifelse(combined_column == "FLT3", 3,
                         ifelse(combined_column == "DNMT3A", 4, 
                                ifelse(combined_column == "FLT3, NPM1", 5,
                                       ifelse(combined_column == "FLT3, DNMT3A", 6,
                                              ifelse(combined_column == "DNMT3A, NPM1", 7,
                                                     ifelse(combined_column == "FLT3, DNMT3A, NPM1", 8,
                                                            ifelse(combined_column == "", 1,NA)))))))))

#saving the values into 2 dataframes for the plots

mutation_names <- c("WT", "NPM1", "FLT3", "DNMT3A", "FLT3_NPM1",
                    "FLT3_DNMT3A", "NPM1_DNMT3A", "all_mutation")

PB_vs_network <- data.frame(PB_BLAST = numeric(0), Network_Score = numeric(0))
BM_vs_network <- data.frame(BM_BLAST = numeric(0), Network_Score = numeric(0))

for (i in 1:8) {
  #filter merged_df for the current group value and pull labIds
  group_id <- merged_df %>%
    filter(group_value == i) %>%
    pull(labId)
  
  #calculate mean for the current group ID
  PB_mean <- mean(subset(merged_df, labId %in% group_id)$PB_Blast, na.rm = TRUE)
  BM_mean <- mean(subset(merged_df, labId %in% group_id)$BM_Blast, na.rm = TRUE)
  
  #get the network score (from previous section)
  network_score <- get(paste0("network_score_", letters[i]))
  
  #append the values to the dataframes
  PB_vs_network <- rbind(PB_vs_network, data.frame(PB_BLAST = PB_mean, Network_Score = network_score))
  BM_vs_network <- rbind(BM_vs_network, data.frame(BM_BLAST = BM_mean, Network_Score = network_score))
}

#add rownames
rownames(PB_vs_network) <- mutation_names
rownames(BM_vs_network) <- mutation_names

print(PB_vs_network)
print(BM_vs_network)

#######making the scatterplots to compare the scores and clinical data####### 
#plotting the PB Blast vs Network Score scatter plot
PB_network_scatterplot <- ggplot(data = PB_vs_network, aes(x = PB_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(PB_vs_network)), vjust = -1) +
  labs(title = "% PB Blast vs Network Phenotypes", x = "% PB Blast", y = "Network Phenotype Score") +
  theme_minimal()
PB_network_scatterplot + coord_cartesian(xlim = c(20,85), ylim = c(-4, 10))

BM_network_scatterplot <- ggplot(data = BM_vs_network, aes(x = BM_BLAST, y = Network_Score)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(BM_vs_network)), vjust = -1) +
  labs(title = "% BM Blast vs Network Phenotypes", x = "% BM Blast", y = "Network Phenotype Score") +
  theme_minimal()
BM_network_scatterplot + coord_cartesian(xlim = c(20,85), ylim = c(-4, 10))

#doing a p test 
PB_diagram <- PB_vs_network$PB_BLAST
Network_diagram <- PB_vs_network$Network_Score
BM_diagram <- BM_vs_network$BM_BLAST

t.test(PB_diagram, Network_diagram, paired = FALSE)
t.test(BM_diagram, Network_diagram, paired = FALSE)

#doing a pearson correlation
PB_network_pearson <- cor(PB_vs_network$PB_BLAST, PB_vs_network$Network_Score)
print(PB_network_pearson)

BM_network_pearson <- cor(BM_vs_network$BM_BLAST, BM_vs_network$Network_Score)
print(BM_network_pearson)
