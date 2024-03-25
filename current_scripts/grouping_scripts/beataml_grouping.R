### Created by Tazein on 1/8/24
### Last modified on 1/9/2024
### This script with take the BeatAML data, create the Wilcoxon rank sum test which will show correlation of specific genes with PB and BM blast %, and based on those genes, create mutation groups of patients and show the PB and BM blast % distribution of those genes in boxplots

#load necessary packages
library(dplyr)
library(stringr)
require(ggplot2)
library(tidyr)
library(ggsci)
library(ggrepel) 
library(readxl)
library(writexl)
library(data.table)

#import the data
s7_combined_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s7_data_combined.xlsx", col_names = c("Patient_Id", "Symbol"))
s5_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s5_table.xlsx")

colnames(s5_df)[1] ="Patient_Id"
lab_ID_common <- intersect(s5_df$Patient_Id,s7_combined_df$Patient_Id)
s5_common_df <- subset(s5_df, Patient_Id %in% c(lab_ID_common))

#### getting the top 100 genes and doing the rank sum test ###
#making a new dataframe with Patient_Id, PB, BM, and Symbol
combined_df <- data.frame(Patient_Id=lab_ID_common,
                          PB_Blast = as.numeric(s5_common_df$`%.Blasts.in.PB`),
                          BM_Blast = as.numeric(s5_common_df$`%.Blasts.in.BM`),
                          Symbol = s7_combined_df$Symbol)

combined_df <- na.omit(combined_df)

#split the character vector into individual genes (have to do it twice?)
individual_genes <- unlist(strsplit(combined_df$Symbol, ", "))
individual_genes <- unlist(strsplit(individual_genes, ","))

gene_counts_dt <- data.table(Gene = individual_genes)
gene_counts_dt <- gene_counts_dt[, .N, by = Gene]

#sort the data table by count in descending order to get the most common genes
setorder(gene_counts_dt, -N)

#select the top 100 most common genes
top_genes <- gene_counts_dt[1:100]
print(top_genes)
top_genes <- top_genes$Gene

#with combined_df, split the symbol column up into separate columns
combined_df <- combined_df %>%
  separate_rows(Symbol, convert = TRUE) %>%
  group_by(Patient_Id) %>%
  mutate(col = row_number()) %>%
  spread(col, Symbol, sep = "_") %>%
  ungroup()

#making the rest file for the loop
reset_df <- combined_df

gene_name <- NA
pval_PB <- NA
pval_BM <- NA
wilcox_df <- data.frame(gene_name,
                        pval_PB,
                        pval_BM)

#we are going to loop through the top_genes
for(gene in top_genes){
  print(gene)
  
  #making all values that are NOT the gene we want to compare NA
  combined_df <- combined_df %>%
    mutate(across(starts_with("col_"), ~ifelse(. %in% gene, ., NA)))
  
  #Recombining all the columns, and the NA values "fall out"
  combined_df <- combined_df %>%
    unite(combined_column, starts_with("col_"), sep = ",", na.rm = TRUE)
  
  #grouping into two groups: mutated and wild_type
  combined_df <- combined_df %>%
    mutate(group_value = ifelse(combined_column == "", "wild_type","mutated"))
  
  print("seperating into groups")
  
  combined_df <- na.omit(combined_df)
  
  #making a vector with the PB values of mutated and wild_type
  WT_PB <- as.numeric(na.omit(combined_df$PB_Blast[combined_df$group_value == "wild_type"]))
  M_PB <- as.numeric(na.omit(combined_df$PB_Blast[combined_df$group_value == "mutated"]))
  
  WT_BM <- as.numeric(na.omit(combined_df$BM_Blast[combined_df$group_value == "wild_type"]))
  M_BM <- as.numeric(na.omit(combined_df$BM_Blast[combined_df$group_value == "mutated"]))
  
  print("finish making vectors")
  
  #doing the wilcox test
  PB_wilcox = wilcox.test(WT_PB, M_PB, alternative = "two.sided")
  BM_wilcox = wilcox.test(WT_BM, M_BM, alternative = "two.sided")
  
  wilcox_df[(which(top_genes == gene)),2] <- rbind(PB_wilcox$p.value)
  wilcox_df[(which(top_genes == gene)),3] <- rbind(BM_wilcox$p.value)
  
  print("added to dataframe")
  
  combined_df <-reset_df
  
}

wilcox_df$gene_name <- top_genes

#adding the pval analysis as well (so I can remember)
wilcox_df$PB_significant <- ifelse(wilcox_df$pval_PB < 0.05, "YES", "NO")
wilcox_df$BM_significant <- ifelse(wilcox_df$pval_BM < 0.05, "YES", "NO")


### now that we have the Wilcoxon rank sum test results, we need to make groups ###
#getting the genes from the wilcoxon test 
PB_genes <- (subset(wilcox_df, PB_significant == 'YES')$gene_name)
BM_genes <- (subset(wilcox_df, BM_significant == 'YES')$gene_name)

combined_genes <- union(PB_genes, BM_genes)

#this makes all the values that are NOT in the PB_genes and BM_genes NA
combined_df <- combined_df %>%
  mutate_at(vars(col_1:col_131), ~ifelse(. %in% c(combined_genes), ., NA))

#recombining all the columns, and the NA values "fall out" 
combined_df <- combined_df %>%
  unite(genes, col_1:col_131, sep = ",", na.rm = TRUE)

#test_df <- combined_df
#combined_df <- test_df

#removing duplicates
for (i in 1:(nrow(combined_df))){
  mutation_profile <- as.character(combined_df[i,4]) #get the column as a character variable
  mutation_profile <- strsplit(mutation_profile,split = ',') #split the character into separate characters
  mutation_profile <- lapply(mutation_profile,unique) #remove duplicates 
  mutation_profile <- paste0(mutation_profile[[1]], collapse = ",") #puts it back together
  
  combined_df[i,4] <- mutation_profile #put it back into the dataframe
}

#replacing any NA mutation group with 'None'
combined_df$genes <- gsub("^$", "None", combined_df$genes)

#getting all the unique mutation profiles 
mutation_groups <- as.vector(unique(combined_df$genes))

#getting the PB and BM Blast % for each mutation groups
PB_Blast <- list()
BM_Blast <- list()

for (i in 1:length(mutation_groups)) {
  PB_Blast[[mutation_groups[i]]] <- (subset(combined_df, genes == mutation_groups[i])$PB_Blast)
  BM_Blast[[mutation_groups[i]]] <- (subset(combined_df, genes == mutation_groups[i])$BM_Blast)
  
}

#making boxplot plots for the mutation_groups

#removing any groups with less than three values
PB_Blast <- Filter(function(gene_values) length(gene_values) > 3, PB_Blast)
BM_Blast <- Filter(function(gene_values) length(gene_values) > 3, BM_Blast)

#names(PB_Blast) <- c("NPM1","FLT3","FLT3_NPM1","None","FLT3_IDH1_NPM1","NRAS","IKZF1","STAG2","TP53","NRAS_NPM1","IDH1","NRAS_TP53")

PB_boxplot <- boxplot(PB_Blast,
                        las = 2,
                        par(mar = c(11, 5, 3, 2)+ 0.05),
                        main = "PB Blast for Mutation Groups",
                        xlab = "",
                        ylab = "PB Blast %",
                        col = "royalblue2",
                        border = "black",
                        horizontal = FALSE,
                        notch = FALSE)


BM_boxplot <- boxplot(BM_Blast,
                      las = 2,
                      par(mar = c(11, 5, 3, 2)+ 0.05),
                      main = "BM Blast for Mutation Groups",
                      xlab = "",
                      ylab = "BM Blast %",
                      col = "orange",
                      border = "black",
                      horizontal = FALSE,
                      notch = FALSE)
