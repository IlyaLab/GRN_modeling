### Created by Tazein Shah
### Last modified on 1/1/2024
### This code gets the gene counts (of every unique gene) in the s7_combined_df

#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(data.table)
library(writexl)

#import the data
s7_combined_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s7_data_combined.xlsx", col_names = c("Patient_Id", "Symbol"))
s5_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s5_table.xlsx")

colnames(s5_df)[1] ="Patient_Id"
lab_ID_common <- intersect(s5_df$Patient_Id,s7_combined_df$Patient_Id)
s5_common_df <- subset(s5_df, Patient_Id %in% c(lab_ID_common))

#making a new dataframe with Patient_Id, PB, BM, and Symbol
combined_df <- data.frame(Patient_Id=lab_ID_common,
                          PB_Blast = as.numeric(s5_common_df$`%.Blasts.in.PB`),
                          BM_Blast = as.numeric(s5_common_df$`%.Blasts.in.BM`),
                          Symbol = s7_combined_df$Symbol)

combined_df <- na.omit(combined_df)

#split the character vector into individual genes (have to do it twice?)
individual_genes <- unlist(strsplit(combined_df$Symbol, ", "))
individual_genes <- unlist(strsplit(individual_genes, ","))

#create a data table from the individual_genes vector
gene_counts_dt <- data.table(Gene = individual_genes)

#use count the occurrences of each gene
gene_counts_dt <- gene_counts_dt[, .N, by = Gene]

#sort the data table by count in descending order to get the most common genes
setorder(gene_counts_dt, -N)

#select the top 40 most common genes
top_40_genes <- gene_counts_dt[1:40]

#print the top 40 genes
print(top_40_genes)
top_40 <- top_40_genes$Gene

#write_xlsx(gene_counts_dt,"C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/gene_count.xlsx", col_names = TRUE, format_headers = TRUE)
