### Created by Tazein Shah
### Last modified on 10/23/2023
### This script does the Wilcox rank sum test b/w two groups of patients: mutant and wild type group. It then compares the blast % of those two groups. 

#load packages
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
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

#with combined_df, split the symbol column up into separate columns
combined_df <- combined_df %>%
  separate_rows(Symbol, convert = TRUE) %>%
  group_by(Patient_Id) %>%
  mutate(col_num = row_number()) %>%
  spread(col_num, Symbol, sep = "_") %>%
  ungroup()

#making the rest file for the loop
reset_df <- combined_df

#making a vector with all the genes included in the figure 3 network (for now, can expand on this later)
top_40_genes <- c("FLT3","DNMT3A","IDH2", "TET2","NPM1",
                "NRAS","RUNX1", "TP53", "SRSF2", "IDH1",
                "WT1", "STAG2", "PTPN11", "CEBPA","SF3B1",
                "ASXL1", "U2AF1","BCOR","KRAS","GATA2","JAK2",
                "PHF6", "NF1", "ACAN", "RAD21", "CSF3R", "SMC1A",
                "IKZF1", "SMC3", "SETBP1", "KIT", "JAK3","PUF60",
                "BCORL1", "HUNK", "OTOP2", "SUZ12","PDS5B","EZH2","COL19A1")

#top_40_genes <- c("DNMT3A", "FLT3", "NPM1") #this is the test, which is going to rewrite the previous vector
#gene <- "DNMT3A"

gene_name <- NA
pval_PB <- NA
pval_BM <- NA
wilcox_df <- data.frame(gene_name,
                        pval_PB,
                        pval_BM)

#we are going to loop through the top_40_genes
for(gene in top_40_genes){
  print(gene)
  
  #making all values that are NOT the gene we want to compare NA
  combined_df <- combined_df %>%
    mutate(across(starts_with("col_num_"), ~ifelse(. %in% gene, ., NA)))
  
  #Recombining all the columns, and the NA values "fall out"
  combined_df <- combined_df %>%
    unite(combined_column, starts_with("col_num_"), sep = ",", na.rm = TRUE)
  
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
  
  wilcox_df[(which(top_40_genes == gene)),2] <- rbind(PB_wilcox$p.value)
  wilcox_df[(which(top_40_genes == gene)),3] <- rbind(BM_wilcox$p.value)
  
  print("added to dataframe")
  
  combined_df <-reset_df
  
}

wilcox_df$gene_name <- top_40_genes

#adding the pval analysis as well (so I can remember)
wilcox_df$PB_significant <- ifelse(wilcox_df$pval_PB < 0.05, "YES", "NO")
wilcox_df$BM_significant <- ifelse(wilcox_df$pval_BM < 0.05, "YES", "NO")

write_xlsx(wilcox_df,"C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/wilcox_pvalue.xlsx", col_names = TRUE, format_headers = TRUE)

