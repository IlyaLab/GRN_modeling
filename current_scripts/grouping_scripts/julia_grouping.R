### Created by Tazein Shah
### Last modified on 1/2/2024
### This script will take the julia_patient_scores excel sheet and group the patients based off of their mutation profile (for now only DNMT3A-NPM1-FLT3)

#load necessary packages
library(dplyr)
library(stringr)
require(ggplot2)
library(tidyr)
library(ggsci)
library(ggrepel) 
library(readxl)

julia_scores <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/julia_patient_scores.xlsx") #has all the scores
s5_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s5_table.xlsx") #has the BM and PB Blast % 

#renaming both patient_id columns into labId
colnames(julia_scores)[colnames(julia_scores) == "Patient_Id"] ="labId"

#to make things easier, remove all the NA or non-numeric values from s5_df
s5_df$`%.Blasts.in.BM` <- as.numeric(s5_df$`%.Blasts.in.BM`)
s5_df$`%.Blasts.in.PB` <- as.numeric(s5_df$`%.Blasts.in.PB`)
s5_df <- na.omit(s5_df)

#getting what's in common between both s5 and s7 labId's
lab_ID_common <- intersect(s5_df$labId,julia_scores$labId)

#making new df with only the labId that overlap
s5_common_df <- subset(s5_df, labId %in% c(lab_ID_common))
julia_common_df <- subset(julia_scores, labId %in% c(lab_ID_common))

#make a merged df of the s5_common_df and s5_common_df
common_merge_df <- julia_common_df %>% full_join(s5_common_df)

#with the combined symbol column, separate the symbols into separate columns
common_merge_df <- common_merge_df %>%
  separate_rows(Mutation_Profile, convert = TRUE) %>%
  group_by(labId) %>%
  mutate(col = row_number()) %>%
  spread(col, Mutation_Profile, sep = "_") %>%
  ungroup()

#this makes all the values that are NOT DNMT3A, NPM1, FLT3 NA
common_merge_df <- common_merge_df %>%
  mutate_at(vars(col_1:col_131), ~ifelse(. %in% c("NPM1", "DNMT3A","FLT3"), ., NA))

#recombining all the columns, and the NA values "fall out" 
common_merge_df <- common_merge_df %>%
  unite(genes, col_1:col_131, sep = ",", na.rm = TRUE)

#removing duplicates
for (i in 1:(nrow(common_merge_df))){
  mutation_profile <- as.character(common_merge_df[i,8]) #get the column as a character variable
  mutation_profile <- strsplit(mutation_profile,split = ',') #split the character into separate characters
  mutation_profile <- lapply(mutation_profile,unique) #remove duplicates 
  mutation_profile <- paste0(mutation_profile[[1]], collapse = ",") #puts it back together
  
  common_merge_df[i,8] <- mutation_profile #put it back into the dataframe
}

#getting the means for every mutation group 
means <- list()
values <- c("Mutation", "PB_Blast", "BM_Blast", "Proliferation", "Differentiation", "Apoptosis", "Final")

gene_mutations <- c("", "FLT3,DNMT3A,NPM1", 
                    "DNMT3A,NPM1","FLT3,DNMT3A",
                    "FLT3,NPM1","DNMT3A", 
                    "FLT3", "NPM1")

for (i in 1:length(gene_mutations)) {
  means[[values[1]]][i] <- gene_mutations[i]
  means[[values[2]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$`%.Blasts.in.PB`)
  means[[values[3]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$`%.Blasts.in.BM`)
  means[[values[4]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$Proliferation)
  means[[values[5]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$Differentation)
  means[[values[6]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$Apoptosis)
  means[[values[7]]][i] <- mean(subset(common_merge_df, genes == gene_mutations[i])$Final)
}

mean_names <- c("None", "All", "DNMT3A_NPM1", "FLT3_DNMT3A","FLT3_NPM1","DNMT3A","FLT3", "NPM1")

#making the dataframes and plotting them
#final score
PB_final <- data.frame(PB_Blast = means$PB_Blast,
                       Final_Score = means$Final)

BM_final <- data.frame(BM_Blast = means$BM_Blast,
                       Final_Score = means$Final)

#plotting the PB Blast vs Network Score scatter plot
PB_network_scatterplot <- ggplot(data = PB_final, aes(x = PB_Blast, y = Final_Score)) +
  geom_point() +
  geom_text_repel(aes(label = mean_names), vjust = -1) +
  labs(title = "% PB Blast vs Final Score Julia: Using DNMT3A-NPM1-FLT3 Grouping", x = "% PB Blast", y = "Final Score Julia") +
  theme_minimal()
PB_network_scatterplot 

BM_network_scatterplot <- ggplot(data = BM_final, aes(x = BM_Blast, y = Final_Score)) +
  geom_point() +
  geom_text_repel(aes(label = mean_names), vjust = -1) +
  labs(title = "% BM Blast vs Final Score Julia: Using DNMT3A-NPM1-FLT3 Grouping", x = "% BM Blast", y = "Final Score Julia") +
  theme_minimal()
BM_network_scatterplot 

#doing a p test 
t.test(means$PB_Blast, means$Final, paired = FALSE)
t.test(means$BM_Blast, means$Final, paired = FALSE)

#doing a pearson correlation
PB_pearson <- cor(means$PB_Blast, means$Final)
print(PB_pearson)

BM__pearson <- cor(means$BM_Blast, means$Final)
print(BM__pearson)

#creating subsets with the mutations to see the score distribution within each mutation group
scores <- list()

for (i in 1:length(mean_names)) {
  scores[[mean_names[i]]] <- subset(common_merge_df, genes == gene_mutations[i])
}

#making the box plots

none_boxplot <- boxplot(scores$None[5],
                        main = "Distribution of the Final Scores - None",
                        xlab = "Patients",
                        ylab = "Scores",
                        col = "red",
                        border = "black",
                        horizontal = FALSE,
                        notch = TRUE)

all_boxplot <- boxplot(scores$All[5],
                        main = "Distribution of the Final Scores - All",
                        xlab = "Patients",
                        ylab = "Scores",
                        col = "orange",
                        border = "black",
                        horizontal = FALSE,
                        notch = TRUE)

DNMT3A_NPM1_boxplot <- boxplot(scores$DNMT3A_NPM1[5],
                       main = "Distribution of the Final Scores - DNMT3A_NPM1",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "yellow",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)

FLT3_DNMT3A_boxplot <- boxplot(scores$FLT3_DNMT3A[5],
                       main = "Distribution of the Final Scores - FLT3_DNMT3A",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "green",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)

FLT3_NPM1_boxplot <- boxplot(scores$FLT3_NPM1[5],
                       main = "Distribution of the Final Scores - FLT3_NPM1",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "blue",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)

DNMT3A_boxplot <- boxplot(scores$DNMT3A[5],
                       main = "Distribution of the Final Scores - DNMT3A",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "purple",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)

FLT3_boxplot <- boxplot(scores$FLT3[5],
                       main = "Distribution of the Final Scores - FLT3",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "pink",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)

NPM1_boxplot <- boxplot(scores$NPM1[5],
                       main = "Distribution of the Final Scores - NPM1",
                       xlab = "Patients",
                       ylab = "Scores",
                       col = "gray",
                       border = "black",
                       horizontal = FALSE,
                       notch = TRUE)
  
  
