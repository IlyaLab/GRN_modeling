# load necessary packages 
library(readxl)
library(GSVA)
library(ggplot2)
library(dplyr)

# importing the excel sheet (CPM)
s8_original_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s8_table.xlsx")
s8_df <- s8_original_df[,!names(s8_original_df) %in% c("Gene", "Symbol")]

#fix dataframe to convert to matrix
colnames(s8_df) <- paste0("s", 1:ncol(s8_df))
rownames(s8_df) <- paste0("g", 1:nrow(s8_df))
s8_matrix <- as.matrix(s8_df)

#finding out the corresponding gene like (FLT3) to the matrix (g5)
#which(s8_original_df$Symbol == "DPM1")

#create three gene sets as character vectors (based off fig 3 equations)
proliferation <- c("g5073","g10378","g6394","g5815","g14242","g7345","g3318",
                   "g8201","g5899","g10484","g3110","g4875", "g4574", "g11852") 

differentiation <- c("g871","g4373","g2048","g8811","g19209","g6715","g7803")

apoptosis <- c("g7017","g430","g1538","g19","g3287", "g7994", "g4261", "g12786", "g14475") #used all the FOXO family 

#combine the character vectors into a list
gs <- list(
  gs_prolif = proliferation,
  gs_diff = differentiation,
  gs_apop = apoptosis
)

#getting the GSVA scores mx.diff=TRUE, so calculates scores dependent of each other
gsva_test <- gsva(s8_matrix, gs, mx.diff = TRUE, verbose = TRUE)

#tranposing the resukts and then turning it into a dataframe
data <- (data.frame(t(gsva_test)))

#calculating the network_score for each patient
data$network_score <- data$gs_prolif - ((data$gs_diff) + (data$gs_apop))
