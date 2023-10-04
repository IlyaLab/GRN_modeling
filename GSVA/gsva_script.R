# load necessary packages 
library(readxl)
library(GSVA)

# importing the excel sheet (CPM)
s8_original_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/s8_table.xlsx")
s8_df <- s8_original_df[,!names(s8_original_df) %in% c("Gene", "Symbol")]

#fix dataframe to convert to matrix
colnames(s8_df) <- paste0("s", 1:ncol(s8_df))
rownames(s8_df) <- paste0("g", 1:nrow(s8_df))
s8_matrix <- as.matrix(s8_df)

#creating 100 sample gene sets 
p = nrow(s8_df) #number of genes
n = ncol(s8_df) #number of samples
gs <- as.list(sample(10:100, size=100, replace=TRUE))

gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))

#getting the GSVA scores mx.diff=TRUE, so calculates scores dependent of each other
gsva_test <- gsva(s8_matrix, gs, mx.diff= TRUE, verbose = FALSE)
