### Created by Tazein Shah 
### Last modified on 10/11/2023
### This code takes the s8 table from the clinical data and runs a GSVA (gene set variation analysis), then it compares the results to the PB and BM Blast from s5 

#load necessary packages 
library(readxl)
library(GSVA)
library(ggplot2)

#importing the excel sheet (CPM)
s8_original_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s8_table.xlsx")
s8_df <- s8_original_df[,!names(s8_original_df) %in% c("Gene", "Symbol")]

#fix dataframe to convert to matrix
#colnames(s8_df) <- paste0("s", 1:ncol(s8_df)) #not using bc we want to see the patient ID
rownames(s8_df) <- paste0("g", 1:nrow(s8_df))
s8_matrix <- as.matrix(s8_df)

#finding out the corresponding gene like (FLT3) to the matrix (g5)
#which(s8_original_df$Symbol == "DPM1")
#!is.na(match(prolif_up, s8_original_df$Symbol))
#paste0("g", which(s8_original_df$Symbol %in% apop_up), collapse = ",")

#signor gsva: using the values from the signor data
prolif_up <- c("g5073", "g10378", "g6394","g5815","g3110","g14242", "g7345", "g3318","g5899", "g10484", "g3110", "g4875", "g11852", "g4574") #("STAT5A","CTNNB1","MYC","CCNA1","MTOR","MEIS1","CCND1","MEN1","STAT3","UBTF","SOX4", "JUN", "PIK3CA") #"PML-RARalpha"
prolif_down <- c("g19","g3287","g4261","g7017","g7754","g7994","g10507","g12602","g12786","g12855","g14475") #("TP53","CDKN2A","WT1","TET2","BCOR","ZBTB16","BAD", "FOXO1", "FOXO3", "FOXO4","FOXO6") 
diff_up <- c("g871","g2048","g4373","g6715","g7803","g8811","g19209") #("SPI1","DNMT3A","EP300","RUNX1","CEBPA","ETV6","NOTCH1")
diff_down <-c("g1285","g4875","g7345","g12602")#("SOX4","HOXA9","MEIS1","BCOR") #"CBFbeta-MYH11"
apop_up <- c("g19","g430","g1538","g3287","g4261","g7017","g7994","g12786","g14475") #("TP53","BAK1","BAX","BAD","ZBTB16", "FOXO1", "FOXO3", "FOXO4","FOXO6")
apop_down <- c("g5073","g6424","g7314","g7345","g10978","g11017") #("BCL2","MEIS1","STAT5A","BCL2L1","PARP1", "PIM1")

#figure 3 network gsva: using the values from figure 3 network
prolif_up_network <- c("g3110","g3318","g4574","g5073","g5815","g6394","g7345","g10378","g11852") #("STAT5A","CTNNB1","MYC","CCNA1","MEIS1","CCND1","UBTF","JUN", "PIK3CA")
prolif_down_network <- c("g4261","g7754","g7994","g12602","g12786","g12855","g14475") #("CDKN2A","WT1","BCOR","FOXO1", "FOXO3", "FOXO4","FOXO6") 
diff_up_network <- c("g871","g2048","g6715","g8811","g19209") #("SPI1","EP300","RUNX1","CEBPA","ETV6")
diff_down_network <-c("g4875") #("SOX4")
apop_up_network <- c("g1538","g4261","g7017","g7994","g12786","g14475") #("TP53","BAX", "FOXO1", "FOXO3", "FOXO4","FOXO6")
apop_down_network <- c("g6424","g11017") # ("BCL2", "PIM1")

#combine the character vectors into a list
gs <- list(prolif_up=prolif_up,
           prolif_down=prolif_down,
           diff_up=diff_up,
           diff_down=diff_down,
           apop_up=apop_up,
           apop_down=apop_down)

#combine the character vectors into a list
gs_network <- list(prolif_up_network=prolif_up_network,
                   prolif_down_network=prolif_down_network,
                   diff_up_network=diff_up_network,
                   diff_down_network=diff_down_network,
                   apop_up_network=apop_up_network,
                   apop_down_network=apop_down_network)

#getting the GSVA scores mx.diff=TRUE, so calculates scores dependent of each other
gsva_test <- gsva(s8_matrix, gs, mx.diff = TRUE, verbose = TRUE)
gsva_network <- gsva(s8_matrix, gs_network, mx.diff = TRUE, verbose = TRUE)

###working with the gsva results ###
#there are two results: the gsva_test and the gsva_network, which have different genes

#tranposing the results and then turning it into a dataframe
data <- (data.frame(t(gsva_test)))
data_network <- (data.frame(t(gsva_network)))

#calculating proliferation, differentiation, and apoptosis
data$prolif_total <- (data$prolif_up - data$prolif_down)
data$diff_total <- (data$diff_up - data$diff_down)
data$apop_total <- (data$apop_up - data$apop_down)
data$network_score <- data$prolif_total - ((data$diff_total) + (data$apop_total))

data_network$prolif_total <- (data_network$prolif_up_network - data_network$prolif_down_network)
data_network$diff_total <- (data_network$diff_up_network - data_network$diff_down_network)
data_network$apop_total <- (data_network$apop_up_network - data_network$apop_down_network)
data_network$network_score <- data_network$prolif_total - ((data_network$diff_total) + (data_network$apop_total))

#importing the clinical data from table 5
s5_df <- read_excel("C:/Users/15167/OneDrive/Documents/ISB/AML-DT-BNM/raw_data/s5_table.xlsx")

#adding the data from table 5 to the gsva results
samples <- rownames(data)
data$BM_Blast <- (subset(s5_df, labId %in% samples)$`%.Blasts.in.BM`)
data$PB_Blast <- (subset(s5_df, labId %in% samples)$`%.Blasts.in.PB`)

samples <- rownames(data_network)
data_network$BM_Blast <- (subset(s5_df, labId %in% samples)$`%.Blasts.in.BM`)
data_network$PB_Blast <- (subset(s5_df, labId %in% samples)$`%.Blasts.in.PB`)

#turning into numerics by removing characters
data$BM_Blast <- (gsub('[^0-9.]', NA, data$BM_Blast))
data$PB_Blast <- (gsub('[^0-9.]', NA, data$PB_Blast))

data_network$BM_Blast <- (gsub('[^0-9.]', NA, data_network$BM_Blast))
data_network$PB_Blast <- (gsub('[^0-9.]', NA, data_network$PB_Blast))

#testing with NA removed to see if that helps correlation
data <- na.omit(data)
data$BM_Blast <- as.numeric(data$BM_Blast)
data$PB_Blast <- as.numeric(data$PB_Blast)

data_network <- na.omit(data)
data_network$BM_Blast <- as.numeric(data_network$BM_Blast)
data_network$PB_Blast <- as.numeric(data_network$PB_Blast)

#making plots to observe the correlation
#making the scatterplots
PB_scatterplot <- ggplot(data = data, aes(x = PB_Blast, y = network_score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "GSVA scores vs PB Blast", x = "% PB Blast", y = "GSVA Network Score") +
  theme_minimal()
PB_scatterplot 

PB_pearson <- cor(data$PB_Blast, data$network_score)
print(PB_pearson)

BM_scatterplot <- ggplot(data = data, aes(x = BM_Blast, y = network_score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "GSVA scores vs BM Blast", x = "% BM Blast", y = "GSVA Network Score") +
  theme_minimal()
BM_scatterplot 

BM_pearson <- cor(data$BM_Blast, data$network_score)
print(BM_pearson)


#now we make the network (fig 3) scatterplots
PB_network_scatterplot <- ggplot(data = data_network, aes(x = PB_Blast, y = network_score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "GSVA scores Fig 3 vs PB Blast", x = "% PB Blast", y = "GSVA Network Score") +
  theme_minimal()
PB_network_scatterplot 

PB_pearson <- cor(data_network$PB_Blast, data_network$network_score)
print(PB_pearson)

BM_network_scatterplot <- ggplot(data = data_network, aes(x = BM_Blast, y = network_score)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "GSVA scores Fig 3 vs BM Blast", x = "% BM Blast", y = "GSVA Network Score") +
  theme_minimal()
BM_network_scatterplot 

BM_pearson <- cor(data_network$BM_Blast, data_network$network_score)
print(BM_pearson)
