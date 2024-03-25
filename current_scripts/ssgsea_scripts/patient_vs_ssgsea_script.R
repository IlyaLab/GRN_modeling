### Created by Tazein Shah (using the ssgsea scores that Guangrong Qin gave to me) 
### Last modified on 10/22/2023
### This code compares the julia patient scores to the ssgsea scores that Guangrong gave me 

#loading the packages
library(readxl)
library(ggplot2)

#import the data
julia_scores <- read_excel("/path/julia_patient_scores.xlsx")
ssGSEA_scores <- read_excel("/path/matrix_factor_beatAML.xlsx")
colnames(ssGSEA_scores)[1] ="Patient_Id"

#they don't have the same number of samples, so I need to get the samples they have in common
#julia score has 608, ssGSEA has 451
lab_ID_common <- intersect(julia_scores$Patient_Id,ssGSEA_scores$Patient_Id)
julia_scores <- subset(julia_scores, Patient_Id %in% c(lab_ID_common))
ssGSEA_scores <- subset(ssGSEA_scores, Patient_Id %in% c(lab_ID_common))
#now they both have 399 samples

### Proliferation ###
#creating a new dataframe with Patient_Id, julia_prolif, and ssGSEA_prolif
proliferation_df <- data.frame(Patient_Id =julia_scores$Patient_Id, 
                               julia_prolif = julia_scores$Proliferation, 
                               ssGSEA_prolif =ssGSEA_scores$Proliferation_ssGSEA)

#comparing the patient proliferation scores to proliferation ssGSEA scores
prolif_scatterplot <- ggplot(data = proliferation_df, aes(x = julia_prolif, y = ssGSEA_prolif)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Comparing the julia proliferation scores to ssGSEA proliferation", x = "Julia Proliferation Score", y = "ssGSEA Proliferation") +
  theme_minimal()
prolif_scatterplot 

prolif_pearson <- cor(proliferation_df$julia_prolif, proliferation_df$ssGSEA_prolif)
print(prolif_pearson)

### Apoptosis ###
#creating a new dataframe with Patient_Id, julia_prolif, and ssGSEA_prolif
apoptosis_df <- data.frame(Patient_Id =julia_scores$Patient_Id, 
                               julia_apoptosis = julia_scores$Apoptosis, 
                               ssGSEA_apoptosis =ssGSEA_scores$Apoptosis_ssGSEA)

#comparing the patient proliferation scores to proliferation ssGSEA scores
apop_scatterplot <- ggplot(data = apoptosis_df, aes(x = julia_apoptosis, y = ssGSEA_apoptosis)) +
  geom_point() +
  #geom_text_repel(aes(label = id), vjust = -1) +
  labs(title = "Comparing the julia apoptosis scores to ssGSEA apoptosis", x = "Julia Apoptosis Score", y = "ssGSEA Apoptosis") +
  theme_minimal()
apop_scatterplot 

apop_pearson <- cor(apoptosis_df$julia_apoptosis, apoptosis_df$ssGSEA_apoptosis)
print(apop_pearson)
