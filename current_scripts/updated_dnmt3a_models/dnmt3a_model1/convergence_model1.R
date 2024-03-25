### Created by Tazein Shah (uses FLT3_NPM1_DNMT3A_model1.xlsx from Guangrong Qin)
### Last modified on 11/12/2023
### Using the FLT3_NPM1_DNMT3A_model1 from Guangrong, this code does a checks the convergence of the time series (before simulation)

#testing the convergence of the dnmt3a_model1 series

#load needed packages
library(BoolNet)
library(reshape2)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggrepel)
library(dplyr)

n.obs <- 100000 # number of observation we want to test up to  
p <- 25 # number of variables

#this is a threshold to flip the node 
chance.to.flip <- 0.01

# 22 time-series of noise for 22 variables
zeta <- matrix(runif(n.obs * p, 0, 1), ncol = p)
zeta <- t(zeta) # adjust zeta to a wide-format
time.series <- matrix(rep(0, p * n.obs), nrow =p , ncol =  n.obs)
time.series[,1] <- runif(p,0, 1) # initial value of the three variables were randomly chosen
time.series[,1] <-ifelse(time.series[,1] >= .5, TRUE,FALSE) # if the randomly generated number is greater than .5, then we convert the value to TRUE, otherwise to FALSE

for (col in 2:n.obs){
    # simulate each time step based on the previous time step and Boolean functions
    time.series[1,col] = !(time.series[2,col-1] | time.series[8,col-1]) & (time.series[18,col-1] | time.series[2,col-1] | time.series[11,col-1]) #AKT1, !(BRAF|GSK3B) & (BCL2|BRAF|MAPK1)
    time.series[2,col] = time.series[10,col-1] #BRAF, MAP2K2
    time.series[3,col] = !(time.series[13,col-1]) & time.series[17,col-1] #CDKN2A, !(NPM1) & (TP53)
    time.series[4,col] = !(time.series[12,col-1] | time.series[19,col-1]) #CEBPA, !(MYC|SOX4)
    time.series[5,col] = !(time.series[20,col-1] | time.series[3,col-1] | time.series[9,col-1] | time.series[21,col-1]) #DNMT3A, !(CCND1|CDKN2A|HOXA9|MEIS1)
    time.series[6,col] = !(time.series[12,col-1]) #FBXW7, !MYC
    time.series[7,col] = !(time.series[4,col-1]) & (time.series[1,col-1] | time.series[12,col-1] | time.series[15,col-1] | time.series[16,col-1]) #FLT3, !(CEBPA) & (AKT1|MYC|PTPN11|STAT5A)
    time.series[8,col] = !(time.series[12,col-1] | time.series[20,col-1]) & (time.series[4,col-1] | time.series[17,col-1])	#GSK3B, !(MYC|CCND1) & (CEBPA|TP53)
    time.series[9,col] = time.series[21,col-1] #HOXA9, MEIS1
    time.series[10,col] = time.series[11,col] #MAP2K2, MAPK1
    time.series[11,col] = !(time.series[2,col-1] | time.series[4,col-1] | time.series[22,col-1] |time.series[8,col-1]) & (time.series[18,col-1] | time.series[12,col-1] | time.series[16,col-1]) #MAPK1, !(BRAF|CEBPA|ETV6|GSK3B) & (BCL2|MYC|STAT5A)
    time.series[12,col] = !(time.series[3,col-1]) & (time.series[20,col-1] | time.series[5,col-1]) #MYC, !(CDKN2A) & (CCND1|DNMT3A)
    time.series[13,col] = !(time.series[9,col-1]) & (time.series[3,col-1] | time.series[6,col-1]) #NPM1, !(HOXA9) & (CDKN2A|FBXW7)
    time.series[14,col] = time.series[2,col-1] #NRAS, BRAF
    time.series[15,col] = time.series[14,col-1] #PTPN11, NRAS
    time.series[16,col] = !(time.series[5,col-1]) & time.series[5,col-1] #STAT5A, !(DNMT3A) & (DNMT3A)	
    time.series[17,col] = !(time.series[18,col-1]) #TP53, !(BCL2) 
    time.series[18,col] = time.series[18,col-1] #BCL2, BCL2
    time.series[19,col] = time.series[19,col-1] #SOX4, SOX4
    time.series[20,col] = time.series[20,col-1] #CCND1, CCND1
    time.series[21,col] = time.series[21,col-1]  #MEIS1, MEIS1
    time.series[22,col] = time.series[22,col-1] #ETV6, ETV6
    time.series[23,col] =  !(time.series[3,col-1]) & (time.series[16,col-1] | time.series[12,col-1]| time.series[21,col-1] | time.series[20,col-1]) #Proliferation
    time.series[24,col] =  !(time.series[19,col-1]) & (time.series[4,col-1] |  time.series[22,col-1]) #Differentiation
    time.series[25,col] =  !(time.series[18,col-1]) & time.series[17,col-1] #Apoptosis
    
    
    #adding noise 
    time.series[1,col] <- ifelse(zeta[1,col] < chance.to.flip, 
                                 !time.series[1,col],
                                 time.series[1,col] )
    time.series[2,col] <- ifelse(zeta[2,col] < chance.to.flip, 
                                 !time.series[2,col],
                                 time.series[2,col] )
    time.series[3,col] <- ifelse(zeta[3,col] < chance.to.flip, 
                                 !time.series[3,col],
                                 time.series[3,col])
    time.series[4,col] <- ifelse(zeta[4,col] < chance.to.flip, 
                                 !time.series[4,col],
                                 time.series[,col] )
    time.series[5,col] <- ifelse(zeta[5,col] < chance.to.flip, 
                                 !time.series[5,col],
                                 time.series[5,col] )
    time.series[6,col] <- ifelse(zeta[6,col] < chance.to.flip, 
                                 !time.series[6,col],
                                 time.series[6,col])
    time.series[7,col] <- ifelse(zeta[7,col] < chance.to.flip, 
                                 !time.series[7,col],
                                 time.series[7,col] )
    time.series[8,col] <- ifelse(zeta[8,col] < chance.to.flip, 
                                 !time.series[8,col],
                                 time.series[8,col] )
    time.series[9,col] <- ifelse(zeta[9,col] < chance.to.flip, 
                                 !time.series[9,col],
                                 time.series[9,col])
    time.series[10,col] <- ifelse(zeta[10,col] < chance.to.flip, 
                                  !time.series[10,col],
                                  time.series[10,col])
    time.series[11,col] <- ifelse(zeta[11,col] < chance.to.flip, 
                                  !time.series[11,col],
                                  time.series[11,col])
    time.series[12,col] <- ifelse(zeta[12,col] < chance.to.flip, 
                                  !time.series[12,col],
                                  time.series[12,col])
    time.series[13,col] <- ifelse(zeta[13,col] < chance.to.flip, 
                                  !time.series[13,col],
                                  time.series[13,col])
    time.series[14,col] <- ifelse(zeta[14,col] < chance.to.flip, 
                                  !time.series[14,col],
                                  time.series[14,col])
    time.series[15,col] <- ifelse(zeta[15,col] < chance.to.flip, 
                                  !time.series[15,col],
                                  time.series[15,col])
    time.series[16,col] <- ifelse(zeta[16,col] < chance.to.flip, 
                                  !time.series[16,col],
                                  time.series[16,col])
    time.series[17,col] <- ifelse(zeta[17,col] < chance.to.flip, 
                                  !time.series[17,col],
                                  time.series[17,col])
    time.series[18,col] <- ifelse(zeta[18,col] < chance.to.flip, 
                                  !time.series[18,col],
                                  time.series[18,col])
    time.series[19,col] <- ifelse(zeta[19,col] < chance.to.flip, 
                                  !time.series[19,col],
                                  time.series[19,col])
    time.series[20,col] <- ifelse(zeta[20,col] < chance.to.flip, 
                                  !time.series[20,col],
                                  time.series[20,col] )
    time.series[21,col] <- ifelse(zeta[21,col] < chance.to.flip, 
                                  !time.series[21,col],
                                  time.series[21,col] )
    time.series[22,col] <- ifelse(zeta[22,col] < chance.to.flip, 
                                  !time.series[22,col],
                                  time.series[22,col])
    time.series[23,col] <- ifelse(zeta[23,col] < chance.to.flip, 
                                  !time.series[23,col],
                                  time.series[23,col])
    time.series[24,col] <- ifelse(zeta[24,col] < chance.to.flip, 
                                  !time.series[24,col],
                                  time.series[24,col])
    time.series[25,col] <- ifelse(zeta[25,col] < chance.to.flip, 
                                  !time.series[25,col],
                                  time.series[25,col])
  }
  
  # adjust the time-series to long-format
  time.series <- t(time.series)
  time.series <- data.frame(time.series)
  names(time.series) <- c("AKT1","BRAF", "CDKN2A", "CEBPA", "DNMT3A", "FBXW7", "FLT3", "GSK3B", "HOXA9", "MAP2K2", "MAPK1",
                          "MYC", "NPM1", "NRAS", "PTPN11", "STAT5A", "TP53", "BCL2", "SOX4", "CCND1", 
                          "MEIS1", "ETV6", "Proliferation", "Differentiation", "Apoptosis")
  
  #want the running value of n.obs rows, with the past 1000 values being considered
  time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))
  
  window = 1000
  time.series <- time.series %>%
    mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))
  
  upper_bound <- 50000  
  lower_bound <- upper_bound-1000   
  
############# testing the convergence  #############

###final scores ###
#getting the relative error points to check if the score converges
s1 <- mean(time.series$running_mean[500:1000])  #for 1000
s2 <- mean(time.series$running_mean[500:2000])  #for 2000
s3 <- mean(time.series$running_mean[500:4000])  #for 4000
s4 <- mean(time.series$running_mean[500:8000])  #for 8000
s5 <- mean(time.series$running_mean[500:16000]) #for 16000
s6 <- mean(time.series$running_mean[500:32000]) #for 32000
s7 <- mean(na.omit(time.series$running_mean[500:64000]))#for 64000
s8 <- mean(na.omit(time.series$running_mean[500:100000]))#for 100000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
point_6 <- ((s7-s6)/s6)
point_7 <- ((s8-s7)/s7)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6,point_7)
observation_intervals <- c(2000, 4000, 8000, 16000, 32000, 64000,100000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(observation_intervals, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = observation_intervals, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot for Final Scores (DNMT3A Model 1)", x = "observation intervals", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,110000), ylim = c(-2, 1))

###proliferation ###
#getting the relative error points to check if the score converges
s1 <- mean(time.series$Proliferation[500:1000])  #for 1000
s2 <- mean(time.series$Proliferation[500:2000])  #for 2000
s3 <- mean(time.series$Proliferation[500:4000])  #for 4000
s4 <- mean(time.series$Proliferation[500:8000])  #for 8000
s5 <- mean(time.series$Proliferation[500:16000]) #for 16000
s6 <- mean(time.series$Proliferation[500:32000]) #for 32000
s7 <- mean(na.omit(time.series$Proliferation[500:64000]))#for 64000
s8 <- mean(na.omit(time.series$Proliferation[500:100000]))#for 100000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
point_6 <- ((s7-s6)/s6)
point_7 <- ((s8-s7)/s7)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6,point_7)
observation_intervals <- c(2000, 4000, 8000, 16000, 32000, 64000,100000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(observation_intervals, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = observation_intervals, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot for Proliferation (DNMT3A Model 1)", x = "observation intervals", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,110000), ylim = c(-2, 1))

### differentiation ###
#getting the relative error points to check if the score converges
s1 <- mean(time.series$Differentiation[500:1000])  #for 1000
s2 <- mean(time.series$Differentiation[500:2000])  #for 2000
s3 <- mean(time.series$Differentiation[500:4000])  #for 4000
s4 <- mean(time.series$Differentiation[500:8000])  #for 8000
s5 <- mean(time.series$Differentiation[500:16000]) #for 16000
s6 <- mean(time.series$Differentiation[500:32000]) #for 32000
s7 <- mean(na.omit(time.series$Differentiation[500:64000]))#for 64000
s8 <- mean(na.omit(time.series$Differentiation[500:100000]))#for 100000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
point_6 <- ((s7-s6)/s6)
point_7 <- ((s8-s7)/s7)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6,point_7)
observation_intervals <- c(2000, 4000, 8000, 16000, 32000, 64000,100000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(observation_intervals, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = observation_intervals, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot for Differentiation (DNMT3A Model 1)", x = "observation intervals", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,110000), ylim = c(-2, 1))

### differentiation ###
#getting the relative error points to check if the score converges
s1 <- mean(time.series$Differentiation[500:1000])  #for 1000
s2 <- mean(time.series$Differentiation[500:2000])  #for 2000
s3 <- mean(time.series$Differentiation[500:4000])  #for 4000
s4 <- mean(time.series$Differentiation[500:8000])  #for 8000
s5 <- mean(time.series$Differentiation[500:16000]) #for 16000
s6 <- mean(time.series$Differentiation[500:32000]) #for 32000
s7 <- mean(na.omit(time.series$Differentiation[500:64000]))#for 64000
s8 <- mean(na.omit(time.series$Differentiation[500:100000]))#for 100000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
point_6 <- ((s7-s6)/s6)
point_7 <- ((s8-s7)/s7)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6,point_7)
observation_intervals <- c(2000, 4000, 8000, 16000, 32000, 64000,100000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(observation_intervals, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = observation_intervals, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot for Differentiation (DNMT3A Model 1)", x = "observation intervals", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,110000), ylim = c(-2, 1))

### apoptosis ###
#getting the relative error points to check if the score converges
s1 <- mean(time.series$Apoptosis[500:1000])  #for 1000
s2 <- mean(time.series$Apoptosis[500:2000])  #for 2000
s3 <- mean(time.series$Apoptosis[500:4000])  #for 4000
s4 <- mean(time.series$Apoptosis[500:8000])  #for 8000
s5 <- mean(time.series$Apoptosis[500:16000]) #for 16000
s6 <- mean(time.series$Apoptosis[500:32000]) #for 32000
s7 <- mean(na.omit(time.series$Apoptosis[500:64000]))#for 64000
s8 <- mean(na.omit(time.series$Apoptosis[500:100000]))#for 100000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
point_6 <- ((s7-s6)/s6)
point_7 <- ((s8-s7)/s7)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5,point_6,point_7)
observation_intervals <- c(2000, 4000, 8000, 16000, 32000, 64000,100000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(observation_intervals, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = observation_intervals, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot for Apoptosis (DNMT3A Model 1)", x = "observation intervals", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,110000), ylim = c(-2, 1))
