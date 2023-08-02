#load needed packages
library(BoolNet)
library(reshape2)
library(dplyr)
library(zoo)
library(Dict)
library(ggplot2)

#create the dictonary
genes <- Dict$new(  
  FLT3 = "time.series[1,col] <- 1" ,
  CEBPA = "time.series[2,col] <- 0" ,
  RUNX1 = "time.series[3,col] <- 1" ,			
  RAD21 = "time.series[4,col] <- 0" ,
  #AKT 
  NRAS = "time.series[6,col] <- 1",
  BRAF = "time.series[7,col] <- 1",
  #AMPK,
  NPM1 = "time.series[9,col] <- 0",
  MEK = "time.series[10,col] <- 0",
  FOXO = "time.series[11,col] <- 1",
  FBXW7 = "time.series[12,col] <- 1",  
  #ERK,
  MYC = "time.series[14,col] <- 1",
  IDH2 = "time.series[15,col] <- 1", 
  IDH1 = "time.series[16,col] <- 1",
  #OXO2,
  CDKN2A = "time.series[18,col] <- 0", 
  TET2 = "time.series[19,col] <- 0",
  MDM2 = "time.series[20,col] <- 1", 
  WT1 = "time.series[21,col] <- 1",			
  TP53 = "time.series[22,col] <- 0",				
  BCL2 = "time.series[23,col] <- 1",
  PTPN11 = "time.series[24,col] <- 1",	
  .class = "character",
  .overwrite = TRUE)

#genes["FLT3"]= time.series[j,col] <- 1

#create the function geneswitch()
geneswitch <- function(mutation_profile) {
  
  #seperate the values by the commas
  mutation_profile <- unlist(strsplit(mutation_profile,","))
  
  #remove repeated values
  mutation_profile <- unique(mutation_profile)
  
  #save as a list
  mutation_list <- as.list(mutation_profile)
  
  #so now it should look like this (example)
  #List values: 
  #"FLT3"
  #"DNMT3A"
  #"NPM1"
  
  #now we loop through the list 
  for (p in mutation_list) {
    cat(genes[p], sep="\n")
  }}


#combine the table 7 rows into columns, so one row is equal to one patient 
patient_profile_df <- s7_df %>%
  dplyr::group_by(labId) %>%
  dplyr::summarise(symbol = paste(symbol, collapse = ","))

#FOR NOW IM MAKING A TEST DATAFRAME
patient_profile_test_df <- patient_profile_df[1:5, ]


############# making the time series #############

n.obs <- 40500 # number of observation aka steps 
n.half <- (n.obs/2) #half the steps (for plot making)
p <- 24 # number of variables

# randomly generated number from 0 to 1 (uniformly distributed) for each time-step, so if we have a probability for flipping the node, 
#we can compare this random number with the probablity, e.g., .2, and when the random number is less than .2, we will flip the node. 
#need to definately fix this part
zeta <- cbind(runif(n.obs,0, 1), #1
              runif(n.obs,0, 1), #2
              runif(n.obs,0, 1), #3
              runif(n.obs,0, 1), #4
              runif(n.obs,0, 1), #5
              runif(n.obs,0, 1), #6
              runif(n.obs,0, 1), #7
              runif(n.obs,0, 1), #8
              runif(n.obs,0, 1), #9 
              runif(n.obs,0, 1), #10
              runif(n.obs,0, 1), #11
              runif(n.obs,0, 1), #12
              runif(n.obs,0, 1), #13
              runif(n.obs,0, 1), #14
              runif(n.obs,0, 1), #15
              runif(n.obs,0, 1), #16
              runif(n.obs,0, 1), #17
              runif(n.obs,0, 1), #18
              runif(n.obs,0, 1), #19
              runif(n.obs,0, 1), #20
              runif(n.obs,0, 1), #21
              runif(n.obs,0, 1), #22
              runif(n.obs,0, 1), #23
              runif(n.obs,0, 1)) #24
# 24 time-series of noise for 24 variables

zeta <- t(zeta) # adjust zeta to a wide-format

time.series <- matrix(rep(0, p * n.obs), nrow =p , ncol =  n.obs)

time.series[,1] <- runif(p,0, 1) # initial value of the three variables were randomly chosen

time.series[,1] <-ifelse(time.series[,1] >= .5, TRUE,FALSE) # if the randomly generated number is greater than .5, then we convert the value to TRUE, otherwise to FALSE

# this is a threshold to flip the node 
chance.to.flip <- 0.01

for (col in 2:n.obs){
    # simulate each time step based on the previous time step and Boolean functions
    time.series[1,col] <-time.series[1,col-1] 									            #FLT3= FLT3
    time.series[2,col] <-!time.series[1,col-1] 									            #CEBPA = !FLT3
    time.series[3,col] <-time.series[1,col-1] 									            #RUNX1 = FLT3
    time.series[4,col] <-!time.series[3,col-1] 									            #RAD21 = !RUNX1
    time.series[5,col] <-time.series[5,col-1] 									            #AKT = AKT
    time.series[6,col] <-time.series[24,col-1] 									            #NRAS = PTPN11
    time.series[7,col] <-time.series[6,col-1] 									            #BRAF = NRAS
    time.series[8,col] <-time.series[8,col-1] 									            #AMPK = AMPK
    time.series[9,col] <-time.series[9,col-1] 								            	#NPM1 = NPM1
    time.series[10,col] <-time.series[7,col-1] 							            		#MEK = BRAF
    time.series[11,col] <-!time.series[5,col-1] & time.series[8,col-1] 			#FOXO = !AKT & AMPK
    time.series[12,col] <-time.series[9,col-1] 									            #FBXW7 = NPM1
    time.series[13,col] <-time.series[10,col-1] 									          #ERK = MEK
    time.series[14,col] <-!(time.series[3,col-1] | time.series[2,col-1] | time.series[12,col-1]) & time.series[13,col-1] 	#MYC = !(RUNX1 | CEBPA | FBXW7) & ERK
    time.series[15,col] <-time.series[15,col-1] 									          #IDH2 = IDH2
    time.series[16,col] <-time.series[11,col-1] 									          #IDH1 = FOXO
    time.series[17,col] <-time.series[16,col-1] & time.series[15,col-1]			#OXO2 = IDH1 & IDH2
    time.series[18,col] <-!time.series[14,col-1]  & time.series[9,col-1]		#CDKN2A = !MYC & NPM1
    time.series[19,col] <-time.series[17,col-1] & time.series[8,col-1]			#TET2 = OXO2 & AMPK
    time.series[20,col] <- !time.series[18,col-1] & time.series[22,col-1] 	#MDM2 = !CDKN2A & TP53
    time.series[21,col] <-time.series[19,col-1] 									          #WT1 = TET2
    time.series[22,col] <-!time.series[20,col-1] 									          #TP53 = !MDM2
    time.series[23,col] <-!time.series[22,col-1] & time.series[13,col-1]		#BCL2 = !TP53 & ERK
    time.series[24,col] <-time.series[1,col-1] 								            	#PTPN11 = FLT3
  
  # add noise, if the random number in zeta is less than .2, then flip the node by applying  the NOT
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
                                time.series[20,col])
  time.series[21,col] <- ifelse(zeta[21,col] < chance.to.flip, 
                                !time.series[21,col],
                                time.series[21,col])
  time.series[22,col] <- ifelse(zeta[22,col] < chance.to.flip, 
                                !time.series[22,col],
                                time.series[22,col])
  time.series[23,col] <- ifelse(zeta[23,col] < chance.to.flip, 
                                !time.series[23,col],
                                time.series[23,col])
  time.series[24,col] <- ifelse(zeta[24,col] < chance.to.flip, 
                                !time.series[24,col],
                                time.series[24,col])
  
#now we add the personalized patient mutations
#for (i in 1:nrow(patient_profile_test_df)){
  #stores the mutation profile as a vector
    #mutation_profile <- as.vector(patient_profile_test_df$symbol[i]) 
  
  #get the values from the dictionary and print them
  #cat(geneswitch(mutation_profile))}
  
}

############# dataframe and rolling mean #############

# adjust the time-series to long-format
time.series <- t(time.series)
time.series <- data.frame(time.series)
names(time.series) <- c("FLT3", "CEBPA", "RUNX1", "RAD21", "AKT", "NRAS", "BRAF", "AMPK", "NPM1", 
                        "MEK", "FOXO", "FBXW7", "ERK", "MYC", "IDH2", "IDH1", "OXO2", "CDKN2A", "TET2",
                        "MDM2", "WT1", "TP53", "BCL2", "PTPN11")

#getting the network score and saving it as a vector
time.series$Proliferation <- with(time.series, ((-1 * FOXO) + (-1 * WT1) + (-1 * CDKN2A) + (-1 * TP53) + MYC))
time.series$Differentiation <- with(time.series, (CEBPA + RUNX1))
time.series$Apoptosis <- with(time.series, ((-1 * BCL2) + TP53))
time.series$Score <- with(time.series, (Proliferation - (Apoptosis + Differentiation)))

#want the running value of n.obs rows, with the past 1000 values being considered
window = 1000
time.series <- time.series %>%
  mutate(running_mean = rollmean(time.series$Score, k = window, fill = NA))

#rolling median of the means 
rolling_median <- rollmedian(time.series$running_mean, k = window)

############# scatterplots and histograms #############

#diagrams
hist(time.series[window:n.obs,28],
     xlim = c(-8,2),
     ylim = c(0, 8000),
     main = "Histogram of running mean",
     xlab = "running mean"
) # plot the distribution of the running mean 

hist(rolling_median,
     xlim = c(-5,-2),
     ylim = c(0, 3000),
     main = "Histogram of rolling median",
     xlab = "rolling median"
) # plot the distribution of the rolling median 

#plot(seq(1,n.obs), time.series$running_mean[upper_bound:lower_bound], xlab = "time points", ylab = "running mean")

#getting the relative error points to check if the score converges
s1 <- mean(time.series$Score[500:1000])  #for 1000
s2 <- mean(time.series$Score[500:2000])  #for 2000
s3 <- mean(time.series$Score[500:4000])  #for 4000
s4 <- mean(time.series$Score[500:8000])  #for 8000
s5 <- mean(time.series$Score[500:16000]) #for 16000
s6 <- mean(time.series$Score[500:32000]) #for 32000
#s7 <- mean(time.series$Score[500:64000]) #for 64000

#getting the relative error points
point_1 <- ((s2-s1)/s1)
point_2 <- ((s3-s2)/s2)
point_3 <- ((s4-s3)/s3)
point_4 <- ((s5-s4)/s4)
point_5 <- ((s6-s5)/s5)
#point_6 <- ((s7-s6)/s6)

#saving into vectors
relative_error <- c(point_1, point_2, point_3, point_4,point_5)
time <- c(1000, 2000, 4000, 8000, 16000)

#create a new dataframe for the relative_error
mean_time_points <- data.frame(time, relative_error)

mean_time_scatterplot <- ggplot(data = mean_time_points, aes(x = time, y = relative_error)) +
  geom_point() +
  labs(title = "Relative Error Scatterplot", x = "time interval", y = "relative error") +
  theme_minimal()
mean_time_scatterplot + coord_cartesian(xlim = c(800,40000), ylim = c(-2, 1))

############# final score #############

#final score based on window value of running mean 
upper_bound <- 20000 
lower_bound <- 40000    
#the last 500 values (without NA) are considered in the mean 

final_score <- mean(time.series$running_mean[upper_bound:lower_bound], na.rm=TRUE)
print(final_score)
