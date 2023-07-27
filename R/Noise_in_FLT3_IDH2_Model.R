# load needed packages
library(BoolNet)
library(reshape2)

n.obs <- 20000 # number of observation 
p <- 24 # number of variables

#randomly generated number from 0 to 1 (uniformly distributed) for each time-step, 
#so if we have a probability for flipping the node, we can compare this random number with the probablity, e.g., .2, 
#and when the random number is less than .2, we will flip the node. 
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

time.series[,1] <- runif(p,0, 1) # initial value of the 24 variables are randomly chosen

time.series[,1] <-ifelse(time.series[,1] >= .5, TRUE,FALSE) # if the randomly generated number is greater than .5, then we convert the value to TRUE, otherwise to FALSE

# this is a threshold to flip the node 
chance.to.flip <- 0.1

for (col in 2:n.obs){
  # simulate each time step based on the previous time step and Boolean functions
  time.series[1,col] <-time.series[1,col-1] 								       	        #FLT3= FLT3
  time.series[2,col] <-!time.series[1,col-1] 									              #CEBPA = !FLT3
  time.series[3,col] <-time.series[1,col-1] 									              #RUNX1 = FLT3
  time.series[4,col] <-!time.series[3,col-1] 									              #RAD21 = !RUNX1
  time.series[5,col] <-time.series[5,col-1] 								            	  #AKT = AKT
  time.series[6,col] <-time.series[24,col-1] 									              #NRAS = PTPN11
  time.series[7,col] <-time.series[6,col-1] 									              #BRAF = NRAS
  time.series[8,col] <-time.series[8,col-1] 									              #AMPK = AMPK
  time.series[9,col] <-time.series[9,col-1] 									              #NPM1 = NPM1
  time.series[10,col] <-time.series[7,col-1] 							            		  #MEK = BRAF
  time.series[11,col] <-!time.series[5,col-1] & time.series[8,col-1] 			  #FOXO = !AKT & AMPK
  time.series[12,col] <-time.series[9,col-1] 									              #FBXW7 = NPM1
  time.series[13,col] <-time.series[10,col-1] 								              #ERK = MEK
  time.series[14,col] <-!(time.series[3,col-1] | time.series[2,col-1] | time.series[12,col-1]) & time.series[13,col-1] 	#MYC = !(RUNX1 | CEBPA | FBXW7) & ERK
  time.series[15,col] <-time.series[15,col-1] 								              #IDH2 = IDH2
  time.series[16,col] <-time.series[11,col-1] 							              	#IDH1 = FOXO
  time.series[17,col] <-time.series[16,col-1] & time.series[15,col-1]		    #OXO2 = IDH1 & IDH2
  time.series[18,col] <-!time.series[14,col-1]  & time.series[9,col-1]			#CDKN2A = !MYC & NPM1
  time.series[19,col] <-time.series[17,col-1] & time.series[8,col-1]		   	#TET2 = OXO2 & AMPK
  time.series[20,col] <- !time.series[18,col-1] & time.series[22,col-1] 		#MDM2 = !CDKN2A & TP53
  time.series[21,col] <-time.series[19,col-1] 									            #WT1 = TET2
  time.series[22,col] <-!time.series[20,col-1] 								            	#TP53 = !MDM2
  time.series[23,col] <-!time.series[22,col-1] & time.series[13,col-1]			#BCL2 = !TP53 & ERK
  time.series[24,col] <-time.series[1,col-1] 							              		#PTPN11 = FLT3
  
  #these did not work
    #time.series[25,col] <-!(time.series[11,col-1] | time.series[21,col-1] | time.series[18,col-1] | time.series[22,col-1]) & time.series[14,col-1]    #Proliferation = !(FOXO | WT1 | CDKN2A | TP53) & MYC
    #time.series[26,col] <-time.series[3,col-1] & time.series[2,col-1] #Differentiation = RUNX1 & CEBPA
    #time.series[27,col] <-!time.series[23,col-1] & time.series[22,col-1] #Apoptosis = !BCL2 & TP53
    
  #add noise, if the random number in zeta is less than .2, then flip the node by applying  the NOT

  
  time.series[1,col] <- ifelse(zeta[1, col] < chance.to.flip, 
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
}

#adjust the time-series to long-format
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

#score_vector <- time.series$Score[1000:2000]  #saving the values from 1000 to 2000
score_vector <- time.series$Score              #saving the all the values

# plot the distribution of the simulated results
par(mfrow = c(3,2))
hist(time.series[15000:n.obs,28], 
     xlim = c(-8,2), 
     #ylim = c(0, 5000),
     main = "Histgram of network score [150000:200000]",
     xlab = "Network score"
     ) # plot the distribution of score_vector

hist(time.series[0:n.obs,28], 
     xlim = c(-8,2), 
     ylim = c(0, 5000),
     main = "Histgram of network score [0:20000]",
     xlab = "Network score"
) # plot the distribution of score_vector

hist(time.series[1:n.obs,27], 
    # xlim = c(-8,2), 
     #ylim = c(0, 50000),
     main = "Histgram of network score [15000:20000]",
     xlab = "Apoptosis"
) # plot the distribution of score_vector

hist(time.series[1:n.obs,26], 
     # xlim = c(-8,2), 
     #ylim = c(0, 50000),
     main = "Histgram of network score [15000:20000]",
     xlab = "Differentiation"
) # plot the distribution of score_vector

hist(time.series[1:n.obs,25], 
     # xlim = c(-8,2), 
     #ylim = c(0, 50000),
     main = "Histgram of network score [15000:20000]",
     xlab = "Proliferation"
) # plot the distribution of score_vector


plot(seq(1,n.obs), time.series[1:n.obs,28], xlab = "time points", ylab = "Network score")

