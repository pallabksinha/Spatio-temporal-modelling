#DAta is from S1C2 
#Bath\paper purpose\simulation\Simulation with vaccination effect\simulation and predicted\S1C2\Without vaccination in linear predictor\X3 with V X1 and X2 without V


rm(list = ls())
set.seed(25)
#data generate 
#same initial condition same but 1 and 3 very close in W matrix

library(plotly)
#library(dplyr)
no_of_simulation<-10000

infected_first_simulation <- replicate(no_of_simulation,rep(0,156),simplify = F)
infected_second_simulation <- replicate(no_of_simulation,rep(0,156),simplify = F)
infected_third_simulation <- replicate(no_of_simulation,rep(0,156),simplify = F)

s_first_category <- c()
s_first_category[1]<- 10^6
s_second_category <- c()
s_second_category[1]<- 10^6
s_third_category<-c()
s_third_category[1]<-10^6

#beta1<- runif(no_of_simulation,.4,.5) ; beta2<- runif(no_of_simulation,.55,.56)
beta1<- runif(no_of_simulation,.45,.5) ; beta2<- beta1
nstep<-26*6

B <- 10^3
# v1_3<- seq(.89,.9,length.out=nstep) #3rd cat. e vaccination besi .7 theke .9
# v2_3<- seq(.9,.95,length.out=nstep)
v1_3<- seq(.7,.9,length.out=nstep) #3rd cat. e vaccination besi .7 theke .9
v2_3<- seq(.4,.7,length.out=nstep)
length(v1_3);length(v2_3)
X3 <- B*(1-0.85*v1_3*(1-v2_3)-.99*v1_3*v2_3)
# X3<-c()
# for(t in 1:nstep){X3[t] <- B}          #1st and 2nd cat. e vaccination nei
# length(X3)


# v1_2<- seq(.4,.4,length.out=nstep) 
# v2_2<- seq(.4,.4,length.out=nstep)
# length(v1_2);length(v2_2)
X2<-c()
for(t in 1:nstep){X2[t] <- B}          #1st and 2nd cat. e vaccination nei
length(X2)


# v1_1<- seq(.4,.4,length.out=nstep)
# v2_1<- seq(.4,.4,length.out=nstep)
# length(v1_1);length(v2_1)
X1<-c()
for(t in 1:nstep){X1[t] <- B}               
length(X1)



for(j in 1:no_of_simulation){
  beta_3_11=beta1[j]/(s_first_category[1])
  beta_3_22=beta1[j]/(s_second_category[1])
  beta_3_33=beta2[j]/(s_third_category[1])
  W=matrix(c(1,exp(-5),exp(-5),exp(-5),1,exp(-.005),exp(-5),exp(-.005),1),nrow=3,ncol=3,byrow=T)
  infected_first_category <- c()
  infected_first_category[1]<-6
  infected_second_category <-c()
  infected_second_category[1]<- 6
  infected_third_category<-c()
  infected_third_category[1]<-6
  
  N1<-s_first_category[1]+infected_first_category[1]
  N2<-s_second_category[1]+infected_second_category[1]
  N3<-s_third_category[1]+infected_third_category[1]
  
  # first we just want to see how the simulated data set looks like for one simulation
  nstep<-26*6 #time points
  #simulations-100
  
  gamma <- 1
  
  
  # Define functions for the differential equations
  AA1 <- function(s1, i1, i2, i3, beta1, beta2, beta3)  - s1 * (beta1 * i1 + beta2 * i2 + beta3 * i3)
  BB1 <- function(s1, i1, i2, i3, beta1, beta2, beta3) s1 * (beta1 * i1 + beta2 * i2 + beta3 * i3)  - gamma * i1
  AA2 <- function(s2, i1, i2, i3, beta1, beta2, beta3)  - s2 * (beta1 * i1 + beta2 * i2 + beta3 * i3)
  BB2 <- function(s2, i1, i2, i3, beta1, beta2, beta3) s2 * (beta1 * i1 + beta2 * i2 + beta3 * i3)  - gamma * i2
  AA3 <- function(s3, i1, i2, i3, beta1, beta2, beta3)  - s3 * (beta1 * i1 + beta2 * i2 + beta3 * i3)
  BB3 <- function(s3, i1, i2, i3, beta1, beta2, beta3) s3 * (beta1 * i1 + beta2 * i2 + beta3 * i3)  - gamma * i3
  
  
  TT=26
  for(t in 2:nstep)
  {
    beta11 <- beta_3_11*(1 + 0.5 * sin(2 * pi * t / TT))*W[1,1]
    beta22 <- beta_3_22*(1 + 0.5 * sin(2 * pi * t / TT))*W[2,2]
    beta33 <- beta_3_33*(1 + 0.5 * sin(2 * pi * t / TT))*W[3,3]
    beta12 <- beta11 * W[1,2]
    beta13 <- beta11 * W[1,3]
    beta21 <- beta22 * W[2,1]
    beta23 <- beta22 * W[2,3]
    beta31 <- beta33 * W[3,1]
    beta32 <- beta33 * W[3,2]
    
    
    s_first_category[t]<- s_first_category[t-1] +  AA1(s_first_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta11, beta12, beta13)+X1[t]
    infected_first_category[t]<- infected_first_category[t-1] +  BB1(s_first_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta11, beta12, beta13)
    infected_first_simulation[[j]][t]<-infected_first_category[t]
    
    s_second_category[t]<- s_second_category[t-1] +  AA2(s_second_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta21, beta22, beta23)+X2[t]
    infected_second_category[t]<- infected_second_category[t-1] +  BB2(s_second_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta21, beta22, beta23)
    infected_second_simulation[[j]][t]<-infected_second_category[t]
    
    s_third_category[t]<- s_third_category[t-1] +  AA3(s_third_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta31, beta32, beta33)+X3[t]
    infected_third_category[t]<- infected_third_category[t-1] +  BB3(s_third_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta31, beta32, beta33)
    infected_third_simulation[[j]][t]<-infected_third_category[t]
  }
}


# simulation from "simulation beta hist"

#######################################################################################################################################
######################################################################################################################################
#1st category

library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Bath/Simulation")
data<-read_excel("finaldata_S1_C2.xlsx")


training_upto<-length(infected_first_simulation[[1]])
length(infected_first_category)
nstep<-length(infected_first_category)

mae<- replicate(no_of_simulation,rep(0,156),simplify = F)
mean_ab_err_for_model<-c() ; mse<-c()
mean_err_for_model<-c()
coefficients1<-c() ; coefficients2<-c() ; coefficients3<-c() ; coefficients4<-c(); coefficients5<-c() ; coefficients6<-c()
N<-1000000
S1<-1000000

# infected_second_simulation[x]<-NULL 
# infected_first_simulation[x]<-NULL
# infected_third_simulation[x]<-NULL

for(j in 1:(no_of_simulation)){
  
  infected_first_category <- infected_first_simulation[[j]]
  infected_second_category<- infected_second_simulation[[j]]
  infected_third_category<-infected_third_simulation[[j]]
  
  cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
  cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
  cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
  
  ar_order<-1
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
  
  log_infected_first_category<-infected_first_category
  log_infected_first_category[log_infected_first_category == 0] <- 1
  log_infected_first_category<-log(log_infected_first_category)
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       "lag2"=dplyr::lag(log_infected_first_category, 2) )
  
  
  log_infected_second_category<-infected_second_category
  log_infected_second_category[log_infected_second_category == 0] <- 1
  log_infected_second_category<-log(log_infected_second_category)
  log_infected_second_category<-  exp(-5)*log_infected_second_category    #this is w[i,j]*infected_2nd_category
  
  log_infected_third_category<-infected_third_category
  log_infected_third_category[log_infected_third_category == 0] <- 1
  log_infected_third_category<-log(log_infected_third_category)
  log_infected_third_category<-exp(-5)*log_infected_third_category       #this is w[i,j]*infected_3rd_category
  
  
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       
                                                       "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                       "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
  
  infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")
  
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}
  
  
  
  data_sf<-cbind(ceiling(infected_first_category_success),ceiling(infected_first_category_failure))
  data_sf<-matrix(c(ceiling(infected_first_category_success),ceiling(infected_first_category_failure)),nrow = length(infected_first_category_success))
  
  
  t<- (ar_order+1):training_upto
  length(t)
  seasonality1<-cos(2*pi*t/52)
  seasonality2<-cos(2*pi*t/26)
  
  
  data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                    seasonality2,t)
  data<-cbind(data_sf,data)
  
  colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time")
  head(data,20)
  
  model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time, 
               family = binomial(link = "logit"), data = data)
  
  
  coefficients1[j]<-as.numeric(coefficients(model)[1])
  coefficients2[j]<-as.numeric(coefficients(model)[2])
  coefficients3[j]<-as.numeric(coefficients(model)[3])
  coefficients4[j]<-as.numeric(coefficients(model)[4])
  
  predicted<-model$fitted.values*(data$S+data$F)

coe2<-scale(coefficients2)
par(mfrow=c(1,3))


dens <- density(coe2,bw=2.2)
# Scale the y-axis to represent relative frequency
dens$y <- dens$y / sum(dens$y)
# Plot the density with y-axis as relative frequency
plot(dens, 
     ylab = "Relative Frequency",  type = "l",xlab=expression("Distribution of " * phi[1]^1),main="First Category")



#######################################################################################################################################################

#2nd category


library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Bath/Simulation")
data<-read_excel("finaldata_S1_C2.xlsx")


training_upto<-length(infected_first_simulation[[1]])
length(infected_first_category)
nstep<-length(infected_first_category)

mae<- replicate(no_of_simulation,rep(0,156),simplify = F)
mean_ab_err_for_model<-c() ; mse<-c()
mean_err_for_model<-c()
coefficients1<-c() ; coefficients2<-c() ; coefficients3<-c() ; coefficients4<-c(); coefficients5<-c() ; coefficients6<-c()
N<-1000000
S1<-1000000



for(j in 1:(no_of_simulation)){
  
  infected_first_category <- infected_second_simulation[[j]]
  infected_second_category<- infected_first_simulation[[j]]
  infected_third_category<-infected_third_simulation[[j]]
  
  cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
  cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
  cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
  
  ar_order<-1
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
  
  log_infected_first_category<-infected_first_category
  log_infected_first_category[log_infected_first_category == 0] <- 1
  log_infected_first_category<-log(log_infected_first_category)
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       "lag2"=dplyr::lag(log_infected_first_category, 2) )
  
  
  log_infected_second_category<-infected_second_category
  log_infected_second_category[log_infected_second_category == 0] <- 1
  log_infected_second_category<-log(log_infected_second_category)
  log_infected_second_category<-  exp(-5)*log_infected_second_category    #this is w[i,j]*infected_2nd_category
  
  log_infected_third_category<-infected_third_category
  log_infected_third_category[log_infected_third_category == 0] <- 1
  log_infected_third_category<-log(log_infected_third_category)
  log_infected_third_category<-exp(-.005)*log_infected_third_category       #this is w[i,j]*infected_3rd_category
  
  
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       
                                                       "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                       "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
  
  infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")
  
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}
  
  
  
  data_sf<-cbind(ceiling(infected_first_category_success),ceiling(infected_first_category_failure))
  data_sf<-matrix(c(ceiling(infected_first_category_success),ceiling(infected_first_category_failure)),nrow = length(infected_first_category_success))
  
  
  t<- (ar_order+1):training_upto
  length(t)
  seasonality1<-cos(2*pi*t/52)
  seasonality2<-cos(2*pi*t/26)
  
  
  data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                    seasonality2,t)
  data<-cbind(data_sf,data)
  
  colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time")
  head(data,20)
  
  model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time, 
               family = binomial(link = "logit"), data = data)
  
  
  coefficients1[j]<-as.numeric(coefficients(model)[1])
  coefficients2[j]<-as.numeric(coefficients(model)[2])
  coefficients3[j]<-as.numeric(coefficients(model)[3])
  coefficients4[j]<-as.numeric(coefficients(model)[4])

  predicted<-model$fitted.values*(data$S+data$F)

coe2<-scale(coefficients2)

# Generate density estimates
dens <- density(coe2,bw=1.1)
# Scale the y-axis to represent relative frequency
dens$y <- dens$y / sum(dens$y)
# Plot the density with y-axis as relative frequency
plot(dens, 
     ylab = "Relative Frequency", type = "l",main="Second Category",xlab=expression("Distribution of " * phi[1]^2))



############################################################################################################################################################
###############################################################################################################################################################
#3rd Category

#Data from "simulation 500" file

library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Bath/Simulation")
data<-read_excel("finaldata_S1_C2.xlsx")


training_upto<-length(infected_first_simulation[[1]])
length(infected_first_category)
nstep<-length(infected_first_category)

mae<- replicate(no_of_simulation,rep(0,156),simplify = F)
mean_ab_err_for_model<-c() ; mse<-c()
mean_err_for_model<-c()
coefficients1<-c() ; coefficients2<-c() ; coefficients3<-c() ; coefficients4<-c(); coefficients5<-c() ; coefficients6<-c()
N<-1000000
S1<-1000000


for(j in 1:(no_of_simulation)){
  
  infected_first_category <- infected_third_simulation[[j]]
  infected_second_category<- infected_second_simulation[[j]]
  infected_third_category<-infected_first_simulation[[j]]
  
  cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
  cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
  cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
  
  ar_order<-1
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
  
  log_infected_first_category<-infected_first_category
  log_infected_first_category[log_infected_first_category == 0] <- 1
  log_infected_first_category<-log(log_infected_first_category)
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       "lag2"=dplyr::lag(log_infected_first_category, 2) )
  
  
  log_infected_second_category<-infected_second_category
  log_infected_second_category[log_infected_second_category == 0] <- 1
  log_infected_second_category<-log(log_infected_second_category)
  log_infected_second_category<-  exp(-.005)*log_infected_second_category    #this is w[i,j]*infected_2nd_category
  
  log_infected_third_category<-infected_third_category
  log_infected_third_category[log_infected_third_category == 0] <- 1
  log_infected_third_category<-log(log_infected_third_category)
  log_infected_third_category<-exp(-5)*log_infected_third_category       #this is w[i,j]*infected_3rd_category
  
  
  infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                       
                                                       "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                       "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
  
  infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")
  
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}
  
  
  
  data_sf<-cbind(ceiling(infected_first_category_success),ceiling(infected_first_category_failure))
  data_sf<-matrix(c(ceiling(infected_first_category_success),ceiling(infected_first_category_failure)),nrow = length(infected_first_category_success))
  
  
  t<- (ar_order+1):training_upto
  length(t)
  seasonality1<-cos(2*pi*t/52)
  seasonality2<-cos(2*pi*t/26)
  
  
  data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                    seasonality2,t)
  data<-cbind(data_sf,data)
  
  colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time")
  head(data,20)
  
  model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time, 
               family = binomial(link = "logit"), data = data)
  
  
  coefficients1[j]<-as.numeric(coefficients(model)[1])
  coefficients2[j]<-as.numeric(coefficients(model)[2])
  coefficients3[j]<-as.numeric(coefficients(model)[3])
  coefficients4[j]<-as.numeric(coefficients(model)[4])
  # coefficients5[j]<-as.numeric(coefficients(model)[5])
  # coefficients6[j]<-as.numeric(coefficients(model)[6])
  
  
  
  predicted<-model$fitted.values*(data$S+data$F)
coe2<-scale(coefficients2)

# Generate density estimates
dens <- density(coe2,bw=1.1)
# Scale the y-axis to represent relative frequency
dens$y <- dens$y / sum(dens$y)
# Plot the density with y-axis as relative frequency
plot(dens, 
     ylab = "Relative Frequency", type = "l",main="Third Category",xlab=expression("Distribution of " * phi[1]^3))


############################################################################################################################################################
