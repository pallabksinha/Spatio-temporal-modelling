rm(list = ls())
set.seed(25)
#data generate 
#same initial condition same but 1 and 3 very close in W matrix

library(plotly)
#library(dplyr)
no_of_simulation<-1
nstep<-26*6
infected_first_simulation <- replicate(no_of_simulation,rep(0,nstep),simplify = F)
infected_second_simulation <- replicate(no_of_simulation,rep(0,nstep),simplify = F)
infected_third_simulation <- replicate(no_of_simulation,rep(0,nstep),simplify = F)

s_first_category <- c()
s_first_category[1]<- 10^6
s_second_category <- c()
s_second_category[1]<- 10^6
s_third_category<-c()
s_third_category[1]<-10^6

#beta1<- runif(no_of_simulation,.4,.5) ; beta2<- runif(no_of_simulation,.55,.56)
beta1<- .5 ; beta2<- .5

nstep<-26*6

B <- 10^3
v1_3<- seq(.7,.9,length.out=nstep) #3rd cat. e vaccination besi .7 theke .9
v2_3<- seq(.4,.7,length.out=nstep)
length(v1_3);length(v2_3)
X3 <- B*(1-0.85*v1_3*(1-v2_3)-.99*v1_3*v2_3)
length(X3)


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
    
    
    s_first_category[t]<- s_first_category[t-1] +  AA1(s_first_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta11, beta12, beta13) + X1[t]
    infected_first_category[t]<- infected_first_category[t-1] +  BB1(s_first_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta11, beta12, beta13)
    infected_first_simulation[[j]][t]<-infected_first_category[t]
    
    s_second_category[t]<- s_second_category[t-1] +  AA2(s_second_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta21, beta22, beta23) + X2[t]
    infected_second_category[t]<- infected_second_category[t-1] +  BB2(s_second_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta21, beta22, beta23)
    infected_second_simulation[[j]][t]<-infected_second_category[t]
    
    s_third_category[t]<- s_third_category[t-1] +  AA3(s_third_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta31, beta32, beta33) + X3[t]
    infected_third_category[t]<- infected_third_category[t-1] +  BB3(s_third_category[t-1], infected_first_category[t-1], infected_second_category[t-1], infected_third_category[t-1], beta31, beta32, beta33)
    infected_third_simulation[[j]][t]<-infected_third_category[t]
  }
  infected_first_simulation[[j]][1]<-infected_first_category[1]
  infected_second_simulation[[j]][1]<-infected_second_category[1]
  infected_third_simulation[[j]][1]<-infected_third_category[1]
}

time <- seq(1, nstep, by = 1)

# Simulate some disease counts for three regions
region1_count <- infected_first_simulation[[1]]
region2_count <- infected_second_simulation[[1]]
region3_count <- infected_third_simulation[[1]]

# Create region IDs
region <- factor(c(rep(1, nstep), rep(2, nstep), rep(3, nstep)))

# Combine data into a dataframe
data <- data.frame(
  time = rep(time, 3), 
  region = region, 
  count = c(region1_count, region2_count, region3_count)
)
plot_ly(data, 
        x = ~time, 
        y = ~region, 
        z = ~count, 
        type = 'scatter3d', 
        mode = 'lines', 
        color = ~region, 
        line = list(width = 4)) %>%
  layout(scene = list(
    xaxis = list(title = 'Time (in biweeks)'),
    yaxis = list(title = 'Regions', tickvals = c(1, 2, 3), ticktext = c('Region-1', 'Region-2', 'Region-3')),
    zaxis = list(title = 'Disease Count',range = c(0, 500))
  ))

########################################################################################################################
##########################################################################################################################



#################################################################################################################
#################################################################################################################
#Prediction of first region


set.seed(25)
library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
training_upto<-nstep
infected_first_category <- infected_first_simulation[[1]][1:training_upto]
plot(infected_first_category,type = "l",col="red",ylab="Case Count",xlab = "Time(biweekly)")
acf(infected_first_category,xlab="Lag",main="Acf plot for counts",ci=0.99)
pacf(infected_first_category,   xlab="Lag",main="PAcf plot for counts",ci=0.99)

length(infected_first_category)
nstep<-length(infected_first_category)

N<-10^6
S1<-10^6

cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+Instep
length(cum)
cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+Instep)
length(cum)
cum<-cum[-nstep]   #S_nstep is deleted because upto S_545 is required
length(cum)



ar_order<-1
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

# t<-1:length(infected_first_category)
#infected_first_category_minus_last_term<-data$Bath[-nstep]  #prob. (ie pie) ber korar somoi I1,..,,I597 lagbe.. pie1 is functn of I1,.., pie2 is functn of I2 etc.. so last er ta delete kore dewoa holo

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


#distance between [leicester,coventry]=52955.13,[Birmingham,coventry]=55266.08 (in meters)
#distance between [leicester,coventry]=52.95513,[Birmingham,coventry]=55.26608 (in kms)
#so exp(-52.95513) = 1.004337e-23 ; exp(-55.26608) = 9.959704e-25


infected_second_category<-infected_second_simulation[[1]][1:training_upto]
infected_third_category<-infected_third_simulation[[1]][1:training_upto]

log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-W[1,2]*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<- W[1,3]*log_infected_third_category       #this is w[i,j]*infected_3rd_category

infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     
                                                     "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))


# birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
# xn<-as.numeric(birth_data[15, 3:((training_upto/26)+2)])
# divided_vector <- xn / 26
# biweekly_birth_data<- rep(divided_vector, each = 26)
# length(biweekly_birth_data)
# biweekly_birth_data_in_inverselogit1<- biweekly_birth_data[-length(biweekly_birth_data)]
# #biweekly_birth_data_in_inverselogit<- biweekly_birth_data_in_inverselogit[-c(1:(ar_order-1))]
# if (ar_order == 1) {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1
# } else {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1[-c(1:(ar_order-1))]
# }
# length(biweekly_birth_data_in_inverselogit)
t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)

# baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
# baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
# if (ar_order == 1) {
#   baby_boom_effect <- baby_boom_effect1
# } else {
#   baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
# }
# length(baby_boom_effect)






head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)

#length(biweekly_birth_data_in_inverselogit)
length(seasonality1)
length(seasonality2)
length(t)
#length(baby_boom_effect)

length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)

nrow(data_sf)

# data<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),
#                infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
#                biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect
# ),
# nrow = length(infected_first_category_in_inverselogit$lag1))

X_t <- X1[-c(1:(ar_order-1))]
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  seasonality2,t,X_t)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time","birth_cohort")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)

# 
# 
# acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.50)
# pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.50)
# 
# 
predicted<-model$fitted.values*(data$S+data$F)
# data.frame( predicted,infected_first_category_success)


plot(infected_first_category_success, type = "l", col = "black", ylim = c(0, max(infected_first_category, predicted)),xlab = "Time",ylab = "Counts")
lines(predicted, type = "l", col = "red")  # Add the line for predicted_infection
points(predicted, col = "red", pch = 8, cex = .25)  # Add the star-shaped points
legend("topright",legend = c("Real","Estimated"),col=c("black","red"),lty =c(1,1),cex=.4 )

# Box_ljung_test<- Box.test(model$residuals,lag=500,type = "Ljung-Box")
# Box_ljung_test



#################################################################################################################
#################################################################################################################
#Prediction of second region

set.seed(25)
library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
training_upto<-nstep
infected_first_category <- infected_second_simulation[[1]][1:training_upto]
# plot(infected_first_category,type = "l",col="red",ylab="Case Count",xlab = "Time(biweekly)")
# acf(infected_first_category,xlab="Lag",main="Acf plot for counts",ci=0.99)
# pacf(infected_first_category,   xlab="Lag",main="PAcf plot for counts",ci=0.99)

length(infected_first_category)
nstep<-length(infected_first_category)

N<-10^6
S1<-10^6

cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+Instep
length(cum)
cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+Instep)
length(cum)
cum<-cum[-nstep]   #S_nstep is deleted because upto S_545 is required
length(cum)



ar_order<-1
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

# t<-1:length(infected_first_category)
#infected_first_category_minus_last_term<-data$Bath[-nstep]  #prob. (ie pie) ber korar somoi I1,..,,I597 lagbe.. pie1 is functn of I1,.., pie2 is functn of I2 etc.. so last er ta delete kore dewoa holo

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


#distance between [leicester,coventry]=52955.13,[Birmingham,coventry]=55266.08 (in meters)
#distance between [leicester,coventry]=52.95513,[Birmingham,coventry]=55.26608 (in kms)
#so exp(-52.95513) = 1.004337e-23 ; exp(-55.26608) = 9.959704e-25


infected_second_category<-infected_first_simulation[[1]][1:training_upto]
infected_third_category<-infected_third_simulation[[1]][1:training_upto]

log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-W[2,1]*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<- W[2,3]*log_infected_third_category       #this is w[i,j]*infected_3rd_category

infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     
                                                     "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))


# birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
# xn<-as.numeric(birth_data[15, 3:((training_upto/26)+2)])
# divided_vector <- xn / 26
# biweekly_birth_data<- rep(divided_vector, each = 26)
# length(biweekly_birth_data)
# biweekly_birth_data_in_inverselogit1<- biweekly_birth_data[-length(biweekly_birth_data)]
# #biweekly_birth_data_in_inverselogit<- biweekly_birth_data_in_inverselogit[-c(1:(ar_order-1))]
# if (ar_order == 1) {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1
# } else {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1[-c(1:(ar_order-1))]
# }
# length(biweekly_birth_data_in_inverselogit)
t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)

# baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
# baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
# if (ar_order == 1) {
#   baby_boom_effect <- baby_boom_effect1
# } else {
#   baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
# }
# length(baby_boom_effect)






head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)

#length(biweekly_birth_data_in_inverselogit)
length(seasonality1)
length(seasonality2)
length(t)
#length(baby_boom_effect)

length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)

nrow(data_sf)

# data<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),
#                infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
#                biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect
# ),
# nrow = length(infected_first_category_in_inverselogit$lag1))

X_t <- X2[-c(1:(ar_order-1))]
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  seasonality2,t,X_t)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time","birth_cohort")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)

# 
# 
# acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.50)
# pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.50)
# 
# 
predicted2<-model$fitted.values*(data$S+data$F)
# data.frame( predicted,infected_first_category_success)


plot(infected_first_category_success, type = "l", col = "black", ylim = c(0, max(infected_first_category, predicted2)),xlab = "Time",ylab = "Counts")
lines(predicted2, type = "l", col = "red")  # Add the line for predicted_infection
points(predicted2, col = "red", pch = 8, cex = .25)  # Add the star-shaped points
legend("topright",legend = c("Real","Estimated"),col=c("black","red"),lty =c(1,1),cex=.4 )

# Box_ljung_test<- Box.test(model$residuals,lag=500,type = "Ljung-Box")
# Box_ljung_test


#################################################################################################################
#################################################################################################################
#Prediction of third region


set.seed(25)
library(readxl)
library(glarma)
library(psych)
library(xts)
library(readxl)
library(forecast)
library(seastests)
library(Rmpfr)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
training_upto<-nstep
infected_first_category <- infected_third_simulation[[1]][1:training_upto]
# plot(infected_first_category,type = "l",col="red",ylab="Case Count",xlab = "Time(biweekly)")
# acf(infected_first_category,xlab="Lag",main="Acf plot for counts",ci=0.99)
# pacf(infected_first_category,   xlab="Lag",main="PAcf plot for counts",ci=0.99)

length(infected_first_category)
nstep<-length(infected_first_category)

N<-10^6
S1<-10^6

cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+Instep
length(cum)
cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+Instep)
length(cum)
cum<-cum[-nstep]   #S_nstep is deleted because upto S_545 is required
length(cum)



ar_order<-1
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

# t<-1:length(infected_first_category)
#infected_first_category_minus_last_term<-data$Bath[-nstep]  #prob. (ie pie) ber korar somoi I1,..,,I597 lagbe.. pie1 is functn of I1,.., pie2 is functn of I2 etc.. so last er ta delete kore dewoa holo

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


#distance between [leicester,coventry]=52955.13,[Birmingham,coventry]=55266.08 (in meters)
#distance between [leicester,coventry]=52.95513,[Birmingham,coventry]=55.26608 (in kms)
#so exp(-52.95513) = 1.004337e-23 ; exp(-55.26608) = 9.959704e-25


infected_second_category<-infected_second_simulation[[1]][1:training_upto]
infected_third_category<-infected_first_simulation[[1]][1:training_upto]

log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-W[3,2]*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<- W[3,1]*log_infected_third_category       #this is w[i,j]*infected_3rd_category

infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     
                                                     "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))


# birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
# xn<-as.numeric(birth_data[15, 3:((training_upto/26)+2)])
# divided_vector <- xn / 26
# biweekly_birth_data<- rep(divided_vector, each = 26)
# length(biweekly_birth_data)
# biweekly_birth_data_in_inverselogit1<- biweekly_birth_data[-length(biweekly_birth_data)]
# #biweekly_birth_data_in_inverselogit<- biweekly_birth_data_in_inverselogit[-c(1:(ar_order-1))]
# if (ar_order == 1) {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1
# } else {
#   biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1[-c(1:(ar_order-1))]
# }
# length(biweekly_birth_data_in_inverselogit)
t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)

# baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
# baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
# if (ar_order == 1) {
#   baby_boom_effect <- baby_boom_effect1
# } else {
#   baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
# }
# length(baby_boom_effect)






head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)

#length(biweekly_birth_data_in_inverselogit)
length(seasonality1)
length(seasonality2)
length(t)
#length(baby_boom_effect)

length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)

nrow(data_sf)

# data<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),
#                infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
#                biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect
# ),
# nrow = length(infected_first_category_in_inverselogit$lag1))

X_t <- X3[-c(1:(ar_order-1))]
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  seasonality2,t,X_t)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","time","birth_cohort")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)

# 
# 
# acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.50)
# pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.50)
# 
# 
predicted3<-model$fitted.values*(data$S+data$F)
# data.frame( predicted,infected_first_category_success)


plot(infected_first_category_success, type = "l", col = "black", ylim = c(0, max(infected_first_category, predicted3)),xlab = "Time",ylab = "Counts")
lines(predicted3, type = "l", col = "red")  # Add the line for predicted_infection
points(predicted3, col = "red", pch = 8, cex = .25)  # Add the star-shaped points
legend("topright",legend = c("Real","Estimated"),col=c("black","red"),lty =c(1,1),cex=.4 )


#################################################################################################################
#################################################################################################################
predicted<-as.numeric(c(0,predicted))
length(predicted)

if (length(predicted) != nstep) {
  stop("Length of predicted vector does not match nstep!")
}


predicted2<-ceiling(as.numeric(c(0,predicted2)))
length(predicted2)

if (length(predicted2) != nstep) {
  stop("Length of predicted vector does not match nstep!")
}

predicted3<-ceiling(as.numeric(c(0,predicted3)))
length(predicted3)

if (length(predicted3) != nstep) {
  stop("Length of predicted vector does not match nstep!")
}


time <- seq(1, nstep, by = 1)
# Simulate some disease counts for three regions
region1_count <- predicted
region2_count <- predicted2
region3_count <- predicted3

# Create region IDs
region <- factor(c(rep(1, nstep), rep(2, nstep), rep(3, nstep)))


print(length(region1_count))  # Should be nstep
print(length(region2_count))  # Should be nstep
print(length(region3_count))  # Should be nstep


data <- data.frame(
  time = rep(time, 3), 
  region = region, 
  count = c(region1_count, region2_count, region3_count)
)

# plot_ly(data, 
#         x = ~time, 
#         y = ~region, 
#         z = ~count, 
#         type = 'scatter3d', 
#         mode = 'lines', 
#         color = ~region, 
#         line = list(width = 4)) %>%
#   layout(scene = list(
#     xaxis = list(title = 'Time (in biweeks)'),
#     yaxis = list(title = 'Regions', tickvals = c(1, 2, 3), ticktext = c('Region-1', 'Region-2', 'Region-3')),
#     zaxis = list(title = 'Disease Count', range = c(0, 1500))
#   ))


#################################################################################################################
#################################################################################################################

#Without box


plot_ly(colors = c("red", "blue", "green")) %>%
  add_trace(data = data, 
            x = ~time, 
            y = ~region, 
            z = ~count, 
            type = 'scatter3d', 
            mode = 'lines', 
            color = as.factor(data$region),  # Convert to factor to ensure proper labeling
            name = ~paste0("Region ", region, " - Predicted Count"), # Assign custom names
            line = list(width = 4)) %>%
  add_trace(
    x = time, 
    y = rep(1, nstep),  # Overlay on Region-1
    z = infected_first_simulation[[1]], 
    type = 'scatter3d', 
    mode = 'markers',
    name = 'Region 1 - Actual Count',  # Rename trace in legend
    marker = list(
      color = 'red', 
      size = 1,  
      symbol = 'star'  
    )
  ) %>%
  add_trace(
    x = time, 
    y = rep(2, nstep),  # Region-2
    z = infected_second_simulation[[1]], 
    type = 'scatter3d', 
    mode = 'markers',  
    name = 'Region 2 - Actual Count',  # Rename trace in legend
    marker = list(
      color = 'blue',  
      size = 2,  
      symbol = 'star'
    )
  ) %>%
  add_trace(
    x = time, 
    y = rep(3, nstep),  # Region-3
    z = infected_third_simulation[[1]], 
    type = 'scatter3d', 
    mode = 'markers',  
    name = 'Region 3 - Actual Count',  # Rename trace in legend
    marker = list(
      color = 'green',  
      size = 2,  
      symbol = 'star'
    )
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Time (in biweeks)', showgrid = FALSE,
                   titlefont = list(size = 17)),
      yaxis = list(title = 'Regions', tickvals = c(1, 2, 3), ticktext = c('Region-1', 'Region-2', 'Region-3'), showgrid = FALSE,
                   titlefont = list(size = 17)),
      zaxis = list(title = 'Disease Count', range = c(0, 1500), showgrid = FALSE, showgrid = FALSE,
                   tickvals = seq(200, 1400, by = 200),
                   titlefont = list(size = 17))
    ),showlegend = FALSE,
    # Move legend inside the top-right corner
    shapes = list(
      list(
        type = "rect",
        x0 = 0, x1 = 1, y0 = 0, y1 = 1,  # Cover entire plot area
        xref = "paper", yref = "paper",  # Use relative positioning
        line = list(color = "black", width = 0.5)  # Define border color & thickness
      )
    )
  )




#################################################################################################################
#################################################################################################################

#With box

# library(plotly)
# 
# plot_ly() %>%
#   add_trace(data = data, 
#             x = ~time, 
#             y = ~region, 
#             z = ~count, 
#             type = 'scatter3d', 
#             mode = 'lines', 
#             color = as.factor(data$region),  
#             name = ~paste0("Region ", region, " - Actual Count"),  
#             line = list(width = 4)) %>%
#   add_trace(
#     x = time, 
#     y = rep(1, nstep),  
#     z = infected_first_simulation[[1]], 
#     type = 'scatter3d', 
#     mode = 'markers',
#     name = 'Region 1 - Actual Count',  
#     marker = list(
#       color = 'black', 
#       size = 1,  
#       symbol = 'star'  
#     )
#   ) %>%
#   add_trace(
#     x = time, 
#     y = rep(2, nstep),  
#     z = infected_second_simulation[[1]], 
#     type = 'scatter3d', 
#     mode = 'markers',  
#     name = 'Region 2 - Actual Count',  
#     marker = list(
#       color = 'black',  
#       size = 2,  
#       symbol = 'star'
#     )
#   ) %>%
#   add_trace(
#     x = time, 
#     y = rep(3, nstep),  
#     z = infected_third_simulation[[1]], 
#     type = 'scatter3d', 
#     mode = 'markers',  
#     name = 'Region 3 - Actual Count',  
#     marker = list(
#       color = 'black',  
#       size = 2,  
#       symbol = 'star'
#     )
#   ) %>%
#   layout(
#     scene = list(
#       xaxis = list(title = 'Time (in biweeks)', showgrid = FALSE),
#       yaxis = list(title = 'Regions', tickvals = c(1, 2, 3), ticktext = c('Region-1', 'Region-2', 'Region-3'), showgrid = FALSE),
#       zaxis = list(title = 'Disease Count', range = c(0, 1500), showgrid = FALSE)
#     ),
#     legend = list(x = 0.6, y = 0.9,  # Move legend inside the plot (top-right)
#                   font = list(size = 6)),
#     #Box around the plot using shapes
#     shapes = list(
#       list(
#         type = "rect",
#         x0 = 0, x1 = 1, y0 = 0, y1 = 1,  # Cover entire plot area
#         xref = "paper", yref = "paper",  # Use relative positioning
#         line = list(color = "black", width = 0.5)  # Define border color & thickness
#       )
#     )
#   )
