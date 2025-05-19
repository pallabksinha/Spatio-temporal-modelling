#Simulation

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
beta1<- .5 ; beta2<- beta1

nstep<-26*6

B <- 10^3
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



ar_order<-2
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


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
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1),
                                                     "lag2_2nd"=dplyr::lag(log_infected_second_category, 2),
                                                     "lag2_3rd"=dplyr::lag(log_infected_third_category, 2))

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag2_3rd")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))



t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)



head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)


length(seasonality1)
length(seasonality2)
length(t)


length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)
length(infected_first_category_in_inverselogit$lag2_2nd)
length(infected_first_category_in_inverselogit$lag2_3rd)

nrow(data_sf)


X_t <- X3[-c(1:(ar_order))]
length(X_t)
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  seasonality2,seasonality1,t,X_t,infected_first_category_in_inverselogit$lag2_2nd,infected_first_category_in_inverselogit$lag2_3rd)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","Seasonalit1","time","birth_cohort","lag2_2nd","lag2_3rd")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+Seasonalit1+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)
predicted<-model$fitted.values*(data$S+data$F)


predicted3<-model$fitted.values*(data$S+data$F)

#############################################3
n_trials <- rowSums(data_sf)
y <- infected_first_category_success
length(y)

predict_prob<- model$fitted.values


residuals_rq <- numeric(length(y))
for (i in 1:length(y)) {
  
  # CDF for y - 1 (lower bound)
  p_lower <- pbinom(round(y[i] - 1), size = round(n_trials[i]), prob = predict_prob[i])
  
  # CDF for y (upper bound)
  p_upper <- pbinom(round(y[i]), size = round(n_trials[i]), prob = predict_prob[i])
  
  # Generate a random uniform value between the lower and upper CDF
  u <- runif(1, p_lower, p_upper)
  
  # Transform to the normal quantile
  residuals_rq[i] <- qnorm(u)
}



residuals_rq<-residuals_rq[is.finite(residuals_rq)]
residuals_rq<- scale(residuals_rq)


ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )
ks.test(residuals_rq, "pnorm",mean = 0, sd = 1 )

Box_ljung_test_rq<- Box.test(residuals_rq,lag=500,type = "Ljung-Box")
Box_ljung_test_rq

############################################################################################
############################################################################################
############################################################################################

par(mfrow = c(2, 2))
acf(residuals_rq ,xlab="Lag",main="ACF plot",ci=0.995,ylim=c(-1,1))
#acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.99)
pacf(residuals_rq,   xlab="Lag",main="PACF plot",ci=0.995,ylim=c(-1,1))
#pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.99)
plot(residuals_rq,ylim = c(-3,3),ylab = "Residuals",main="Residuals")
abline(h=0,col="red",lty=3)
# qqnorm(residuals_rq, main = "QQ Plot for Normality")
# abline(0, 1, col = "red", lwd = 0.2) 
qqnorm(residuals_rq, main = "QQ Plot for Normality")
qqline(residuals_rq, col = "red", lwd = 0.5)  # Add a reference line
par(mfrow = c(1, 1))


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
infected_second_category<- infected_second_simulation[[1]][1:training_upto]
infected_third_category<-infected_third_simulation[[1]][1:training_upto]



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



ar_order<-2
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

# t<-1:length(infected_first_category)
#infected_first_category_minus_last_term<-data$Bath[-nstep]  #prob. (ie pie) ber korar somoi I1,..,,I597 lagbe.. pie1 is functn of I1,.., pie2 is functn of I2 etc.. so last er ta delete kore dewoa holo

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-W[1,2]*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<- W[1,3]*log_infected_third_category       #this is w[i,j]*infected_3rd_category

infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2),
                                                     "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1),
                                                     "lag2_2nd"=dplyr::lag(log_infected_second_category, 2),
                                                     "lag2_3rd"=dplyr::lag(log_infected_third_category, 2))

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag2_3rd")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))



t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)



head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)


length(seasonality1)
length(seasonality2)
length(t)


length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)
length(infected_first_category_in_inverselogit$lag2_2nd)
length(infected_first_category_in_inverselogit$lag2_3rd)

nrow(data_sf)

# data<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),
#                infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
#                biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect
# ),
# nrow = length(infected_first_category_in_inverselogit$lag1))

X_t <- X3[-c(1:(ar_order))]
length(X_t)
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag2+infected_first_category_in_inverselogit$lag1_2nd+
                    infected_first_category_in_inverselogit$lag1_3rd+infected_first_category_in_inverselogit$lag2_2nd+infected_first_category_in_inverselogit$lag2_3rd,
                  seasonality2,seasonality1,t,X_t)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag2_and_Lag1_2nd_and_LAg1_3rd_and_lag2_2nd_and_lag2_3rd","Seasonalit2","Seasonalit1","time","birth_cohort")
head(data,20)

model <- glm(cbind(S,F) ~Lag1_and_Lag2_and_Lag1_2nd_and_LAg1_3rd_and_lag2_2nd_and_lag2_3rd+Seasonalit2+Seasonalit1+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)
predicted<-model$fitted.values*(data$S+data$F)
data.frame( predicted,infected_first_category_success)

predicted3<-model$fitted.values*(data$S+data$F)


#############################################3
n_trials <- rowSums(data_sf)
y <- infected_first_category_success
length(y)

predict_prob<- model$fitted.values

residuals_rq <- numeric(length(y))
for (i in 1:length(y)) {
  
  # CDF for y - 1 (lower bound)
  p_lower <- pbinom(round(y[i] - 1), size = round(n_trials[i]), prob = predict_prob[i])
  
  # CDF for y (upper bound)
  p_upper <- pbinom(round(y[i]), size = round(n_trials[i]), prob = predict_prob[i])
  
  # Generate a random uniform value between the lower and upper CDF
  u <- runif(1, p_lower, p_upper)
  
  # Transform to the normal quantile
  residuals_rq[i] <- qnorm(u)
}

residuals_rq<-residuals_rq[is.finite(residuals_rq)]
residuals_rq<- scale(residuals_rq)

ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )
ks.test(residuals_rq, "pnorm",mean = 0, sd = 1 )

Box_ljung_test_rq<- Box.test(residuals_rq,lag=500,type = "Ljung-Box")
Box_ljung_test_rq

############################################################################################
############################################################################################
############################################################################################

par(mfrow = c(2, 2))
acf(residuals_rq ,xlab="Lag",main="ACF plot",ci=0.995,ylim=c(-1,1))
#acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.99)
pacf(residuals_rq,   xlab="Lag",main="PACF plot",ci=0.995,ylim=c(-1,1))
#pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.99)
plot(residuals_rq,ylim = c(-3,3),ylab = "Residuals",main="Residuals")
abline(h=0,col="red",lty=3)
# qqnorm(residuals_rq, main = "QQ Plot for Normality")
# abline(0, 1, col = "red", lwd = 0.2) 
qqnorm(residuals_rq, main = "QQ Plot for Normality")
qqline(residuals_rq, col = "red", lwd = 0.5)  # Add a reference line
par(mfrow = c(1, 1))

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


infected_first_category <- infected_second_simulation[[j]]
infected_second_category<- infected_first_simulation[[j]]
infected_third_category<-infected_third_simulation[[j]]

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



ar_order<-3
infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke I598 obdhi lagbe . So first one deleted
length(infected_first_category_success)

# t<-1:length(infected_first_category)
#infected_first_category_minus_last_term<-data$Bath[-nstep]  #prob. (ie pie) ber korar somoi I1,..,,I597 lagbe.. pie1 is functn of I1,.., pie2 is functn of I2 etc.. so last er ta delete kore dewoa holo

log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-W[2,1]*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<- W[2,3]*log_infected_third_category       #this is w[i,j]*infected_3rd_category

infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2),
                                                     "lag3"=dplyr::lag(log_infected_first_category, 3),
                                                     "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                     "lag1_3rd"=dplyr::lag(log_infected_third_category, 1),
                                                     "lag2_2nd"=dplyr::lag(log_infected_second_category, 2),
                                                     "lag2_3rd"=dplyr::lag(log_infected_third_category, 2),
                                                     "lag3_2nd"=dplyr::lag(log_infected_second_category, 3),
                                                     "lag3_3rd"=dplyr::lag(log_infected_third_category, 3)
)

infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag3")


head(infected_first_category_in_inverselogit,20)

cum_fromI2<-(cumsum(infected_first_category[-1]))
infected_first_category_failure<-S1-cum_fromI2
if(ar_order>1){infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]}

length(infected_first_category_failure)

data_sf<-cbind(infected_first_category_success,infected_first_category_failure)
data_sf<-matrix(c(infected_first_category_success,infected_first_category_failure),nrow = length(infected_first_category_success))



t<- (ar_order+1):training_upto
length(t)
seasonality1<-cos(2*pi*t/52)
seasonality2<-cos(2*pi*t/26)
length(seasonality1)



head(infected_first_category_in_inverselogit,n=20)
#length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)


length(seasonality1)
length(seasonality2)
length(t)


length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)
length(infected_first_category_in_inverselogit$lag2_2nd)
length(infected_first_category_in_inverselogit$lag2_3rd)

nrow(data_sf)

X_t <- X3[-c(1:(ar_order))]
length(X_t)
data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  seasonality2,seasonality1,t,X_t,infected_first_category_in_inverselogit$lag2_2nd,infected_first_category_in_inverselogit$lag2_3rd)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","Seasonalit2","Seasonalit1","time","birth_cohort","lag2_2nd","lag2_3rd")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+Seasonalit2+Seasonalit1+time+birth_cohort, 
             family = binomial(link = "logit"), data = data)

summary(model)
predicted<-model$fitted.values*(data$S+data$F)


predicted3<-model$fitted.values*(data$S+data$F)


#############################################3
n_trials <- rowSums(data_sf)
y <- infected_first_category_success
length(y)

predict_prob<- model$fitted.values

residuals_rq <- numeric(length(y))
for (i in 1:length(y)) {
  
  # CDF for y - 1 (lower bound)
  p_lower <- pbinom(round(y[i] - 1), size = round(n_trials[i]), prob = predict_prob[i])
  
  # CDF for y (upper bound)
  p_upper <- pbinom(round(y[i]), size = round(n_trials[i]), prob = predict_prob[i])
  
  # Generate a random uniform value between the lower and upper CDF
  u <- runif(1, p_lower, p_upper)
  
  # Transform to the normal quantile
  residuals_rq[i] <- qnorm(u)
}

residuals_rq<-residuals_rq[is.finite(residuals_rq)]
residuals_rq<- scale(residuals_rq)

ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )
ks.test(residuals_rq, "pnorm",mean = 0, sd = 1 )

Box_ljung_test_rq<- Box.test(residuals_rq,lag=500,type = "Ljung-Box")
Box_ljung_test_rq

############################################################################################
############################################################################################
############################################################################################

par(mfrow = c(2, 2))
acf(residuals_rq ,xlab="Lag",main="ACF plot",ci=0.995,ylim=c(-1,1))
#acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.99)
pacf(residuals_rq,   xlab="Lag",main="PACF plot",ci=0.995,ylim=c(-1,1))
#pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.99)
plot(residuals_rq,ylim = c(-3,3),ylab = "Residuals",main="Residuals")
abline(h=0,col="red",lty=3)
# qqnorm(residuals_rq, main = "QQ Plot for Normality")
# abline(0, 1, col = "red", lwd = 0.2) 
qqnorm(residuals_rq, main = "QQ Plot for Normality")
qqline(residuals_rq, col = "red", lwd = 0.5)  # Add a reference line
par(mfrow = c(1, 1))





