rm(list = ls())
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
training_upto<-546
infected_first_category <- data$Coventry[1:training_upto]
plot(infected_first_category,type = "l",col="red",ylab="Case Count",xlab = "Time(biweekly)")
acf(infected_first_category,xlab="Lag",main="Acf plot for counts",ci=0.99)
pacf(infected_first_category,   xlab="Lag",main="PAcf plot for counts",ci=0.99)

length(infected_first_category)
nstep<-length(infected_first_category)

N<-265950
S1<-265950

cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
length(cum)
cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
length(cum)
cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
length(cum)



ar_order<-1
infected_first_category_success<-infected_first_category[-c(1:ar_order)] 
log_infected_first_category<-infected_first_category
log_infected_first_category[log_infected_first_category == 0] <- 1
log_infected_first_category<-log(log_infected_first_category)
infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                     "lag2"=dplyr::lag(log_infected_first_category, 2) )


#distance between [leicester,coventry]=52955.13,[Birmingham,coventry]=55266.08 (in meters)
#distance between [leicester,coventry]=52.95513,[Birmingham,coventry]=55.26608 (in kms)
#so exp(-52.95513) = 1.004337e-23 ; exp(-55.26608) = 9.959704e-25


infected_second_category<-data$Birmingham[1:training_upto]
infected_third_category<-data$Leicester[1:training_upto]

log_infected_second_category<-infected_second_category
log_infected_second_category[log_infected_second_category == 0] <- 1
log_infected_second_category<-log(log_infected_second_category)
log_infected_second_category<-9.959704e-25*log_infected_second_category    #this is w[i,j]*infected_2nd_category

log_infected_third_category<-infected_third_category
log_infected_third_category[log_infected_third_category == 0] <- 1
log_infected_third_category<-log(log_infected_third_category)
log_infected_third_category<-1.004337e-23*log_infected_third_category       #this is w[i,j]*infected_3rd_category

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


birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
xn<-as.numeric(birth_data[15, 3:((training_upto/26)+2)])
divided_vector <- xn / 26
biweekly_birth_data<- rep(divided_vector, each = 26)
length(biweekly_birth_data)
biweekly_birth_data_in_inverselogit1<- biweekly_birth_data[-length(biweekly_birth_data)]
#biweekly_birth_data_in_inverselogit<- biweekly_birth_data_in_inverselogit[-c(1:(ar_order-1))]
if (ar_order == 1) {
  biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1
} else {
  biweekly_birth_data_in_inverselogit <- biweekly_birth_data_in_inverselogit1[-c(1:(ar_order-1))]
}
length(biweekly_birth_data_in_inverselogit)
t<- (ar_order+1):training_upto
length(t)
seasonality1<-sin(2*pi*t/52)
seasonality2<-sin(2*pi*t/26)
length(seasonality1)

baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
if (ar_order == 1) {
  baby_boom_effect <- baby_boom_effect1
} else {
  baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
}
length(baby_boom_effect)






head(infected_first_category_in_inverselogit,n=20)
length(biweekly_birth_data_in_inverselogit)
length(infected_first_category_in_inverselogit$lag1)
length(infected_first_category_in_inverselogit$lag2)

length(biweekly_birth_data_in_inverselogit)
length(seasonality1)
length(seasonality2)
length(t)
length(baby_boom_effect)

length(infected_first_category_in_inverselogit$lag1_2nd)
length(infected_first_category_in_inverselogit$lag1_3rd)

nrow(data_sf)

data<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),
               infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
               biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect
),
nrow = length(infected_first_category_in_inverselogit$lag1))

data<- data.frame(infected_first_category_in_inverselogit$lag1+infected_first_category_in_inverselogit$lag1_2nd+infected_first_category_in_inverselogit$lag1_3rd,
                  biweekly_birth_data_in_inverselogit,seasonality2,t,baby_boom_effect)
data<-cbind(data_sf,data)

colnames(data) <- c("S","F","Lag1_and_Lag1_2nd_and_LAg1_3rd","biweeklybirth","Seasonalit2","time","baby_boom")
head(data,20)

model <- glm(cbind(S,F) ~ Lag1_and_Lag1_2nd_and_LAg1_3rd+biweeklybirth+Seasonalit2+time+baby_boom, 
             family = binomial(link = "logit"), data = data)


acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.50)
pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.50)


predicted<-model$fitted.values*(data$S+data$F)
data.frame( predicted,infected_first_category_success)


plot(infected_first_category_success, type = "l", col = "black", ylim = c(0, max(infected_first_category, predicted)),xlab = "Time",ylab = "Counts")
lines(predicted, type = "l", col = "red")  # Add the line for predicted_infection
points(predicted, col = "red", pch = 8, cex = .25)  # Add the star-shaped points
legend("topright",legend = c("Real","Estimated"),col=c("black","red"),lty =c(1,1),cex=.4 )

Box_ljung_test<- Box.test(model$residuals,lag=500,type = "Ljung-Box")
Box_ljung_test
