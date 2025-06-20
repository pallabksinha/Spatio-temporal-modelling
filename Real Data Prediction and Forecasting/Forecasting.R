#Coventry

rm(list = ls())
set.seed(25)
library(readxl)
library(glarma)
library(psych)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
N<-265950
S1<-265950
ar_order<-2
birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
xn<-as.numeric(birth_data[15, 3:  25])
divided_vector <- xn / 26
biweekly_birth_data<- rep(divided_vector, each = 26)


start_value_of_k<-50
end_value_of_k<-100


forecasting_loop_starts_from <-571
forecasting_loop_ends_in <- 595

probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)    #"probabilities" are the predicted values
predicted_probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)   #"predicted_probabilities" are probabilities of "probabilities" 
predicted_probabilities_unlist<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)
stored<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)  #stored will be unlisting "probabilities"

no_iterations<- 1000000  #for HDR intervals
HDR_interval_pts <- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_probabilities<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_pts_to_keep<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)

lower_interval_pts<-c()
upper_interval_pts<-c()

#glarmamod<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,list(NA),simplify = F)
glarmamod<- replicate(end_value_of_k-start_value_of_k+1,list(NA),simplify = F)
alphaa<-.05

start_time <- Sys.time() 
for(training_upto in forecasting_loop_starts_from:forecasting_loop_ends_in){
  infected_first_category <- data$Coventry[1:(training_upto+1)]
  infected_second_category<-data$Birmingham[1:(training_upto+1)]
  infected_third_category<-data$Leicester[1:(training_upto+1)]
  
  nstep<-length(infected_first_category)
  
  cum<-cumsum(infected_first_category[-1]) 
  cum<-c(S1,S1-cum) 
  cum<-cum[-nstep]   
  
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom / x.. binomial e I2 theke..
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]
  
  data_sf<-cbind(infected_first_category_success[-c(length(infected_first_category_success))],
                 infected_first_category_failure[-c(length(infected_first_category_success))])
  data_sf<-matrix(c(infected_first_category_success[-c(length(infected_first_category_success))],
                    infected_first_category_failure[-c(length(infected_first_category_success))]),
                  nrow = length(infected_first_category_success[-c(length(infected_first_category_success))]))
  
  
  biweekly_birth_data_in_inverselogit<- biweekly_birth_data[ar_order:(training_upto)]
  t<- (ar_order+1): (training_upto+1)
  seasonality2<-cos(2*pi*t/26)
  
  for(k in start_value_of_k:end_value_of_k){
    
    new_row <- c(k, S1-(sum(data_sf[,1])+k) )
    
    temp_data_sf <- rbind(data_sf, new_row)
    
    log_infected_first_category<-infected_first_category
    log_infected_first_category<-data$Coventry[1:(training_upto+1)]
    log_infected_first_category<-log_infected_first_category[-c(length(log_infected_first_category))]
    log_infected_first_category<- c(log_infected_first_category,k)
    log_infected_first_category[log_infected_first_category == 0] <- 1
    log_infected_first_category<-log(log_infected_first_category)
    
    log_infected_second_category<-infected_second_category
    log_infected_second_category<-data$Birmingham[1:(training_upto+1)]
    log_infected_second_category<-log_infected_second_category[-c(length(log_infected_second_category))]
    log_infected_second_category<- c(log_infected_second_category,k)
    log_infected_second_category[log_infected_second_category == 0] <- 1
    log_infected_second_category<-log(log_infected_second_category)
    log_infected_second_category<-1.809428e-05*log_infected_second_category
    
    log_infected_third_category<-infected_third_category
    log_infected_third_category<-data$Leicester[1:(training_upto+1)]
    log_infected_third_category<-log_infected_third_category[-c(length(log_infected_third_category))]
    log_infected_third_category<- c(log_infected_third_category,k)
    log_infected_third_category[log_infected_third_category == 0] <- 1
    log_infected_third_category<-log(log_infected_third_category)
    log_infected_third_category<-1.888391e-05*log_infected_third_category 
    
    infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),"lag2"=dplyr::lag(log_infected_first_category, 2),
                                                         
                                                         "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                         "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
    
    infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag2")
    
    data_covariates<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag2)),infected_first_category_in_inverselogit$lag1,
                              infected_first_category_in_inverselogit$lag2,
                              infected_first_category_in_inverselogit$lag1_2nd,
                              infected_first_category_in_inverselogit$lag1_3rd,
                              biweekly_birth_data_in_inverselogit,seasonality2,t),
                            nrow = length(infected_first_category_in_inverselogit$lag2))
    colnames(data_covariates) <- c("Intercept","Lag1","Lag2","Lag1_2nd","Lag1_3rd","biweeklybirth","Seasonalit","time")
    
    glarmamod[[k-(start_value_of_k-1)]] <- glarma(temp_data_sf, data_covariates,thetaLags = c(1,2), type = "Bin", method = "FS",
                                                  residuals = "Identity", maxit = 100, grad = 1e-6)
    probabilities[[k-(start_value_of_k-1)]]<- glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]
    predicted_probabilities[[k-(start_value_of_k-1)]]<- dbinom(x=floor(glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]),
                                                               size = S1-(sum(data_sf[,1])+k),
                                                               prob = logistic(as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[3])*infected_first_category_in_inverselogit$lag2[length(infected_first_category_in_inverselogit$lag2)]+
                                                                                 
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[4])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+  
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[5])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[6])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[7])*seasonality2[length(seasonality2)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[8])*t[length(t)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[9])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[10])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-2])
                                                               )   )
  }
  
  predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]]<- unlist(predicted_probabilities)[!is.na(unlist(predicted_probabilities))]
  #stored[[training_upto-(forecasting_loop_starts_from-1)]]<-median(unlist(probabilities)[!is.na(unlist(probabilities))])
  stored[[training_upto-(forecasting_loop_starts_from-1)]]<-unlist(probabilities)[!is.na(unlist(probabilities))][which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])]
  
  max_position<- which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])
  
  glarmamod1<-glarmamod[[max_position]]    
  # So there are M possible candidate values for the time point we are trying to forecast.. among them we are choosing the one which produces the highest probability
  # consequently in HDR also, we will use this set of estimated parameters 
  
  HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]]<- rbinom(n=no_iterations,size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                              logistic(as.numeric(glarmamod1$delta[1])+
                                                                                         as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                         as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag2[length(infected_first_category_in_inverselogit$lag2)]+
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                         as.numeric(glarmamod1$delta[5])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[6])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                         as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                         as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                         as.numeric(glarmamod1$delta[9])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                         as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                       
                                                                              ) )
  
  HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]<- dbinom(x=HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]],
                                                                                        size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                                        prob=logistic(as.numeric(glarmamod1$delta[1])+
                                                                                                        as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                                        as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag2[length(infected_first_category_in_inverselogit$lag2)]+
                                                                                                        
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                                        as.numeric(glarmamod1$delta[5])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[6])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                                        as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                                        as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                                        as.numeric(glarmamod1$delta[9])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                                        as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                                      
                                                                                        ) )
  
  
  HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]]<-  
    HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]][which(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]>=
                                                                               quantile(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]],alphaa))]
  lower_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- min(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  upper_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- max(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  
}


end_time <- Sys.time()   
end_time - start_time

HDR_interval_pts
HDR_interval_pts_to_keep

lower_interval_pts
upper_interval_pts
stored<-floor(unlist(stored)[!is.na(unlist(stored))])
stored<-pmax(lower_interval_pts, pmin(stored, upper_interval_pts))

x<-data.frame( "Lower HDR"=lower_interval_pts,"Upper HDR"=upper_interval_pts,"Predicted"=stored, "True value"=data$Coventry[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])
x$"inclusion"<- ifelse(x$True.value<x$Upper.HDR & x$True.value>x$Lower.HDR,1,0)
x

par(mfrow=c(3,1))

plot(x$Predicted,type = "p",ylim = c(min(lower_interval_pts,stored,data$Coventry[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)]),
                                     max(upper_interval_pts,stored,data$Coventry[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])),ylab="HDR Interval",xlab = "Time(biweekly)",main = "95% HDR Interval for Coventry Data"
)
lines(x$True.value,type = "p",col="darkred",pch=16)
# lines(lower_interval_pts,type = "p",pch=3)
lines(lower_interval_pts,type = "l",lty=4)
lines(upper_interval_pts,type = "l",lty=4)
# legend(x = "topright",lty = c(NA,NA,4,4),pch=c(16,1,NA,NA),
#        col= c("darkred","black","black","black"),
#        legend=c( "observed","Predicted","HDR Interval"),cex=0.1)


#########################################################################################################################

#Leicester

set.seed(25)
library(readxl)
library(glarma)
library(psych)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
N<-276550
S1<-276550
ar_order<-1
birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
xn<-as.numeric(birth_data[15, 3:  25])
divided_vector <- xn / 26
biweekly_birth_data<- rep(divided_vector, each = 26)


start_value_of_k<-50
end_value_of_k<-100


forecasting_loop_starts_from <-571
forecasting_loop_ends_in <- 595

probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)    #"probabilities" are the predicted values
predicted_probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)   #"predicted_probabilities" are probabilities of "probabilities" 
predicted_probabilities_unlist<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)
stored<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)  #stored will be unlisting "probabilities"

no_iterations<- 1000000  #for HDR intervals
HDR_interval_pts <- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_probabilities<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_pts_to_keep<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)

lower_interval_pts<-c()
upper_interval_pts<-c()

#glarmamod<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,list(NA),simplify = F)
glarmamod<- replicate(end_value_of_k-start_value_of_k+1,list(NA),simplify = F)
alphaa<-.05

start_time <- Sys.time() 
for(training_upto in forecasting_loop_starts_from:forecasting_loop_ends_in){
  infected_first_category <- data$Leicester[1:(training_upto+1)]
  infected_second_category<-data$Birmingham[1:(training_upto+1)]
  infected_third_category<-data$Coventry[1:(training_upto+1)]
  
  nstep<-length(infected_first_category)
  
  cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I548
  cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
  cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
  
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke..
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]
  
  data_sf<-cbind(infected_first_category_success[-c(length(infected_first_category_success))],
                 infected_first_category_failure[-c(length(infected_first_category_success))])
  data_sf<-matrix(c(infected_first_category_success[-c(length(infected_first_category_success))],
                    infected_first_category_failure[-c(length(infected_first_category_success))]),
                  nrow = length(infected_first_category_success[-c(length(infected_first_category_success))]))
  #two deleted because 547 no entry random dhukbe and 548 no entry ta predict hobe .. 
  
  biweekly_birth_data_in_inverselogit<- biweekly_birth_data[ar_order:(training_upto)]
  t<- (ar_order+1): (training_upto+1)
  
  seasonality1<-cos(2*pi*t/52)
  seasonality2<-cos(2*pi*t/26)
  baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
  baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
  if (ar_order == 1) {
    baby_boom_effect <- baby_boom_effect1
  } else {
    baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
  }
  
  
  for(k in start_value_of_k:end_value_of_k){
    
    new_row <- c(k, S1-(sum(data_sf[,1])+k) )
    
    temp_data_sf <- rbind(data_sf, new_row)
    
    log_infected_first_category<-infected_first_category
    log_infected_first_category<-data$Leicester[1:(training_upto+1)]
    log_infected_first_category<-log_infected_first_category[-c(length(log_infected_first_category))]
    log_infected_first_category<- c(log_infected_first_category,k)
    log_infected_first_category[log_infected_first_category == 0] <- 1
    log_infected_first_category<-log(log_infected_first_category)
    
    log_infected_second_category<-infected_second_category
    log_infected_second_category<-data$Birmingham[1:(training_upto+1)]
    log_infected_second_category<-log_infected_second_category[-c(length(log_infected_second_category))]
    log_infected_second_category<- c(log_infected_second_category,k)
    log_infected_second_category[log_infected_second_category == 0] <- 1
    log_infected_second_category<-log(log_infected_second_category)
    # log_infected_second_category<-1.809428e-05*log_infected_second_category
    
    log_infected_third_category<-infected_third_category
    log_infected_third_category<-data$Coventry[1:(training_upto+1)]
    log_infected_third_category<-log_infected_third_category[-c(length(log_infected_third_category))]
    log_infected_third_category<- c(log_infected_third_category,k)
    log_infected_third_category[log_infected_third_category == 0] <- 1
    log_infected_third_category<-log(log_infected_third_category)
    # log_infected_third_category<-1.888391e-05*log_infected_third_category 
    
    infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                         
                                                         "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                         "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
    
    infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")
    
    data_covariates<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),infected_first_category_in_inverselogit$lag1,
                              
                              infected_first_category_in_inverselogit$lag1_2nd,
                              infected_first_category_in_inverselogit$lag1_3rd,
                              biweekly_birth_data_in_inverselogit,seasonality1,seasonality2,t,baby_boom_effect),
                            nrow = length(infected_first_category_in_inverselogit$lag1))
    colnames(data_covariates) <- c("Intercept","Lag1","Lag1_2nd","Lag1_3rd","biweeklybirth","Seasonalit1","Seasonalit2","time","baby_boom_effect")
    
    glarmamod[[k-(start_value_of_k-1)]] <- glarma(temp_data_sf, data_covariates,thetaLags = c(1,2), type = "Bin", method = "FS",
                                                  residuals = "Identity", maxit = 100, grad = 1e-6)
    probabilities[[k-(start_value_of_k-1)]]<- glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]
    predicted_probabilities[[k-(start_value_of_k-1)]]<- dbinom(x=floor(glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]),
                                                               size = S1-(sum(data_sf[,1])+k),
                                                               prob = logistic(as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+  
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[6])*seasonality1[length(seasonality1)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[7])*seasonality2[length(seasonality2)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[8])*t[length(t)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[10])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[11])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-2])
                                                               )   )
  }
  
  predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]]<- unlist(predicted_probabilities)[!is.na(unlist(predicted_probabilities))]
  #stored[[training_upto-(forecasting_loop_starts_from-1)]]<-median(unlist(probabilities)[!is.na(unlist(probabilities))])
  stored[[training_upto-(forecasting_loop_starts_from-1)]]<-unlist(probabilities)[!is.na(unlist(probabilities))][which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])]
  
  max_position<- which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])
  # i <- (max_position - 1) %/% length(glarmamod[[1]]) + 1  # List index
  # j <- (max_position - 1) %% length(glarmamod[[1]]) + 1    # Sub-index within each list element
  glarmamod1<-glarmamod[[max_position]]    
  
  
  HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]]<- rbinom(n=no_iterations,size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                              logistic(as.numeric(glarmamod1$delta[1])+
                                                                                         as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                         
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                         as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                         as.numeric(glarmamod1$delta[6])*seasonality1[length(seasonality1)]+
                                                                                         as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                         as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                         as.numeric(glarmamod1$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                         as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                         as.numeric(glarmamod1$delta[11])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                       
                                                                              ) )
  
  HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]<- dbinom(x=HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]],
                                                                                        size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                                        prob=logistic(as.numeric(glarmamod1$delta[1])+
                                                                                                        as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                                        
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                                        as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                                        as.numeric(glarmamod1$delta[6])*seasonality1[length(seasonality1)]+
                                                                                                        as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                                        as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                                        as.numeric(glarmamod1$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                                        as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                                        as.numeric(glarmamod1$delta[11])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                                      
                                                                                        ) )
  
  
  HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]]<-  
    HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]][which(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]>=
                                                                               quantile(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]],alphaa))]
  lower_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- min(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  upper_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- max(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  
}


end_time <- Sys.time()   
end_time - start_time

HDR_interval_pts
HDR_interval_pts_to_keep

lower_interval_pts
upper_interval_pts
stored<-floor(unlist(stored)[!is.na(unlist(stored))])
stored<-pmax(lower_interval_pts, pmin(stored, upper_interval_pts))

x<-data.frame( "Lower HDR"=lower_interval_pts,"Upper HDR"=upper_interval_pts,"Predicted"=stored, "True value"=data$Leicester[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])
x$"inclusion"<- ifelse(x$True.value<x$Upper.HDR & x$True.value>x$Lower.HDR,1,0)
x


plot(x$Predicted,type = "p",ylim = c(min(lower_interval_pts,stored,data$Leicester[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)]),
                                     max(upper_interval_pts,stored,data$Leicester[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])),ylab="HDR Interval",xlab = "Time(biweekly)",main = "95% HDR Interval for Leicester Data"
)
lines(x$True.value,type = "p",col="darkred",pch=16)
# lines(lower_interval_pts,type = "p",pch=3)
lines(lower_interval_pts,type = "l",lty=4)
lines(upper_interval_pts,type = "l",lty=4)
legend(x = "topright",lty = c(NA,NA,4,4),pch=c(16,1,NA,NA),
       col= c("darkred","black","black","black"),
       legend=c( "observed","Predicted","HDR Interval"))



###############################################################################################################
#Birmingham


set.seed(25)
library(readxl)
library(glarma)
library(psych)
setwd("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket")
data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 measles.xlsx")
N<-276550
S1<-276550
ar_order<-1
birth_data<-read_excel("C:/Users/ASUS/Desktop/Current Modelling/data from niket and grenfell/data sent by niket/60 cities.xlsx")
xn<-as.numeric(birth_data[15, 3:  25])
divided_vector <- xn / 26
biweekly_birth_data<- rep(divided_vector, each = 26)


start_value_of_k<-50
end_value_of_k<-100


forecasting_loop_starts_from <-571
forecasting_loop_ends_in <- 595

probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)    #"probabilities" are the predicted values
predicted_probabilities<- replicate(end_value_of_k-start_value_of_k+1,rep(NA,1),simplify = F)   #"predicted_probabilities" are probabilities of "probabilities" 
predicted_probabilities_unlist<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)
stored<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1 ,rep(NA,1),simplify = F)  #stored will be unlisting "probabilities"

no_iterations<- 1000000  #for HDR intervals
HDR_interval_pts <- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_probabilities<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)
HDR_interval_pts_to_keep<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,rep(NA,no_iterations),simplify = F)

lower_interval_pts<-c()
upper_interval_pts<-c()

#glarmamod<- replicate(forecasting_loop_ends_in-forecasting_loop_starts_from+1,list(NA),simplify = F)
glarmamod<- replicate(end_value_of_k-start_value_of_k+1,list(NA),simplify = F)
alphaa<-.01

start_time <- Sys.time() 
for(training_upto in forecasting_loop_starts_from:forecasting_loop_ends_in){
  infected_first_category <- data$Birmingham[1:(training_upto+1)]
  infected_second_category<-data$Leicester[1:(training_upto+1)]
  infected_third_category<-data$Coventry[1:(training_upto+1)]
  
  nstep<-length(infected_first_category)
  
  cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I548
  cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
  cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
  
  infected_first_category_success<-infected_first_category[-c(1:ar_order)] #this will be used in dbinom er x.. binomial e I2 theke..
  
  cum_fromI2<-(cumsum(infected_first_category[-1]))
  infected_first_category_failure<-S1-cum_fromI2
  infected_first_category_failure<-infected_first_category_failure[-c(1:(ar_order-1))]
  
  data_sf<-cbind(infected_first_category_success[-c(length(infected_first_category_success))],
                 infected_first_category_failure[-c(length(infected_first_category_success))])
  data_sf<-matrix(c(infected_first_category_success[-c(length(infected_first_category_success))],
                    infected_first_category_failure[-c(length(infected_first_category_success))]),
                  nrow = length(infected_first_category_success[-c(length(infected_first_category_success))]))
  #two deleted because 547 no entry random dhukbe and 548 no entry ta predict hobe .. 
  
  biweekly_birth_data_in_inverselogit<- biweekly_birth_data[ar_order:(training_upto)]
  t<- (ar_order+1): (training_upto+1)
  
  seasonality1<-cos(2*pi*t/52)
  seasonality2<-cos(2*pi*t/26)
  baby_boom_effect1<- c(rep(0,104),rep(1,182-104),rep(0,nstep-182))
  baby_boom_effect1<-baby_boom_effect1[-length(baby_boom_effect1)]
  if (ar_order == 1) {
    baby_boom_effect <- baby_boom_effect1
  } else {
    baby_boom_effect <- baby_boom_effect1[-c(1:(ar_order-1))]
  }
  
  
  for(k in start_value_of_k:end_value_of_k){
    
    new_row <- c(k, S1-(sum(data_sf[,1])+k) )
    
    temp_data_sf <- rbind(data_sf, new_row)
    
    log_infected_first_category<-infected_first_category
    log_infected_first_category<-data$Birmingham[1:(training_upto+1)]
    log_infected_first_category<-log_infected_first_category[-c(length(log_infected_first_category))]
    log_infected_first_category<- c(log_infected_first_category,k)
    log_infected_first_category[log_infected_first_category == 0] <- 1
    log_infected_first_category<-log(log_infected_first_category)
    
    log_infected_second_category<-infected_second_category
    log_infected_second_category<-data$Leicester[1:(training_upto+1)]
    log_infected_second_category<-log_infected_second_category[-c(length(log_infected_second_category))]
    log_infected_second_category<- c(log_infected_second_category,k)
    log_infected_second_category[log_infected_second_category == 0] <- 1
    log_infected_second_category<-log(log_infected_second_category)
    log_infected_second_category<-1.809428e-05*log_infected_second_category
    
    log_infected_third_category<-infected_third_category
    log_infected_third_category<-data$Coventry[1:(training_upto+1)]
    log_infected_third_category<-log_infected_third_category[-c(length(log_infected_third_category))]
    log_infected_third_category<- c(log_infected_third_category,k)
    log_infected_third_category[log_infected_third_category == 0] <- 1
    log_infected_third_category<-log(log_infected_third_category)
    log_infected_third_category<-1.888391e-05*log_infected_third_category
    
    infected_first_category_in_inverselogit<- data.frame("lag1"=dplyr::lag(log_infected_first_category, 1),
                                                         
                                                         "lag1_2nd"=dplyr::lag(log_infected_second_category, 1),
                                                         "lag1_3rd"=dplyr::lag(log_infected_third_category, 1))
    
    infected_first_category_in_inverselogit<-na.omit(infected_first_category_in_inverselogit, cols = "lag1")
    
    data_covariates<-matrix(c(rep(1,length(infected_first_category_in_inverselogit$lag1)),infected_first_category_in_inverselogit$lag1,
                              
                              infected_first_category_in_inverselogit$lag1_2nd,
                              infected_first_category_in_inverselogit$lag1_3rd,
                              biweekly_birth_data_in_inverselogit,seasonality1,seasonality2,t,baby_boom_effect),
                            nrow = length(infected_first_category_in_inverselogit$lag1))
    colnames(data_covariates) <- c("Intercept","Lag1","Lag1_2nd","Lag1_3rd","biweeklybirth","Seasonalit1","Seasonalit2","time","baby_boom_effect")
    
    glarmamod[[k-(start_value_of_k-1)]] <- glarma(temp_data_sf, data_covariates,thetaLags = c(2,4), type = "Bin", method = "FS",
                                                  residuals = "Identity", maxit = 100, grad = 1e-6)
    probabilities[[k-(start_value_of_k-1)]]<- glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]
    predicted_probabilities[[k-(start_value_of_k-1)]]<- dbinom(x=floor(glarmamod[[k-(start_value_of_k-1)]]$fitted.values[length(glarmamod[[k-(start_value_of_k-1)]]$fitted.values)]),
                                                               size = S1-(sum(data_sf[,1])+k),
                                                               prob = logistic(as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+  
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                 
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[6])*seasonality1[length(seasonality1)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[7])*seasonality2[length(seasonality2)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[8])*t[length(t)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[10])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-1])+
                                                                                 as.numeric(glarmamod[[k-(start_value_of_k-1)]]$delta[11])*(glarmamod[[k-(start_value_of_k-1)]]$residuals[length(glarmamod[[k-(start_value_of_k-1)]]$residuals)-2])
                                                               )   )
  }
  
  predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]]<- unlist(predicted_probabilities)[!is.na(unlist(predicted_probabilities))]
  #stored[[training_upto-(forecasting_loop_starts_from-1)]]<-median(unlist(probabilities)[!is.na(unlist(probabilities))])
  stored[[training_upto-(forecasting_loop_starts_from-1)]]<-unlist(probabilities)[!is.na(unlist(probabilities))][which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])]
  
  max_position<- which.max(predicted_probabilities_unlist[[training_upto-(forecasting_loop_starts_from-1)]])
  # i <- (max_position - 1) %/% length(glarmamod[[1]]) + 1  # List index
  # j <- (max_position - 1) %% length(glarmamod[[1]]) + 1    # Sub-index within each list element
  glarmamod1<-glarmamod[[max_position]]   
  
  
  HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]]<- rbinom(n=no_iterations,size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                              logistic(as.numeric(glarmamod1$delta[1])+
                                                                                         as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                         
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                         as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                         
                                                                                         as.numeric(glarmamod1$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                         as.numeric(glarmamod1$delta[6])*seasonality1[length(seasonality1)]+
                                                                                         as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                         as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                         as.numeric(glarmamod1$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                         as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                         as.numeric(glarmamod1$delta[11])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                       
                                                                              ) )
  
  HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]<- dbinom(x=HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]],
                                                                                        size=floor(S1-(sum(data_sf[,1])+median(start_value_of_k:end_value_of_k))),
                                                                                        prob=logistic(as.numeric(glarmamod1$delta[1])+
                                                                                                        as.numeric(glarmamod1$delta[2])*infected_first_category_in_inverselogit$lag1[length(infected_first_category_in_inverselogit$lag1)]+
                                                                                                        
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[3])*infected_first_category_in_inverselogit$lag1_2nd[length(infected_first_category_in_inverselogit$lag1_2nd)]+
                                                                                                        as.numeric(glarmamod1$delta[4])*infected_first_category_in_inverselogit$lag1_3rd[length(infected_first_category_in_inverselogit$lag1_3rd)]+
                                                                                                        
                                                                                                        as.numeric(glarmamod1$delta[5])*biweekly_birth_data_in_inverselogit[length(biweekly_birth_data_in_inverselogit)]+
                                                                                                        as.numeric(glarmamod1$delta[6])*seasonality1[length(seasonality1)]+
                                                                                                        as.numeric(glarmamod1$delta[7])*seasonality2[length(seasonality2)]+
                                                                                                        as.numeric(glarmamod1$delta[8])*t[length(t)]+
                                                                                                        as.numeric(glarmamod1$delta[9])*baby_boom_effect[length(baby_boom_effect)]+
                                                                                                        as.numeric(glarmamod1$delta[10])*(glarmamod1$residuals[length(glarmamod1$residuals)-1])+
                                                                                                        as.numeric(glarmamod1$delta[11])*(glarmamod1$residuals[length(glarmamod1$residuals)-2])   
                                                                                                      
                                                                                        ) )
  
  
  HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]]<-  
    HDR_interval_pts[[training_upto-(forecasting_loop_starts_from-1)]][which(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]]>=
                                                                               quantile(HDR_interval_probabilities[[training_upto-(forecasting_loop_starts_from-1)]],alphaa))]
  lower_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- min(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  upper_interval_pts[training_upto-(forecasting_loop_starts_from-1)] <- max(HDR_interval_pts_to_keep[[training_upto-(forecasting_loop_starts_from-1)]])
  
}


end_time <- Sys.time()   
end_time - start_time

HDR_interval_pts
HDR_interval_pts_to_keep

lower_interval_pts
upper_interval_pts
stored<-floor(unlist(stored)[!is.na(unlist(stored))])
stored<-pmax(lower_interval_pts, pmin(stored, upper_interval_pts))

x<-data.frame( "Lower HDR"=lower_interval_pts,"Upper HDR"=upper_interval_pts,"Predicted"=stored, "True value"=data$Birmingham[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])
x$"inclusion"<- ifelse(x$True.value<x$Upper.HDR & x$True.value>x$Lower.HDR,1,0)
x

plot(x$Predicted,type = "p",ylim = c(min(lower_interval_pts,stored,data$Birmingham[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)]),
                                     max(upper_interval_pts,stored,data$Birmingham[(forecasting_loop_starts_from+1):(forecasting_loop_ends_in+1)])),ylab="HDR Interval",xlab = "Time(biweekly)",main = "95% HDR Interval for Birmingham Data"
)
lines(x$True.value,type = "p",col="darkred",pch=16)
# lines(lower_interval_pts,type = "p",pch=3)
lines(lower_interval_pts,type = "l",lty=4)
lines(upper_interval_pts,type = "l",lty=4)
# legend(x = "topright",lty = c(NA,NA,4,4),pch=c(16,1,NA,NA),
#        col= c("darkred","black","black","black"),
#        legend=c( "observed","Predicted","HDR Interval"),cex=0.2)






