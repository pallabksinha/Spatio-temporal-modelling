e2<-infected_first_category_success-predicted

eta<-c()
wr<-c()

# for (t in (ar_order+1)  :training_upto){
#   eta[t-ar_order]<-exp(coeff[t-ar_order]) / ((1+exp(coeff[t-ar_order]))^2)
#   wr[t-ar_order]<- (e2[t-ar_order]) / eta[t-ar_order]
# }
# 
# wr
# length(wr)
# mean_wr<-mean((wr[is.finite(wr)]))
# mean_wr

MSE_wr <- mean((model$residuals)^2)
##############################################################################################################################################

# cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
# length(cum)
# cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
# length(cum)
# cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
# length(cum)
# cum<- cum[-c(1:(ar_order-1))]
# length(cum)
n<-n_trials
#aic
x<-infected_first_category[-c(1:ar_order)] # this is I7,I8,...,I546
length(x)

no_of_estimated_parameter<-length(coefficients(model))

probability<-predict_prob

length(probability)


aic<-  model$aic
bic<-  model$aic-2*no_of_estimated_parameter + (no_of_estimated_parameter*log(training_upto))
aic
bic

#pearsonian residual
#pearsonian residual
probability_s<-  predict_prob
length(probability_s)
probability_f<- 1-probability_s
standard_deviation<- sqrt(n*probability_s*probability_f) 
pearsonian_resid <- (e2) / standard_deviation

MSE_pearson<- mean(pearsonian_resid^2)     
# this is same as :  sd(residuals(model,type = "pearson"))




#Chi_squared statistic
variance<- n*probability_s*probability_f
chi_squared_statistic<-mean(   ((e2)^2) / predicted   )
chi_squared_statistic


#deviance statistic
probability_for_deviance<- x/n
length(probability_for_deviance)
D<- sqrt(2*( sum(dbinom(x,n,probability_for_deviance,log=T))  -   sum(dbinom(x,n,probability,log=T)) ))
D

#this is same as sqrt(model$deviance).. 
#in page 26, eq 1.72 mentions to add a sign(e_i).. that has not been added because then we cannot compare


#RMSE
rmse<- sqrt(mean(e2^2))
rmse
#mae
mae<- mean(abs(e2))


gof<-data.frame("Model"=2,"df"=training_upto-no_of_estimated_parameter,"Pearson Residual" =MSE_pearson,
                "WR"=MSE_wr,"AIC"=aic,"D"=D,"box_lunj_pvalue_on raw resid"=Box_ljung_test$p.value,
                "box_lunj_pvalue_on rqr"=Box_ljung_test_rq$p.value,"ks p value"=ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )$p.value) 

gof
############################

library(openxlsx)
library(readxl)
setwd("C:/Users/ASUS/Desktop/Bath/paper purpose/Coventry spatial connection/prediction/set 6(1)/Coventry/set 8")

# Read the existing Excel file
existing_data <- read_excel("GOF statistics.xlsx")

# Create a new row as a data.frame
new_row <- data.frame("Model"=2,"df"=training_upto-no_of_estimated_parameter,"Pearson Residual" =MSE_pearson,
                      "WR"=MSE_wr,"AIC"=aic,"D"=D,"box_lunj_pvalue_on raw resid"=Box_ljung_test$p.value,
                      "box_lunj_pvalue_on rqr"=Box_ljung_test_rq$p.value,"ks p value"=ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )$p.value)  
# Add the new row to the existing data
updated_data <- rbind(existing_data, new_row)
# Write the updated data back to the Excel file
write.xlsx(updated_data, file = "GOF statistics.xlsx")







