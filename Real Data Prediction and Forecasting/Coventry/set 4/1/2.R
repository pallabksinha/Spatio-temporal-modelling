data<-data[,-c(1,2)]
data<-cbind(rep(1,length(data$Seasonalit2)),data)
m<-as.vector(model$coefficients)
coeff <- numeric(nrow(data))
# Loop through each row of the data matrix
for (i in 1:nrow(data)) {
  coeff[i] <- sum(data[i, ] * m)
}
predict_prob<- model$fitted.values
length(predict_prob)
head(coeff,10) 
head(model$linear.predictors,10)


#############################################3
n_trials <- rowSums(data_sf)
y <- infected_first_category_success
length(y)


# cum<-cumsum(infected_first_category[-1]) #this is I2, I2+I3, I2+I3+I4, ... ,I2+...+I546
# length(cum)
# cum<-c(S1,S1-cum) #This is S1, S1-I2, S1-(I2+I3), ... ,S1 - (I2+...+I546)
# length(cum)
# cum<-cum[-nstep]   #S_546 is deleted because upto S_545 is required
# length(cum)
# cum<-cum[-c(1:(ar_order-1))] 
# length(cum)
# n<- cum
# length(n)


residuals_rq <- numeric(length(y))
for (i in 1:length(y)) {
  
  # CDF for y - 1 (lower bound)
  p_lower <- pbinom(y[i] - 1, size = n_trials[i], prob = predict_prob[i])
  
  # CDF for y (upper bound)
  p_upper <- pbinom(y[i], size = n_trials[i], prob = predict_prob[i])
  
  # Generate a random uniform value between the lower and upper CDF
  u <- runif(1, p_lower, p_upper)
  
  # Transform to the normal quantile
  residuals_rq[i] <- qnorm(u)
}

###############################################################################
#Residula_rq is being used as it is

residuals_rq<-residuals_rq[is.finite(residuals_rq)]
plot(residuals_rq,ylim = c(-7,7))
plot(residuals_rq)


qqnorm(residuals_rq, main = "QQ Plot for Normality")
qqline(residuals_rq, col = "red", lwd = 2)  # Add a reference line


qqnorm(residuals_rq, main = "QQ Plot for Normality")
abline(0, 1, col = "red", lwd = 2) 


acf(residuals_rq ,xlab="Lag",main="Acf plot for errors",ci=0.90)
pacf(residuals_rq,   xlab="Lag",main="PAcf plot for errors",ci=0.90)

ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )
ks.test(residuals_rq, "pnorm",mean = 0, sd = 1 )



######################################################################################
#Residula_rq is being used after being normalized


residuals_rq<-residuals_rq[is.finite(residuals_rq)]
residuals_rq<- scale(residuals_rq)

plot(residuals_rq)
plot(residuals_rq,ylim = c(-3,3),ylab = "Residuals")

qqnorm(residuals_rq, main = "QQ Plot for Normality")
qqline(residuals_rq, col = "red", lwd = 2)  # Add a reference line


qqnorm(residuals_rq, main = "QQ Plot for Normality")
abline(0, 1, col = "red", lwd = 0.2) 

acf(residuals_rq ,xlab="Lag",main="Acf plot for errors",ci=0.99)
pacf(residuals_rq,   xlab="Lag",main="PAcf plot for errors",ci=0.99)

ks.test(residuals_rq, "pnorm",mean = mean(residuals_rq), sd = sd(residuals_rq) )
ks.test(residuals_rq, "pnorm",mean = 0, sd = 1 )

Box_ljung_test_rq<- Box.test(residuals_rq,lag=500,type = "Ljung-Box")
Box_ljung_test_rq

############################################################################################
############################################################################################
############################################################################################

par(mfrow = c(2, 2))
acf(model$residuals   ,xlab="Lag",main="Acf plot for errors",ci=0.99)
pacf(model$residuals,   xlab="Lag",main="PAcf plot for errors",ci=0.99)
plot(residuals_rq,ylim = c(-3,3),ylab = "Residuals")
qqnorm(residuals_rq, main = "QQ Plot for Normality")
abline(0, 1, col = "red", lwd = 0.2) 
par(mfrow = c(1, 1))











