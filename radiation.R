library(dplyr)
library(zoo)
library(forecast)
library(tseries)
library(astsa)
library(MASS)
library(gridExtra)
library(ggplot2)

############################# formatting ###########################
radiation  = read.csv("C:\\Users\\lauyi\\Desktop\\math4425\\project\\dataset\\monthly_KP_GSR_ALL.csv", header =TRUE)
radiation$value = as.numeric(radiation$value)
radiation

########################### split testing ##############################
gsr = radiation$value
time = radiation$yyyy_mm

n_gsr = length(gsr)
n_test = round(n_gsr*0.2)
n_train = n_gsr - n_test

d_train = radiation$value[1:n_train]
d_test = radiation$value[(n_train+1):n_gsr]
t_train = radiation$yyyy_mm[1:n_train]
t_test = radiation$yyyy_mm[(n_train+1):n_gsr]


########################## data visualization ##########################
is.not.null <- function(x) !is.null(x)
t_plot = function(values, time, n_year = NULL, tt = "", ylab = ""){
  
  if (is.not.null(n_year)){
    values = values[1:(12*n_year)]
    #time = time[1:(12*n_year)]
  }
  
  n = length(values)
  time = time[1:n]
  
  #### plot the values ####
  #dev.new(width=15, height=4, unit="in")
  plot(values, type = "o", pch=20, col = "blue", xaxt ="n", xlab="",  ylab = ylab, main = tt)
  axis(1, at = 1:n, labels=time, las = 2)
  title(xlab = "year", line = 5)
  
  #### add red lines for marking the end of the year ####
  for (i in 1:n){
    if (i%%12 ==1 ){
      abline(v = i, col = "red", lty =2, lwd = 1)
    }
  }
}

t_plot(d_train, t_train, ylab = "MGSR (MJ/m^2)", tt = "MGSR data from 1992-07 to 2022-03")
t_plot(d_ln_train, t_train, ylab = "ln(MGSR)", tt = "ln(MGSR) data from 1992-07 to 2022-03")

acf(d_ln_train, 50)
pacf(d_ln_train, 50)



t_plot(diff(d_ln_train,12), t_train, ylab = "W_t", tt = "Differenced MGSR data from 1992-07 to 2022-03")

########################### log transformation ##############################
d_ln_train = log(d_train)
d_index = c(1:n_train)
d_ln_test = log(d_test)


########################### difference transformation ##############################

w0 = diff(d_ln_train,12)
acf(w0, 50)
pacf(w0, 50)

t_plot(w0, t_train)

### conclude W_t = (1 - B^12)log(Z_t) 

########################### find p, q, P, Q + model checking ##############################
## ARIMA


## L-JUNG BOX TEST
# refernce: number of lags
# https://robjhyndman.com/hyndsight/ljung-box-test/


#H0: The data are independently distributed (i.e. the correlations in the population from which the sample is taken are 0, so that any observed correlations in the data result from randomness of the sampling process).
#Ha: The data are not independently distributed; they exhibit serial correlation.


#------------------------------- fit arima function ---------------------------------------#
fit_arima = function(data, p,d,q, P = 0 ,D = 0 ,Q = 0 ,S = 12){
  
  ### seasonal differencing ####
  if (D>=1){
      for (i in 1: D){
          data = diff(data , S)
      }
  }
  
  ### arima differencing ####
  if (d>=1){
    for (i in 1: d){
      data = diff(data)
    }
  }
  
  #### p, q, P, Q first fitting ####
  #est = sarima(data, p,d,q,P,D,Q,S, details = F)
  pdq = c(p,0,q); PDQ = c(P,0,Q)
  est=arima(data, pdq, seasonal = list(order = PDQ,  period = S))
  n_coef = length(est$coef); s_error = sqrt(diag(est$var.coef))
  
  #### p, q, P, Q second fitting ####
  set_zero = abs(est$coef) < 2*s_error
  fixed = rep(NA,n_coef); fixed[which(set_zero ==TRUE)] = 0
  #est = sarima(data, p,d,q,P,D,Q,S, details = F, fixed = fixed)
  #print(est)
  est = arima(data, pdq, seasonal = list(order = PDQ, period = S), fixed = fixed, transform.pars = FALSE)
  
  #### L-JUNG BOX TEST ####
  box = Box.test(est$residuals,lag = 24,type="Ljung")
  pv=1-pchisq(box$statistic, box$parameter-n_coef)
  
  return(c(p,d,q,P,D,Q, S, pv, pv<0.05, est$aic))
  
}

get_coef_arima = function(data, p,d,q,P,D,Q, S=12){
  est = sarima(data, p,d,q, P,D,Q, S, details = FALSE)
  n_coef = length(est$fit$coef); s_error = sqrt(diag(est$fit$var.coef))
  
  set_zero = abs(est$fit$coef) < 2*s_error
  fixed = rep(NA,n_coef); fixed[which(set_zero ==TRUE)] = 0
  
  est = sarima(data,p,d,q, P,D,Q,S, details = FALSE, fixed = fixed)
  return(est)
}

get_fixed_arima = function(data, p,d,q,P,D,Q, S=12){
  est = sarima(data, p,d,q, P,D,Q, S, details = FALSE)
  n_coef = length(est$fit$coef); s_error = sqrt(diag(est$fit$var.coef))
  
  set_zero = abs(est$fit$coef) < 2*s_error
  fixed = rep(NA,n_coef); fixed[which(set_zero ==TRUE)] = 0

  return(fixed)
}

plot_arima = function(data, p,d,q,P,D,Q, S=12){
  est = sarima(data, p,d,q, P,D,Q, S, details = FALSE)
  n_coef = length(est$fit$coef); s_error = sqrt(diag(est$fit$var.coef))
  
  set_zero = abs(est$fit$coef) < 2*s_error
  fixed = rep(NA,n_coef); fixed[which(set_zero ==TRUE)] = 0
  
  est = sarima(data,p,d,q, P,D,Q,S, fixed = fixed)
  #return(est)
}

get_pred = function(data, p,d,q,P,D,Q, S=12){
  fixed = get_fixed_arima(data, p,d,q, P,D,Q,S)
  forecast = sarima.for(data,n_test, p,d,q, P,D,Q,S, fixed = fixed)
  return(forecast$pred)
}

#--------------------------------iterative fit------------------------------------------------#

iter_fit = function(data, vp, vd, vq, vP, vD, vQ, S=12){
  
  df = NULL
  
  for (p in vp){
    for (d in vd){
      for (q in vq){
        for (P in vP){
          for (D in vD){
            for (Q in vQ){
              
              tryCatch(
                expr = {
                  df = rbind(df, fit_arima(data, p,d,q,P,D,Q,S))
                },
                
                error = function(e){
                  print(paste0("error in (",p,",",d,",",q,")(",P,",",D,",",Q,")_12"))
                },
                
                warning = function(w){
                  print(paste0("NaN produced in (",p,",",d,",",q,")(",P,",",D,",",Q,")_12"))
                }
                
              )
              
            }
          }
        }
      }
    }
    
  }
  
  colnames(df) = c("p", "d", "q", "P", "D","Q", "S", "L-jung test p-value", "H1", "AIC")
  df = as.data.frame(df)
  df = df[order(df$AIC),]
  
  
  return(df)  
}


#------------------------------ model selection (main of step 3-4) -------------------------------------------#
result = iter_fit(d_ln_train, vp=c(0:1),vd=0,vq=c(0:1),vP=c(0,1),vD=1,vQ=c(0:2),S = 12)
result
rownames(result) = NULL
grid.arrange(tableGrob(result))

est = get_coef_arima(d_ln_train, 0,0,1,1,1,1)
plot_arima(d_ln_train, 0,0,1,1,1,1)
grid.arrange(tableGrob(est$ttable))

###################################### prediction #########################################################


#-------------------------------------- plot log prediction --------------------------------------------#
p_plot = function(d_train, d_test, n_train, n_test, p,d,q, P,D,Q,S, tt = "", ylab = ""){
    
    #combine train and test data
    split = length(d_train)
    ind = c((split - n_train+1): (split+n_test))
    real = d_train[(split-n_train+1): split]
    real = c(real, d_test[1:n_test])
    
    #forecasting
    ind_test = c((split+1):(split+n_test))
    fixed = get_fixed_arima(d_train, p,d,q, P,D,Q,S)
    forecast = sarima.for(d_train,n_test, p,d,q, P,D,Q,S, fixed = fixed)
    
    #get predictions and standard error and U / L
    pred = forecast$pred
    U = pred + 1.96*forecast$se
    L = pred - 1.96*forecast$se
    
    #plotting
    plot(ind, real, type="o", xlab="", xaxt = "n", ylab = ylab,ylim = c(min(L)-0.1,max(U)+0.1) ,lwd = 2,main=tt)
    polygon(c(ind_test, rev(ind_test)), c(L, rev(U)), col = "pink", density = 10)
    lines(ind_test,d_ln_test, type="o")
    lines(ind_test,pred,type="o",col="red", lty = 1, lwd =2, pch = 1)
    
    axis(1, at = ind, labels=time[ind], las = 2)
    
    legend(x="topleft",c("True Value","Prediction","Forecast Interval"),lty=c(1,2,1),pch=c(1,1),col=c("black","red", "pink"))
    
}


#------------------------------------ plot original values ---------------------------------------------------#

o_plot = function(d_train, d_test, n_train, n_test, p,d,q, P,D,Q,S, tt="", ylab = ""){
  
  #combine train and test data
  split = length(d_train)
  ind = c((split - n_train+1): (split+n_test))
  real = d_train[(split-n_train+1): split]
  real = c(real, d_test[1:n_test])
  
  #forecasting
  ind_test = c((split+1):(split+n_test))
  fixed = get_fixed_arima(d_train, p,d,q, P,D,Q,S)
  forecast = sarima.for(d_train,n_test, p,d,q, P,D,Q,S, fixed = fixed)
  
  #get predictions and standard error and U / L
  pred = forecast$pred
  U = pred + 1.96*forecast$se; U = exp(U)
  L = pred - 1.96*forecast$se; L = exp(L)
  pred = exp(pred)
  
  #plotting
  plot(ind, exp(real), type="o", xlab="", xaxt = "n", ylab = ylab, ylim = c(min(L)-0.1,max(U)+0.1) ,lwd = 2,main=tt)
  polygon(c(ind_test, rev(ind_test)), c(L, rev(U)), col = "pink", density = 10)
  lines(ind_test,exp(d_ln_test), type="o")
  lines(ind_test,pred,type="o",col="red", lty = 1, lwd =2, pch = 1)
  
  axis(1, at = ind, labels=time[ind], las = 2)
  
  legend(x="topleft",c("True Value","Prediction","Forecast Interval"),lty=c(1,2,1),pch=c(1,1),col=c("black","red", "pink"))
  
}

#-------------------------------------------main of prediction ---------------------------------##

### with forecast interval ###
p_plot(d_ln_train, d_ln_test, 50,n_test,0,0,1,1,1,1,12, ylab = "Z_t", tt="Forecasting with Z_t = ln(y_t) using Model 1 SARIMA(0,0,1) X (1,1,1)_12")
o_plot(d_ln_train, d_ln_test, 50,n_test,0,0,1,1,1,1,12, ylab = "y_t", tt="Forecasting with Original Data y_t using Model 1 SARIMA(0,0,1) X (1,1,1)_12")

### MSE ###
pred = get_pred(d_ln_train, 0,0,1,1,1,1)
#log version
sum((d_ln_test - pred)**2)/(n_test)
#original version
sum((exp(d_ln_test) - exp(pred))**2)/(n_test)



