#' @title Using arima model to predict water supply time series in Shenzhen
#' @name arima
#' @description Using arima model to predict water supply time series in Shenzhen

#' @param i a number between 2 and 21, means the selected column
#'
#' @return full dataset : hourly
#' @import tseries
#' @import forecast
#' @examples
#' \dontrun{
#' predict(2)
#' }
#' @export
predict <- function(i){
  load(file="../data/test_public.rda")
  load(file="../data/hourly.rda")
  train1 = hourly[(hourly['time']>='2022-01-01 01:00:00')&(hourly['time']<'2022-05-01 01:00:00'),]
  test1 = hourly[(hourly['time']>='2022-05-01 01:00:00')&(hourly['time']<'2022-05-08 01:00:00'),]
  
  train2 = hourly[(hourly['time']>='2022-05-08 01:00:00')&(hourly['time']<'2022-06-01 01:00:00'),]
  test2 = hourly[(hourly['time']>='2022-06-01 01:00:00')&(hourly['time']<'2022-06-08 01:00:00'),]
  
  train3 = hourly[(hourly['time']>='2022-06-08 01:00:00')&(hourly['time']<'2022-07-21 01:00:00'),]
  test3 = hourly[(hourly['time']>='2022-07-21 01:00:00')&(hourly['time']<'2022-07-28 01:00:00'),]
  
  train4 = hourly[(hourly['time']>='2022-07-28 01:00:00')&(hourly['time']<'2022-08-21 01:00:00'),]
  test4 = hourly[(hourly['time']>='2022-08-21 01:00:00')&(hourly['time']<'2022-08-28 01:00:00'),]
  j1 <- which(hourly$train.or.test=="test1")[1]
  j2 <- which(hourly$train.or.test=="test2")[1]
  j3 <- which(hourly$train.or.test=="test3")[1]
  j4 <- which(hourly$train.or.test=="test4")[1]
  j10 <- which(hourly$train.or.test=="test1")[168]
  j20 <- which(hourly$train.or.test=="test2")[168]
  j30 <- which(hourly$train.or.test=="test3")[168]
  j40 <- which(hourly$train.or.test=="test4")[168]
  mx <- numeric(20)
  for (k in 2:21)
    mx[k-1] <- mean(hourly[,k],na.rm = T)
  
  if (is.na(train1[1,i])==TRUE){
    train1[1,i] = 0
  }
  for( j in 2:nrow(train1)){
    if (is.na(train1[j,i])==TRUE & train1$train.or.test[j]=="train"){
      if(j<=72) train1[j,i] = train1[j-1,i] 
      else train1[j,i] = mean(train1[c(j-72,j-48,j-24),i])
    }
  }
  train1[which(train1[,i]>=5*mx[i-1]),i] <- mx[i-1] 
  
  if (is.na(train2[1,i])==TRUE){
    train2[1,i] = 0
  }
  for( j in 2:nrow(train2)){
    if (is.na(train2[j,i])==TRUE & train2$train.or.test[j]=="train"){
      if(j<=72) train2[j,i] = train2[j-1,i ] 
      else train2[j,i] = mean(train2[c(j-72,j-48,j-24),i])
    }
  }
  train2[which(train2[,i]>=5*mx[i-1]),i] <- mx[i-1] 
  
  
  if (is.na(train3[1,i])==TRUE){
    train3[1,i] = 0
  }
  for( j in 2:nrow(train3)){
    if (is.na(train3[j,i])==TRUE & train3$train.or.test[j]=="train"){
      if(j<=72) train3[j,i] = train3[j-1,i ] 
      else train3[j,i] = mean(train3[c(j-72,j-48,j-24),i])
    }
  }
  train3[which(train3[,i]>=5*mx[i-1]),i] <- mx[i-1] 
  
  
  if (is.na(train4[1,i])==TRUE){
    train4[1,i] = 0
  }
  for( j in 2:nrow(train4)){
    if (is.na(train4[j,i])==TRUE & train4$train.or.test[j]=="train"){
      if(j<=72) train4[j,i] = train4[j-1,i ] 
      else train4[j,i] = mean(train4[c(j-72,j-48,j-24),i])
    }
  }
  train4[which(train4[,i]>=5*mx[i-1]),i] <- mx[i-1] 
  
  hourly[1:nrow(train1),] <- train1
  hourly[(j10+1):(j10+nrow(train2)),] <- train2
  hourly[(j20+1):(j20+nrow(train3)),] <- train3
  hourly[(j30+1):(j30+nrow(train4)),] <- train4
  
  for(j in j1:(j1+167)){
    flow_train <- hourly[(j-24*(7:1)),i]
    fit <- auto.arima(flow_train)
    flow_test <- forecast(fit,h=1)$mean
    if (flow_test<0.001) flow_test <- mean(hourly[c(j-(26:24)),i])
    hourly[j,i] <- flow_test
  }
  for(j in j2:(j2+167)){
    flow_train <- hourly[(j-24*(7:1)),i]
    fit <- auto.arima(flow_train)
    flow_test <- forecast(fit,h=1)$mean
    if (flow_test<0.001) flow_test <- mean(hourly[c(j-(26:24)),i])
    hourly[j,i] <- flow_test
  }
  for(j in j3:(j3+167)){
    flow_train <- hourly[(j-24*(4:1)),i]
    fit <- auto.arima(flow_train)
    flow_test <- forecast(fit,h=1)$mean
    if (flow_test<0.001) flow_test <- mean(hourly[c(j-(26:24)),i])
    hourly[j,i] <- flow_test
  }
  for(j in j4:(j4+167)){
    flow_train <- hourly[(j-24*(4:1)),i]
    fit <- auto.arima(flow_train)
    flow_test <- forecast(fit,h=1)$mean
    if (flow_test<0.001) flow_test <- mean(hourly[c(j-(26:24)),i])
    hourly[j,i] <- flow_test
  }
  return (hourly)
}
