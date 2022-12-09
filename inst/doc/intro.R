## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#Sample1
set.seed(114514)
matr = matrix(nrow = 15, ncol = 15)
for(i in 1:15)
  for(j in 1:15)
  {
    matr[i, j] = 0.5 ^ abs(i-j)
  }
mu = rep(0, 15)
Z = MASS::mvrnorm(n = 50, mu, matr)
X = matrix(0, nrow = 50, ncol = 30)
for(i in 1:50)
  for(j in 1:15)
  {
    if(Z[i, j] < qnorm(1/3))
    {
      X[i, 2*j-1] = 1
    }
    else if(Z[i, j] < qnorm(2/3))
    {
      X[i, 2*j] = 1
    }
  }
coef = rep(0, 30)
coef[2] = 1.8
coef[1] = -1.2
coef[6] = 1
coef[5] = 0.5
coef[10] = 1
coef[9] = 1

Y = X %*% coef

sigma = sqrt(var(Y) * 49 / 10 / 50)

Y = Y + rnorm(50, 0, sigma)

qr = qr(scale(X, scale = FALSE))
X = qr.Q(qr)
Y = scale(Y, scale = FALSE)

p <- rep(2,15)
J = 15

## -----------------------------------------------------------------------------
norm_1 <- function(x){
  return (sum(abs(x)))
}
norm_2 <- function(x){
  return (sqrt(t(x) %*% x))
}
norm_K <- function(x,p){
  return (sqrt(p * t(x) %*% x))
}
pos <- function(x){
  return (max(x,0))
}

g_lasso <- function(X, Y, lambda, p, J, index){#index是p的前缀和
  beta <- rep(0, sum(p))
  
  while(1){
    beta1 <- beta
    
  X1 <- X[,1:index[1]]
  S1 <- t(X1) %*% (Y - X %*% c(rep(0,p[1]),beta[(index[1]+1):index[J]]) )
  beta[1:p[1]] <- pos(1 - lambda * sqrt(p[1]) / norm_2(S1)) * S1
  #print( S1)
  
  for (j in 2:J){
    Xj <- X[, (index[j-1]+1):index[j]]
    
    if (j<J) Sj <- t(Xj) %*% (Y - X %*% c(beta[1:index[j-1]],rep(0,p[j]),beta[(index[j]+1):index[J]]) )
    else Sj <- t(Xj) %*% (Y - X %*% c(beta[1:index[j-1]],rep(0,p[j]) ))
    
    beta[(index[j-1]+1) : index[j]] <- pos(1 - lambda * sqrt(p[j]) / norm_2(Sj)) * Sj
  }
  #print(beta)
  eps <- norm_1(beta - beta1)
  if (abs(eps) < 1e-6) break
  }
  return (beta)
}

## -----------------------------------------------------------------------------
index <- numeric(length(p))
for (i in 1:J)
  index[i] <- sum(p[1:i])

X1 <- X[, 1:index[1]]
lambda_upp <- norm_2(t(X1) %*% Y)/sqrt(p[1])
for (j in 2:J){
  Xj <- X[, (index[j-1]+1):index[j]]
  lambda_upp = max(lambda_upp, norm_2(t(Xj) %*% Y)/sqrt(p[1]))
}

M <- 50
lam <- seq(0,lambda_upp,lambda_upp / M) #参数取值范围 10为参数取值数量
beta_s <- matrix(sum(p)*(M+1), sum(p), M+1) #p为pj向量
beta_s[,1] <- rep(0,sum(p)) #beta存储阵 第一列为0

for (j in 1:M){
  betaj <- g_lasso(X,Y,lam[j],p,J, index)
  beta_s[,j+1] <- betaj
}


## -----------------------------------------------------------------------------
beta_s

## -----------------------------------------------------------------------------
load(file="../data/test_public.rda")
load(file="../data/hourly.rda")

## -----------------------------------------------------------------------------
library(tseries)
library(forecast)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
j1 <- which(hourly$train.or.test=="test1")[1]
j2 <- which(hourly$train.or.test=="test2")[1]
j3 <- which(hourly$train.or.test=="test3")[1]
j4 <- which(hourly$train.or.test=="test4")[1]
j10 <- which(hourly$train.or.test=="test1")[168]
j20 <- which(hourly$train.or.test=="test2")[168]
j30 <- which(hourly$train.or.test=="test3")[168]
j40 <- which(hourly$train.or.test=="test4")[168]

c(j1,j2,j3,j4)
c(j10,j20,j30,j40)
c(nrow(train1),nrow(train2),nrow(train3),nrow(train4))
c(nrow(test1),nrow(test2),nrow(test3),nrow(test4))
mx <- numeric(20)
for (k in 2:21)
  mx[k-1] <- mean(hourly[,k],na.rm = T)

## -----------------------------------------------------------------------------
predict <- function(i){
  if (is.na(train1[1,i])==TRUE){
    train1[1,i] = 0
  }
  for( j in 2:nrow(train1)){
    if (is.na(train1[j,i])==TRUE & train1$train.or.test[j]=="train"){
      if(j<=72) train1[j,i] = train1[j-1,i ] #用前一天代替
      else train1[j,i] = mean(train1[c(j-72,j-48,j-24),i])
    }
  }
  train1[which(train1[,i]>=5*mx[i-1]),i] <- mx[i-1] #去掉异常大值

  if (is.na(train2[1,i])==TRUE){
    train2[1,i] = 0
  }
  for( j in 2:nrow(train2)){
    if (is.na(train2[j,i])==TRUE & train2$train.or.test[j]=="train"){
      if(j<=72) train2[j,i] = train2[j-1,i ] #用前一天代替
      else train2[j,i] = mean(train2[c(j-72,j-48,j-24),i])
    }
  }
  train2[which(train2[,i]>=5*mx[i-1]),i] <- mx[i-1] #去掉异常大值


  if (is.na(train3[1,i])==TRUE){
    train3[1,i] = 0
  }
  for( j in 2:nrow(train3)){
    if (is.na(train3[j,i])==TRUE & train3$train.or.test[j]=="train"){
      if(j<=72) train3[j,i] = train3[j-1,i ] #用前一天代替
      else train3[j,i] = mean(train3[c(j-72,j-48,j-24),i])
    }
  }
  train3[which(train3[,i]>=5*mx[i-1]),i] <- mx[i-1] #去掉异常大值


  if (is.na(train4[1,i])==TRUE){
    train4[1,i] = 0
  }
  for( j in 2:nrow(train4)){
    if (is.na(train4[j,i])==TRUE & train4$train.or.test[j]=="train"){
      if(j<=72) train4[j,i] = train4[j-1,i ] #用前一天代替
      else train4[j,i] = mean(train4[c(j-72,j-48,j-24),i])
    }
  }
  train4[which(train4[,i]>=5*mx[i-1]),i] <- mx[i-1] #去掉异常大值
  
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

## -----------------------------------------------------------------------------
predict(3)

## -----------------------------------------------------------------------------
dates(2022,12)

