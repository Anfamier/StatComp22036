## -----------------------------------------------------------------------------
library(knitr)
library(ggplot2)
library(MASS)

## -----------------------------------------------------------------------------
summary(Cars93)

## -----------------------------------------------------------------------------
xtable::xtable(head(Cars93))

## -----------------------------------------------------------------------------
attach(Cars93)
plot(x=Horsepower, y=MPG.highway, main = "Horsepower-MPG.highway")

## -----------------------------------------------------------------------------
hist(Cars93$Weight)

## -----------------------------------------------------------------------------
pie(table(Cars93$Type),main = "Car types in Car93")

## -----------------------------------------------------------------------------
n <- 1000
set.seed(0)
u <- runif(n)
x <- 2/ sqrt(1-u) # F(x) = 1- (2/x)^2, 2<=x
hist(x, prob = TRUE, main = expression(f(x)==frac(8,x^3)), breaks = 100, xlim=c(0,50))
y <- seq(2, 50, .01)
lines(y, 8/(y^3))


## -----------------------------------------------------------------------------
n <- 1e3;j<-k<-0;y <- numeric(n)
ar_beta <- function(a,b){
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1) #random variate from g().
    if (x^(a-1) * (1-x)^(b-1) > u) {
      #we accept x
      k <- k + 1
      y[k] <- x
    }
  }
  return (list(j,y))
}

## -----------------------------------------------------------------------------
ans <- ar_beta(3,2)
ans[[1]]/n# value of c

x <- seq(0,1,0.01)

hist(ans[[2]], prob=T,main = expression(Be(3,2)), xlim=c(0,1), xlab="x")

lines(x, 12*x^2*(1-x))

## -----------------------------------------------------------------------------
n <- 1000; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda) # the length of lambda = n
hist(x, prob = TRUE, breaks = 20, main = "Mixture Distribution", ylim = c(0,2))

## -----------------------------------------------------------------------------
hist(x, prob = TRUE, breaks = 20, main = expression(f(x)==frac(64,(2+x)^5)), ylim = c(0,2))
y <- seq(0, 50, .01)
lines(y, 64/(2+y)^5 )


## -----------------------------------------------------------------------------
quick_sort <- function(x){
  num<-length(x)
  if (num==0||num==1){
    return(x)
  }
  else{
    a<-x[1]
    y<-x[-1]
    lower <- y[y<a]
    upper <-y [y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))
  }
}

## -----------------------------------------------------------------------------
e <- c(10^4, 2 * 10^4, 4 * 10^4, 6 * 10^4, 8 * 10^4)
a <- numeric(5)
t <- e*log(e)
for (i in 1:5){
  sum_t <- 0
  j <- 1
  for (j in 1:100){
    test <- sample(e[i])
    sum_t <- sum_t + system.time(quick_sort(test))[1]
    #system.time(quick_sort(test))[1]
  }
  a[i] <- sum_t/100
}

## -----------------------------------------------------------------------------
fit <- lm(a~t)
summary(fit)
fit$coef
plot(x=t,y=a,main = "time - nlogn")
y <- seq(0,10^6,1)
lines(y,fit$coef[1]+fit$coef[2]*y,col="red")

## -----------------------------------------------------------------------------
MC.Phi <- function(R = 10000, antithetic = FALSE) {
  u <- runif(R/2)
  if (antithetic) v <- 1 - u else v <- runif(R/2)
  u <- c(u, v)
  g <- exp(u) # x*u ~ N(0,x)
  cdf <- mean(g)
  cdf
}

## -----------------------------------------------------------------------------
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1
for (i in 1:m) {
  MC1[i] <- MC.Phi(R = 1000, antithetic = FALSE)
  MC2[i] <- MC.Phi(R = 1000, antithetic = TRUE)
}
c(mean(MC1),mean(MC2))
#theta估计：MC方法 对偶方法
c(sd(MC1)^2,sd(MC2)^2,sd(MC2)^2/sd(MC1)^2)
#MC方法 对偶方法 对偶方法方差/MC方法方差

## -----------------------------------------------------------------------------
pi <- 3.1415926
g <- function(x){
  return (x^2/sqrt(2*pi)*exp(-0.5*x^2))
}
f1 <- function(x){
  return (exp(1-x))
}
f2 <- function(x){
  return (1/(x^2))
  #return (1/sqrt(2*pi)*exp(-0.5*x^2))
}
m <- 1000
set.seed(123)
u <- runif(m)
x1 <- 1-log(1-u) # Inverse transform algorithm
x2 <- 1/(1-u) # Inverse transform algorithm
fg1 <- g(x1) / f1(x1)
fg2 <- g(x2) / f2(x2)
theta <- integrate(g,1,Inf)
theta1_hat <- mean(fg1)
theta2_hat <- mean(fg2)
sd1 <- sd(fg1)
sd2 <- sd(fg2)
c(theta[[1]], theta1_hat, theta2_hat) #(theta_true, theta_estimate_of_f1, theta_estimate_of_f2)
c(sd1,sd2) #(stand error with f1, stand error with f2)

## -----------------------------------------------------------------------------
x <- seq(1,5,0.1)
plot(x,g(x),type="l",ylim=c(0,1),ylab="prob",main="y-x")
lines(x,f1(x),col="red")
lines(x,f2(x),col="blue")
legend("topright",col=c("black","red","blue"),lty=1,legend=c("g","f1","f2"))

## -----------------------------------------------------------------------------
m <- 10000
N <- 50
k <- numeric(5)
theta_s <- numeric(5)
theta <- numeric(N)
var_s <- numeric(5)

set.seed(123)
g <- function(x) {
  exp(-x - log(1+x^2))
}
ab <- c(0,1/5,2/5,3/5,4/5,1)
for (i in 1:5)
  k[i] <- (1-exp(-1)) / (exp(-ab[i]) - exp(-ab[i+1]))


for (j in 1:N){
  for (i in 1:5){
    u <- runif(m)
    x <- -log(exp(-ab[i]) - u * (exp(-ab[i]) - exp(-ab[i+1])))
    tmpfg <- g(x) / (exp(-x) / (exp(-ab[i]) - exp(-ab[i+1])) )
    theta_s[i] <- mean(tmpfg)
    #var_s[i] <- sd(tmpfg)^2
  }
  theta[j] <- sum(theta_s) 
}
mean(theta)
var(theta)

## -----------------------------------------------------------------------------
c(0.5257801,0.0970314^2) #Estimate and variance from Example 5.13.

## -----------------------------------------------------------------------------
mcci <- function(m,n,miu, se){
  lower <- numeric(m)
  upper <- numeric(m)
  for (i in 1:m){
    Lx <- rnorm(n,miu,se)
    lower[i] <- mean(Lx) - qt(0.975, n-1)*sd(Lx)/sqrt(n)
    upper[i] <- mean(Lx) + qt(0.975, n-1)*sd(Lx)/sqrt(n)
  }
  L <- exp(mean(lower))
  R <- exp(mean(upper))
  rm(Lx) #clearing
  rm(upper)
  rm(lower)
  return (c(L, R))
}

## -----------------------------------------------------------------------------
mcci(1000,20,0,1)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
rm(X)
rm(Y)
return(as.integer(max(c(outx, outy)) > 5))
}
Ftest <- function(x,y){
  return (as.integer(var.test(x,y,alternative = "two.sided",conf.level = 0.945)[[3]]<0.055))
}

## -----------------------------------------------------------------------------
sigma1 <- 1
sigma2 <- 1.5
comparision <- function(m,n){
  num1 = num2 =0
  for (i in 1:m){
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    num1 = num1 + count5test(x, y)
  }
  for (i in 1:m){
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    num2 = num2 + Ftest(x, y)
  }
power1 <- num1/m
power2 <- num2/m
rm(x) #clear memory
rm(y) #clear memory
return (c(power1,power2))
}

## -----------------------------------------------------------------------------
print("CountFive Test, F Test")
comparision(1000,20) #small sample size
comparision(1000,200) #medium sample size
comparision(1000,2000) #large sample size

## -----------------------------------------------------------------------------
t <- matrix(c("n11","n21","n12","n22"),2,2)
colnames(t) <- c("reject under H1", "accept under H1")
rownames(t) <- c("reject under H1", "accept under H1")
t

## -----------------------------------------------------------------------------
library(boot)

## -----------------------------------------------------------------------------
x <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
B <- 1e4
set.seed(1)
thetastar <- numeric(B)
theta <- mean(x)
for(b in 1:B){
  xstar <- sample(x,replace=TRUE)
  thetastar[b] <- mean(xstar)
}


## -----------------------------------------------------------------------------
c(thetahat= mean(thetastar), theta= theta)
c(bias=mean(thetastar)-theta, se.boot=sd(thetastar), se.samp=sd(x)/sqrt(length(x)))
rm(thetastar)

## -----------------------------------------------------------------------------
mu <- 1
b <- 1
n <- 1e1
m <- 1e3

set.seed(12345)

boot.median <- function(t,i){
  return (median(t[i]))
}
boot.mean <- function(t,i){
  return (mean(t[i]))
}

ci.norm <- ci.basic <- ci.perc <- ci.bca <- matrix(NA,m,2)
de <- boot(data=x,statistic=boot.median, R = 999)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"), conf=0.95)
print(ci)

## -----------------------------------------------------------------------------
tmp <- rgamma(1000,shape = 12, scale = 1/12*108)
hist(tmp)
rm(tmp)

## -----------------------------------------------------------------------------
bootci_mc <- function(mu, sigma, n = 10, m = 10000){
  ci.norm <- ci.basic <- ci.perc <- matrix(NA, m, 2)
  prob.c <- prob.l <- prob.r <- c(0, 0, 0)
  for(i in 1:m){
    s <- rnorm(n, mu, sigma)
    de <- boot(data = s, statistic = boot.median, R = 999)
    ci <- boot.ci(boot.out = de, type = c("norm", "basic", "perc"))
    ci.norm[i, ] <- ci$norm[2:3];ci.basic[i, ] <- ci$basic[4:5];ci.perc[i, ] <- ci$percent[4:5]
    
    #ms <- mean(s)
    #print(ms)
    #if (ci.norm[, 1] <= ms & ci.norm[, 2] >= ms) prob.c[1] <- prob.c[1] +1
    #if (ci.basic[, 1] <= ms & ci.basic[, 2] >= ms) prob.c[2] <- prob.c[2] +1
    #if (ci.perc[, 1] <= ms & ci.perc[, 2] >= ms) prob.c[3] <- prob.c[3] +1
    
    #if (ci.norm[, 2] <= ms) prob.l[1] <- prob.l[1] +1
    #if (ci.basic[, 2] <= ms) prob.l[2] <- prob.l[2] +1
    #if (ci.perc[, 2] <= ms) prob.l[3] <- prob.l[3] +1
    
    #if (ci.norm[, 1] >= ms) prob.r[1] <- prob.r[1] +1
    #if (ci.basic[, 1] >= ms) prob.r[2] <- prob.r[2] +1
    #if (ci.perc[, 1] >= ms) prob.r[3] <- prob.r[3] +1
  }
  rm(s)
  rm(de)
  return(rbind(c(mean(ci.norm[, 1] <= mu & ci.norm[, 2] >= mu), mean(ci.norm[, 2] <= mu), mean(ci.norm[, 1] >= mu)), 
               c(mean(ci.basic[, 1] <= mu & ci.basic[, 2] >= mu), mean(ci.basic[, 2] <= mu), mean(ci.basic[, 1] >= mu)),
               c(mean(ci.perc[, 1] <= mu & ci.perc[, 2] >= mu), mean(ci.perc[, 2] <= mu), mean(ci.perc[, 1] >= mu))))
  #return (c(prob.c, prob.l, prob.r))
}

## -----------------------------------------------------------------------------
ci <- bootci_mc(mu = 0, sigma = 1, n = 50, m = 1000)
rownames(ci) <- c("normal", "basic", "perc")
colnames(ci) <- c("coverage.prob", "missed.left", "missed.right")

## -----------------------------------------------------------------------------
ci

## -----------------------------------------------------------------------------
library(bootstrap)
data <- bootstrap::scor

## -----------------------------------------------------------------------------
sigma_hat <- cov(data)
theta_hat <- eigen(sigma_hat)$values[1] / sum(eigen(sigma_hat)$values)

## -----------------------------------------------------------------------------
n <- nrow(data)
theta_jack <- numeric(n)
for (i in 1:n){
  tmp_scor <- cov(data[-i,])
  e <- eigen(tmp_scor)
  theta_jack[i] <- e$values[1] / sum(e$values)
}

## -----------------------------------------------------------------------------
rm(data)
rm(tmp_scor)

## -----------------------------------------------------------------------------
sigma_hat
theta_hat

bias <- (n-1) * (mean(theta_jack) - theta_hat)
se_theta <- sqrt((n-1) * mean((theta_jack - mean(theta_jack))^2))
c(bias,se_theta)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)

i <- 0
for (k in 1:(n-1)) {
  for (j in (k+1):n){
    #if (k==j) break
    y <- magnetic[-c(k,j)]
    x <- chemical[-c(k,j)]
    
    J1 <- lm(y ~ x)
    yhat1k <- J1$coef[1] + J1$coef[2] * chemical[k]
    yhat1j <- J1$coef[1] + J1$coef[2] * chemical[j]
    e1[i] <- ((magnetic[k] - yhat1k)^2 + (magnetic[j] - yhat1j)^2)/2
    
    J2 <- lm(y ~ x+ I(x^2))
    yhat2k <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
    yhat2j <- J2$coef[1] + J2$coef[2] * chemical[j] + J2$coef[3] * chemical[j]^2
    e2[i] <- ((magnetic[k] - yhat2k)^2 + (magnetic[j] - yhat2j)^2)/2
    
    J3 <- lm(log(y) ~ x)
    yhat3k <- exp(J3$coef[1] + J3$coef[2] * chemical[k])
    yhat3j <- exp(J3$coef[1] + J3$coef[2] * chemical[j])
    e3[i] <- ((magnetic[k] - yhat3k)^2 + (magnetic[j] - yhat3j)^2)/2

    J4 <- lm(log(y) ~ log(x))
    yhat4k <- exp(J4$coef[1] + J4$coef[2] * log(chemical[k]))
    yhat4j <- exp(J4$coef[1] + J4$coef[2] * log(chemical[j]))
    e4[i] <- ((magnetic[k] - yhat4k)^2 + (magnetic[j] - yhat4j)^2)/2
        
    i <- i+1
  }
}

## -----------------------------------------------------------------------------
c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
rm(e1,e2,e3,e4,yhat1j,yhat1k,yhat2k,yhat2j,yhat3k,yhat3j,yhat4k,yhat4j)

## -----------------------------------------------------------------------------
set.seed(123456)
x <- 1:10
y <- c(sort(rnorm(5)),sort(rnorm(5)))
cor_hat <- cor(x,y,method = "spearman")
p_cor.test <- cor.test(x,y, method = "spearman")[[3]]

## -----------------------------------------------------------------------------
m <- 2000
cor_t <- numeric(m)
for (i in 1:m){
  z <- sample(c(x,y))
  cor_t[i] <- cor(z[1:10],z[11:20],method = "spearman")
}
p_permutation <- mean(abs(cor_t)>=abs(cor_hat))
c(p_permutation, p_cor.test)
rm(x,y,cor_t)

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= exp(-abs(y)+abs(x[i-1])))
      x[i] <- y
    else{x[i] <- x[i-1]; k <- k+1}
  }
  return (list(x=x,k=k))
}

Gelman.Rubin <- function(psi){
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi,1,"var")
  W <- mean(psi.w)*(n-1)/n
  v.hat <- W * (n-1)/n + (B/(n))
  r.hat <- v.hat / W
  return (r.hat)
}

## -----------------------------------------------------------------------------
N <- 2000
sigma <- c(0.05, 0.5, 1, 5)
x0 <- 5
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
reject_rate <- c(rw1$k, rw2$k ,rw3$k, rw4$k)/N
reject_rate

## ----eval=FALSE---------------------------------------------------------------
#  #par(mfrow=c(2,2))
#  plot(rw1$x,type="l",xlab="sigma=0.05",ylab="X",ylim=range(rw1$x))
#  plot(rw2$x,type="l",xlab="sigma=0.5",ylab="X",ylim=range(rw2$x))
#  plot(rw3$x,type="l",xlab="sigma=1",ylab="X",ylim=range(rw3$x))
#  plot(rw4$x,type="l",xlab="sigma=5",ylab="X",ylim=range(rw4$x))

## -----------------------------------------------------------------------------
GR.converge <- function(sigma, k = 10, m = 5000){
  n.converge <- numeric(4)
    x0 <- rnorm(k, 0, 5)
    psi <- data.frame()
    for(j in 1:k)
      psi <- rbind(psi, rw.Metropolis(sigma, x0[j],  m)$x)
  
    for(n in 2:m){
      GR <- Gelman.Rubin(as.matrix(psi)[, 1:n])
      if(GR < 1.2)
        break
    }
  return (n)
}

## -----------------------------------------------------------------------------
set.seed(1)
cat("sigma=0.5:",GR.converge(sigma[2]),"\n")
cat("sigma=1:",GR.converge(sigma[3]),"\n")
cat("sigma=5:",GR.converge(sigma[4]),"\n")

## -----------------------------------------------------------------------------
rm(rw1,rw2,rw3,rw4)

## -----------------------------------------------------------------------------
gibbs <- function(x01, x02, N = 5000){

burn <- 1000 # burn-in length
X <- matrix(0, N, 2) # the chain, a bivariate sample

rho <- 0.9 
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(x01, x02) 
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]
return (X)
}

## -----------------------------------------------------------------------------
burn <- 1000
x <- gibbs(0,0)[(burn+1):N,]
plot(x[,1],type="l",xlab="x[1]",ylab="X",ylim=range(x[,1]))
plot(x[,2],type="l",xlab="x[2]",ylab="X",ylim=range(x[,2]))

## -----------------------------------------------------------------------------
colMeans(x)
cov(x)
cor(x)
plot(x, main="", cex=.5, xlab = bquote(X[1]), ylab = bquote(X[1]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
fit <- lm(x[,2]~x[,1])
summary(fit)

## -----------------------------------------------------------------------------
GR.biconverge <- function(m = 5000){
  k = 20
  x0 <- rnorm(k, 5, 5)
  psi1 <- psi2 <- data.frame()
  for(j in 1:k/2){
    xk <- gibbs(x0[j], x0[j+10],  m)
    psi1 <- rbind(psi1, xk[,1])
    psi2 <- rbind(psi2, xk[,2])
  }
  for(n in 10:m){
    GR1 <- Gelman.Rubin(as.matrix(psi1)[, 1:n])
    GR2 <- Gelman.Rubin(as.matrix(psi2)[, 1:n])
    #print(c(GR1,GR2))
    if(GR1 < 1.2 & GR2< 1.2)
      break
  }
  return (n)
}

## -----------------------------------------------------------------------------
GR.biconverge()

## -----------------------------------------------------------------------------
rm(x)

## -----------------------------------------------------------------------------
generating <- function(alpha, beta, gamma, N){
  x <- rnorm(N)
  eM <- rnorm(N)
  eY <- rnorm(N)
  m <- alpha*x + eM
  y <- beta*m + gamma*x + eY
  return(data.frame(x, m, y))
}

## -----------------------------------------------------------------------------
t_comp <- function(x, m, y){
  data <- data.frame(x, m, y)
  Model_M <- lm(m~x, data)
  Model_Y <- lm(y~m+x, data)
  t_stat <- Model_M$coef[2]*Model_Y$coef[2]/sqrt(Model_M$coef[2]^2*(summary(Model_M)$coefficients[2,2])^2+Model_Y$coef[2]^2*(summary(Model_Y)$coefficients[2,2])^2)
  return(as.numeric(t_stat))
}

permutation_1 <- function(data, N = 20, B = 99){
  reps <- numeric(B)
  for (i in 1:B) {
    k <- sample(1:N, size = N, replace = FALSE)
    x <- data[k, 1] 
    reps[i] <- t_comp(x, data$m, data$y)
  }
  p <- (sum(abs(reps) > abs(t_comp(data$x, data$m, data$y)))+1)/(1+B)
  return (p)
}

permutation_2 <- function(data, N = 20, B = 99){
  reps <- numeric(B)
  for (i in 1:B) {
    k <- sample(1:N, size = N, replace = FALSE)
    y <- data[k, 3]
    reps[i] <- t_comp(data$x, data$m, y)
  }
  p <- (sum(abs(reps) > abs(t_comp(data$x, data$m, data$y)))+1)/(1+B)
  return (p)
}

permutation_3 <- function(data, N = 20, B = 99){
  reps <- numeric(B)
  for (i in 1:B) {
    k <- sample(1:N, size = N, replace = FALSE)
    m <- data[k, 2]
    reps[i] <- t_comp(data$x, m, data$y)
  }
  p <- (sum(abs(reps)>abs(t_comp(data$x, data$m, data$y)))+1)/(1+B)
  return (p)
}

## -----------------------------------------------------------------------------
func <- function(alpha, beta, gamma=1, N=20){
  K <- 100
  p <- matrix(0, nrow = 3, ncol = K)
  for (i in 1:K) {
    data <- generating(alpha, beta, gamma, N)
    p[1,i] <- permutation_1(data)
    p[2,i] <- permutation_2(data)
    p[3,i] <- permutation_3(data)
  }
  return(c(mean(p[1,]<=0.05),mean(p[2,]<=0.05),mean(p[3,]<=0.05)))
}

## -----------------------------------------------------------------------------
set.seed(1)
simulation1 <- func(alpha=0, beta=0)
simulation2 <- func(alpha=0, beta=1)
simulation3 <- func(alpha=1, beta=0)

## -----------------------------------------------------------------------------
Permutation_Test <- data.frame(Model1=simulation1, Model2=simulation2, Model3=simulation3)
row.names(Permutation_Test) <- c("Permutation_Test_1", "Permutation_Test_2", "Permutation_Test_3")
Permutation_Test
#Type 1 error table

## -----------------------------------------------------------------------------
rm(Permutation_Test)

## -----------------------------------------------------------------------------
set.seed(12345)

#N <- 1e6; b1 <- 1; b2 <- 0.5; f0 <- 0.01
work <- function(N=1e6, b1=1, b2=0.5, b3=1, f0=0.01){
x1 <- rpois(N,1); x2 <- rexp(N,1); x3 <- sample(0:1,N,replace=TRUE)
g <- function(alpha){
tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
p <- 1/(1+tmp)
mean(p) - f0
}
solution <- uniroot(g,c(-20,0))
return (round(unlist(solution),5)[1:3])
}

## -----------------------------------------------------------------------------
work(1e6,0,1,-1,0.1)
work(1e6,0,1,-1,0.01)
work(1e6,0,1,-1,0.001)
work(1e6,0,1,-1,0.0001)

## -----------------------------------------------------------------------------
f00 <- c(0.1,0.01,0.001,0.0001)
x <- f00
for (i in 1:length(f00))
  x[i] <- work(1e6,0,1,-1,f00[i])
plot(y=f00,x=x)

## -----------------------------------------------------------------------------
rm(x,f00)

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
n = length(u)

## -----------------------------------------------------------------------------
L1 <- function(lam){
  
  t = 1
  t <- -sum(u) + sum( (v-u) / ( exp(lam*(v-u))-1 ))
  return (t)
}

## -----------------------------------------------------------------------------
uniroot(L1,c(0,30))[1]

## -----------------------------------------------------------------------------
iter <- 0
eps <- 1e-5
lambda_1 <- 1
lambda <- 2
while(1){
  lambda_1 <- lambda_1 - L1(lambda)/n 
  if (abs(1/lambda_1 - lambda) < eps) break
  iter <- iter+1
  lambda <- 1/lambda_1
}
c(iter = iter, lambda = lambda)

## -----------------------------------------------------------------------------
rm(u,v)
rm(lambda,lambda_1)

## -----------------------------------------------------------------------------
dim(c(1,2,3))

## -----------------------------------------------------------------------------
t

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
data1 <- data.frame(x=c(1,2,3),y=c(4,5,6))
data1
lapply(data1, scale01)

## ----warning=FALSE------------------------------------------------------------
library(tidyverse)
data2 <- data.frame(x=c(1,2,3),y=c(4,5,6),z=c("a","b","c"))
lapply(select_if(data2,is.numeric), scale01)

## -----------------------------------------------------------------------------
vapply(data1, sd, FUN.VALUE = numeric(1))

## -----------------------------------------------------------------------------
vapply( data2[,which(vapply(
  data2,is.numeric,FUN.VALUE = logical(1)))]
  , sd, FUN.VALUE = numeric(1))

## -----------------------------------------------------------------------------
rm(data1,data2)

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

## -----------------------------------------------------------------------------
cppFunction('NumericMatrix gibbs(double x01, double x02, int N = 1000){

NumericMatrix x(N+1000, 2);

double rho = 0.9 ;
double mu1 = 0, mu2 = 0;
double sigma1 = 1, sigma2 = 1;
double s1 = sqrt(1-rho*rho)*sigma1 ;
double s2 = sqrt(1-rho*rho)*sigma2 ;

x(0,0) = x01, x(0,1) = x02;

for (int i=1; i<N+1000; i++) {
  double x2 = x(i-1,1);
  double m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2;
  x(i,0) = rnorm(1, m1, s1)[0];
  double x1 = x(i,0);
  double m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1;
  x(i,1) = rnorm(1, m2, s2)[0];
}

return (x);
}')

## -----------------------------------------------------------------------------
gibbsR <- function(x01, x02, N = 2000){

burn <- 1000 # burn-in length
X <- matrix(0, N, 2) # the chain, a bivariate sample

rho <- 0.9 
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(x01, x02) 
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]
return (x)
}

## -----------------------------------------------------------------------------
x <- gibbs(100,10)[c(1001:2000),]
y <- gibbsR(100,10)
qqplot(x[,1], y[,1],main="qqplot of x1", xlab="from C", ylab="from R")
qqplot(x[,2], y[,2],main="qqplot of x2", xlab="from C", ylab="from R")

## -----------------------------------------------------------------------------
ts <- microbenchmark(gibbR=gibbsR(100,10),
gibbC=gibbs(100,10))
summary(ts)[,c(1,3,5,6)]


## -----------------------------------------------------------------------------
rm(x,y)

