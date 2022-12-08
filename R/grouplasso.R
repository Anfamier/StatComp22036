#' @title Data generating for illustration of group lasso 
#' @name g_data
#' @description A dataset used to illustrate the performance of \code{g_lasso}
#' @import MASS
#' @examples
#' \dontrun{
#' 
#' set.seed(114514)
#' matr = matrix(nrow = 15, ncol = 15)
#' for(i in 1:15)
#'   for(j in 1:15)
#'     matr[i, j] = 0.5 ^ abs(i-j)
#' mu = rep(0, 15)
#' Z = mvrnorm(n = 50, mu, matr)
#' X = matrix(0, nrow = 50, ncol = 30)
#' for(i in 1:50)
#'   for(j in 1:15){
#'     if(Z[i, j] < qnorm(1/3))
#'       X[i, 2*j-1] = 1
#'     else if(Z[i, j] < qnorm(2/3))
#'       X[i, 2*j] = 1
#'   }
#' coef = rep(0, 30)
#' coef[2] = 1.8
#' coef[1] = -1.2
#' coef[6] = 1
#' coef[5] = 0.5
#' coef[10] = 1
#' coef[9] = 1
#' 
#' Y = X %*% coef
#' sigma = sqrt(var(Y) * 49 / 10 / 50)
#' Y = Y + rnorm(50, 0, sigma)
#' 
#' qr = qr(scale(X, scale = FALSE))
#' X = qr.Q(qr)
#' Y = scale(Y, scale = FALSE)
#' 
#' p <- rep(2,15)
#' J = 15
#' 
#' index <- numeric(length(p))
#' for (i in 1:J)
#'   index[i] <- sum(p[1:i])
#' 
#' X1 <- X[, 1:index[1]]
#' lambda_upp <- norm_2(t(X1) %*% Y)/sqrt(p[1])
#' for (j in 2:J){
#'   Xj <- X[, (index[j-1]+1):index[j]]
#'   lambda_upp = max(lambda_upp, norm_2(t(Xj) %*% Y)/sqrt(p[1]))
#' }
#' 
#' M <- 50
#' lam <- seq(0,lambda_upp,lambda_upp / M) #参数取值范围 10为参数取值数量
#' beta_s <- matrix(sum(p)*(M+1), sum(p), M+1) #p为pj向量
#' beta_s[,1] <- rep(0,sum(p)) #beta存储阵 第一列为0
#' 
#' for (j in 1:M){
#'   betaj <- g_lasso(X,Y,lam[j],p,J, index)
#'   beta_s[,j+1] <- betaj
#' }
#' }
NULL

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

#' @title group lasso function for variable selection
#' @description group lasso function for variable selection
#' @param X the design matrix with factor form
#' @param Y the response variable vector
#' @param lambda tuning parameter for penalty in group lasso
#' @param p numeric vector for factor sizes
#' @param J number of factors
#' @param index index for element position of factors, computed from p
#' @return estimate of parameter beta in group lasso model
#' @export
g_lasso <- function(X, Y, lambda, p, J, index){
  beta <- rep(0, sum(p))
  while(1){
    beta1 <- beta
    X1 <- X[,1:index[1]]
    S1 <- t(X1) %*% (Y - X %*% c(rep(0,p[1]),beta[(index[1]+1):index[J]]) )
    beta[1:p[1]] <- pos(1 - lambda * sqrt(p[1]) / norm_2(S1)) * S1
    
    for (j in 2:J){
      Xj <- X[, (index[j-1]+1):index[j]]
      
      if (j<J) Sj <- t(Xj) %*% (Y - X %*% c(beta[1:index[j-1]],rep(0,p[j]),beta[(index[j]+1):index[J]]) )
      else Sj <- t(Xj) %*% (Y - X %*% c(beta[1:index[j-1]],rep(0,p[j]) ))
      
      beta[(index[j-1]+1) : index[j]] <- pos(1 - lambda * sqrt(p[j]) / norm_2(Sj)) * Sj
    }
    eps <- norm_1(beta - beta1)
    if (abs(eps) < 1e-6) break
  }
  return (beta)
}

