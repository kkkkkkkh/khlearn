#' Title Partial Least Square Regression
#'
#' @param x matrix
#' @param y vector
#'
#' @return just loadings
#' @import stats
#' @export
#'

PLSR1 <- function(x,y){
  # 用5折交叉验证的方法确定偏最小二乘回归主成分的个数
  # 然后用全体x、y求解对应的w*
  # PLSR1
  x <- scale(x)
  y <- scale(y)
  pis <- round(length(y)/5)
  a <- vector(length = 6)
  a[1] <- 1
  for (i in 1:5){
    a[i+1] <- i*pis
  }
  a[6] <- length(y)
  k <- length(x[1,])
  rmse <- vector(length = k)
  rmse[is.na(rmse)] <- 0
  for (i in 1:5){
    xx <- x[-c((a[i]):(a[i+1]-1)),]
    xtest <- x[c((a[i]):(a[i+1]-1)),]
    yy <- y[-c((a[i]):(a[i+1]-1))]
    ytest <- y[c((a[i]):(a[i+1]-1))]
    ypre <- vector(length = length(ytest))
    ypre[is.na(ypre)] <- 0
    w <- matrix(nrow = k,ncol = k)
    r <- vector(length = k)
    p <- matrix(nrow = k,ncol = k)
    for (j in 1:k){
      xp <- xx
      yp <- yy
      mat <- t(xp)%*%yp%*%t(yp)%*%xp
      en <- eigen(mat)
      w[,j] <- as.numeric(en$vec[,1])
      t <- xp%*%w[,j]
      p[,j] <- (t(xp)%*%t)/sum(t^2)
      r[j] <- t(yp)%*%t/sum(t^2)
      xx <- xp-t%*%t(p[,j])
      yy <- yp-t*(r[j])
    }
    for (e in 1:k){
      if (e==1){
        y_hat <- xtest%*%(w[,c(1:1)]*r[1])
        rmse[e] <- rmse[e] + mean(abs((ytest-y_hat)/ytest))
      }
      if (e>1){
        y_hat <- xtest%*%(w[,c(1:e)]%*%(r[c(1:e)]))
        rmse[e] <- rmse[e] + mean(abs((ytest-y_hat)/ytest))
      }
    }
  }
  rmse <- rmse/5
  or <- which.min(rmse)
  ww <- matrix(nrow = k,ncol = or)
  pp <- matrix(nrow = k,ncol = or)
  rr <- vector(length = or)
  for (d in 1:or){
    mat <- t(x)%*%y%*%t(y)%*%x
    en <- eigen(mat)
    ww[,d] <- as.numeric(en$vec[,1])
    t <- x%*%ww[,d]
    pp[,d] <- (t(x)%*%t)/sum(t^2)
    rr[d] <- t(y)%*%t/sum(t^2)
    x <- x-t%*%t(pp[,d])
    y <- y-t*(rr[d])
  }
  loadings <- ww%*%solve(t(pp)%*%ww)
  return(list(rmse,loadings))
}
