#' Title predictor of KPLS AND guassian kriging
#'
#' @param model a model from KPLS or guassian kriging
#' @param input x need to be create y_hat
#'
#' @return the prediction of x
#' @import stats
#' @export
#'

Predictor <- function(model,input){
  xsam <- model$xsam
  ysam <- model$ysam
  r_hat <- model$r_hat
  beta_hat <- model$beta_hat
  sigma_hat <- model$sigma_hat
  theta_aimm <- model$theta_aimm
  beta_hat_vect <- vector(length = length(ysam))
  for (i in 1:length(ysam)){
    beta_hat_vect[i] <- as.numeric(beta_hat)
  }
  r0_matrix <- matrix(nrow = length(input[,1]),ncol = length(xsam[,1]))
  r0_matrix[is.na(r0_matrix)] <- 0
  for (i in 1:length(input[,1])){
    for (j in 1:length(xsam[,1])){
      r0_matrix[i,][j] <- exp(-1*as.numeric(sum((theta_aimm*(input[i,]-xsam[j,])^2))))
    }
  }
  pre_vect <- vector(length = length(input[,1]))
  pre_vect[is.na(pre_vect)] <- 0
  for (i in 1:length(input[,1])){
    pre_vect[i] <- as.numeric(beta_hat)+t(r0_matrix[i,])%*%solve(r_hat)%*%(as.vector(ysam-beta_hat_vect))
  }
  return(list(pre_vect))
}
