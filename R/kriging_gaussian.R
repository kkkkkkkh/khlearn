# second_Step3: assume p = 2
#' Title a kernel with p=2
#'
#' @param x_sample matrix
#' @param y_sample vector
#'
#' @return list
#' @import stats
#' @export
#'

kriging_gaussian <- function(x_sample,y_sample){
  # x_sample 样本x
  # y_sample 样本y，sample为训练集
  sam_matrix <- x_sample
  sam_y <- y_sample
  k <- length(sam_matrix[1,])
  t1 <- proc.time()
  theta <- matrix(runif(2*k*k,min=0,max=2),ncol = k)
  len_sam <- length(sam_y)
  # step_1抽了n2个样本
  # 这个版本的kriging 模型的核心就是优化上述参数theta。
  # 初始化相关阵
  # 设置循环迭代次数，开始循环之前，确定levy flight。
  # https://blog.csdn.net/zyqblog/article/details/80905019
  Loop <- 50*k #最大迭代次数
  loop <- 1
  sigma_mu <- (gamma(2.5)/2^0.5/gamma(1.25)/1.5/2^0.25)^(2/3) #莱维飞行的随机数
  r_ori <- matrix(nrow = len_sam,ncol = len_sam)
  r_ori[is.na(r_ori)] <- 0#初始化相关阵
  t2 <- sample(1:2*k,1)
  theta_ori <- abs(theta[t2,]) #随机选择theta作为起始迭代theta
  for (i in 1:len_sam){
    for (j in i:len_sam){
      r_ori[i,][j] <- exp(-1*as.numeric(sum(theta_ori*(abs(sam_matrix[i,]-sam_matrix[j,]))^2)))
    }
  }
  r_ori <- r_ori + t(r_ori)
  frac1 <- det(r_ori)^(1/(len_sam-1)) #目标函数（似然函数转化来的）第一部分
  beta <- 1/sum(rowSums(solve(r_ori)))*colSums(solve(r_ori))%*%sam_y # 均值
  betavec <- vector(length = len_sam)
  for (i in 1:len_sam) {
    betavec[i] <- beta
  }
  y_beta <- sam_y-betavec # y_beta <- y-beta
  # 似然函数第二部分
  frac2 <- as.numeric(1/(len_sam-1)*abs(as.numeric(t(y_beta)%*%solve(r_ori)%*%(y_beta))))
  aim_func <- frac1*frac2 #目标函数，迭代的目的是通过最小化aim_func来求得theta
  while (loop<Loop) {
    #levy flight
    levy2 <- vector(length = k) # 莱维飞行
    for (i in 1:k){
      ss <- rnorm(1,mean = 0,sd = sigma_mu)/(abs(rnorm(1,mean = 0,sd = 1))^(2/3))
      levy2[i] <- ss
    }
    t4 <- sample(1:2*k-1,1)
    theta2 <- abs(theta[-t2,][t4,]+levy2) #除了初始选中的一行外，再选一行
    rr <- matrix(nrow = len_sam,ncol = len_sam)
    rr[is.na(rr)] <- 0 # 再次计算相关阵以及目标函数
    for (i in 1:len_sam){
      for (j in i:len_sam){
        rr[i,][j] <- exp(-1*as.numeric(sum(theta2*(abs(sam_matrix[i,]-sam_matrix[j,]))^2)))
      }
    }
    rr <- rr+t(rr)
    frac11 <- abs(det(rr))^(1/(len_sam-1))
    beta2 <- 1/sum(rowSums(solve(rr)))*colSums(solve(rr))%*%sam_y
    betavec2 <- vector(length = len_sam)
    for (i in 1:len_sam) {
      betavec2[i] <- beta2
    }
    y_beta2 <- sam_y-betavec2
    frac22 <- as.numeric(1/(len_sam-1)*abs(as.numeric(t(y_beta2)%*%solve(rr)%*%(y_beta2))))
    aim_func2 <- frac11*frac22
    if (aim_func2<aim_func){
      # 说明目标函数经过刚才的计算有改进，将新的theta(加上levy flight的)
      # 替换掉初始的theta
      theta[t2,] <- theta2
    }
    pa <- 0.7 #以概率0.7更新k个不好的点
    aim_vec <- vector(length = 2*k)
    for (i in 1:(2*k)){
      rl <- matrix(nrow = len_sam,ncol = len_sam)
      theta <- abs(theta)
      theta_l <- theta[i,]
      rl[is.na(rl)] <- 0
      for (o in 1:len_sam){
        for (j in o:len_sam){
          rl[o,][j] <- exp(-1*as.numeric(sum(theta_l*(abs(sam_matrix[o,]-sam_matrix[j,]))^2)))
        }
      }
      rl <- rl+t(rl)
      frac1 <- det(rl)^(1/(len_sam-1))
      beta <- 1/sum(rowSums(solve(rl)))*colSums(solve(rl))%*%sam_y
      betavec <- vector(length = len_sam)
      betavec[is.na(betavec)] <- beta
      y_beta <- sam_y-betavec
      frac2 <- as.numeric(1/(len_sam-1)*abs(as.numeric(t(y_beta)%*%solve(rl)%*%(y_beta))))
      aim_func <- frac1*frac2
      aim_vec[i] <- aim_func
    }
    aim_vec_list <- order(aim_vec)
    theta_aim <- theta[aim_vec_list[1],] # 选取最优点为theta_aim
    for (i in 1:(2*k-1)){
      randu <- runif(1,min = 0,max = 1)
      if (randu<pa){
        ran3 <- 0.1
        gg <- sample(1:(2*k),1)
        o3 <- theta[gg,]
        o4 <- theta[-gg,][sample(1:(2*k-1),1),]
        theta[aim_vec_list[1+i],] <- abs(theta[aim_vec_list[1],]+ran3*(o3-o4))
      }
    }
    aim_func <- aim_vec[aim_vec_list[1]]
    t2 <- sample(1:(2*k),1)
    loop <- loop+1
  }
  t2 <- proc.time()
  theta_aim <- as.vector(theta_aim)
  r_hat <- matrix(nrow = len_sam,ncol = len_sam)
  r_hat[is.na(r_hat)] <- 0
  for (i in 1:len_sam){
    for (j in i:len_sam){
      r_hat[i,][j] <- exp(-1*as.numeric(sum(theta_aim*(sam_matrix[i,]-sam_matrix[j,])^2)))
    }
  }
  r_hat <- r_hat+t(r_hat)
  beta_hat <- 1/sum(rowSums(solve(r_hat)))*colSums(solve(r_hat))%*%sam_y
  betavec_hat <- vector(length = len_sam)
  for (i in 1:len_sam) {
    betavec2[i] <- beta_hat
  }
  y_beta_hat <- sam_y-betavec_hat
  sigma_hat <- as.numeric(1/(len_sam-1)*abs(as.numeric(t(y_beta_hat))%*%solve(r_hat)%*%((y_beta_hat))))
  return(list(r_hat=r_hat,
              beta_hat=beta_hat,
              sigma_hat=sigma_hat,
              theta_aimm=theta_aim,
              time=t2-t1,
              xsam=sam_matrix,
              ysam=sam_y))
}
