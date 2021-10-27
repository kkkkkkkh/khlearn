#' This function is about to build a lhd
#' @title LHD
#' @description This function is about to build a lhd
#' @details no
#' @param num number of design point
#' @param fac k
#' @param method maxmin means maximum the minimum distance
#'               mincorr means minimum the column correlation
#'
#' @return latin hyper design
#'
#' @import stats
#'
#' @export
LHD <- function(num,fac,method=c('maxmin','mincorr')){
  # num 每一个维度要划分的个数
  # fac 因素的个数
  #这个函数目的是形成特定准则下的拉丁超立方体设计
  # maxmin:最大化最小距离准则
  # mincorr：最小化因子相关性准则
  call <- match.call()
  method <- match.arg(method)
  loop <- 1
  fa <- matrix(nrow=num,ncol=fac)
  while (loop<=fac) {
    fa[c(1:num),loop] <- c(1:num)
    fa[c(1:num),loop] <- sample(fa[c(1:num),loop],length(fa[c(1:num),loop]))
    loop <- loop+1
  }
  #MMA算法
  if (method=='maxmin'){
    #迭代循环找出最大最小距离的设计：
    loop <- 1
    Loop <- 1000
    len <- matrix(nrow=num,ncol=num)
    yitalist <- integer(Loop)
    for (loop in 1:Loop) {
      for (roun in 1:num){
        for (loo in 1:roun){
          len[loo,roun] <- sum(abs(fa[loo,1:fac]-fa[roun,1:fac]))
          len[roun,loo] <- sum(abs(fa[loo,1:fac]-fa[roun,1:fac]))
        }
      }
      #计算距离标准yita
      lenwiot0 <- len[len!=0]
      lenwiot0 <- 1/lenwiot0
      yita <- sum(lenwiot0/2)
      #计算最小距离
      maxlen <- which(len==len[which.max(len)],arr.ind=T)
      dmax <- len[maxlen[1,1],maxlen[1,2]]
      for (numb in 1:num){
        len[numb,numb] <- dmax
      }
      minlen <- which(len==len[which.min(len)],arr.ind=T)
      #显示最大化后的最小距离
      lenmin <- len[minlen[1,1],minlen[1,2]]
      #随机交换LHD某一列的某两个数字
      colu <- sample(1:fac,1,FALSE)
      ro <- fa[1:num,colu]
      rosam <- sample(1:num,2,FALSE)
      roo <- ro
      roo[rosam[1]]=ro[rosam[2]]
      roo[rosam[2]]=ro[rosam[1]]
      fa_new <- fa
      fa_new[1:num,colu] <- roo
      #计算新的距离标准
      len_new <- matrix(nrow=num,ncol=num)
      for (roun in 1:num){
        for (loo in 1:roun){
          len_new[loo,roun] <- sum(abs(fa_new[loo,1:fac]-fa_new[roun,1:fac]))
          len_new[roun,loo] <- sum(abs(fa_new[loo,1:fac]-fa_new[roun,1:fac]))
        }
      }
      lenwiot0_new <- len_new[len_new!=0]
      lenwiot0_new <- 1/lenwiot0_new
      yita_new <- sum(lenwiot0_new/2)
      maxlen_new <- which(len_new==len_new[which.max(len_new)],arr.ind=T)
      dmax_new <- len[maxlen_new[1,1],maxlen_new[1,2]]
      for (numb in 1:num){
        len_new[numb,numb] <- dmax_new
      }
      minlen_new <- which(len_new==len_new[which.min(len_new)],arr.ind=T)
      #显示最大化后的最小距离
      lenmin_new <- len_new[minlen_new[1,1],minlen_new[1,2]]
      if (lenmin_new>lenmin){
        fa <- fa_new
        len <- len_new
        yita <- yita_new
        yitalist[loop] <- yita
      }
      else {
        prob <- exp((yita-yita_new)*(200*yita))
        rand <- runif(1)
        if (prob>rand){
          fa <- fa_new
          len <- len_new
          yita <- yita_new
          yitalist[loop] <- yita
        }
        else {
          fa <- fa
          len <- len
          yita <- yita
          yitalist[loop] <- yita
        }
      }
    }
    maxlen <- which(len==len[which.max(len)],arr.ind=T)
    dmax <- len[maxlen[1,1],maxlen[1,2]]
    for (numb in 1:num){
      len[numb,numb] <- dmax
    }
    #找到之间距离最小的两个点
    minlen <- which(len==len[which.min(len)],arr.ind=T)
    #显示最大化后的最小距离
    lenmin <- len[minlen[1,1],minlen[1,2]]
    roumatrix <- matrix(ncol=fac,nrow=fac)
    for (corr in 1:fac) {
      for (cor in 1:corr) {
        #计算相关系数
        rou <- sum((fa[1:num,corr]-num*(num-1)/2)*(fa[1:num,cor]-num*(num-1)/2))/sqrt((sum((fa[1:num,corr]-num*(num-1)/2)^2))*sum((fa[1:num,cor]-num*(num-1)/2)^2))
        roumatrix[corr,cor]=rou^2
        roumatrix[cor,corr]=rou^2
      }
    }
    rouaim <- sum(roumatrix[roumatrix!=1],na.rm=TRUE)/(fac*(fac-1))
    lenmatrix <- len
    lenmatrix2 <- len
    diag(lenmatrix) <- 0
    diag(lenmatrix2) <- NA
    lenmatrix2 <- t(matrix(t(lenmatrix2)[which(!is.na(lenmatrix2))],nrow=num-1,ncol=num))
    len_table <- table(lenmatrix2)/2
    lhd <- fa
    # 返回值 yitalist:目标函数每次迭代值
    # lenmatrix 距离矩阵
    # lenmin 最小距离
    # len_table 距离汇总
    # rouaim 相关系数大小
    # lhd 迭代所得拉丁超立方体设计
    # return (list(lenmatrix,lenmin,len_table,rouaim,lhd))
    return (list(yitalist[Loop],rouaim,lhd))
  }
  if (method=='mincorr'){
    loop <- 1
    Loop <- 10000
    len <- matrix(nrow=num,ncol=num)
    roumatrix <- matrix(ncol=fac,nrow=fac)
    roumatrix_new <- matrix(ncol=fac,nrow=fac)
    roulist <- integer(Loop)
    for (loop in 1:Loop){
      for (corr in 1:fac) {
        for (cor in 1:corr) {
    #计算相关系数
            rou <- sum((fa[1:num,corr]-num*(num-1)/2)*(fa[1:num,cor]-num*(num-1)/2))/sqrt((sum((fa[1:num,corr]-num*(num-1)/2)^2))*sum((fa[1:num,cor]-num*(num-1)/2)^2))
            roumatrix[corr,cor]=rou^2
            roumatrix[cor,corr]=rou^2
        }
      rouaim <- sum(roumatrix[roumatrix!=1],na.rm=TRUE)/(fac*(fac-1))
      #随机交换某列中的两个数字
      colu <- sample(1:fac,1,FALSE)
      ro <- fa[1:num,colu]
      rosam <- sample(1:num,2,FALSE)
      roo <- ro
      roo[rosam[1]]=ro[rosam[2]]
      roo[rosam[2]]=ro[rosam[1]]
      fa_new <- fa
      fa_new[1:num,colu] <- roo
      #再次计算相关系数
      for (corr in 1:fac) {
        for (cor in 1:corr) {
          #计算相关系数
          rou_new <- sum((fa_new[1:num,corr]-num*(num-1)/2)*(fa_new[1:num,cor]-num*(num-1)/2))/sqrt((sum((fa_new[1:num,corr]-num*(num-1)/2)^2))*sum((fa_new[1:num,cor]-num*(num-1)/2)^2))
          roumatrix_new[corr,cor]=rou_new^2
          roumatrix_new[cor,corr]=rou_new^2
        }
      }
      rouaim_new <- sum(roumatrix_new[roumatrix_new!=1],na.rm=TRUE)/(fac*(fac-1))
      #MMA
      if (abs(rouaim)>=abs(rouaim_new)){
        fa <- fa_new
        roumatrix <- roumatrix_new
        roulist[loop] <- rouaim_new
      }
      else {
        prob <- exp(-1*abs(rouaim-rouaim_new)*(2000*rouaim))
        rand <- runif(1)
        if (prob>rand){
          fa <- fa_new
          roumatrix <- roumatrix_new
          roulist[loop] <- rouaim_new
        }
        else {
          fa <- fa
          roumatrix <- roumatrix
          roulist[loop] <- rouaim
              }
            }
      }
    }
    for (roun in 1:num){
      for (loo in 1:roun){
        len[loo,roun] <- sum(abs(fa[loo,1:fac]-fa[roun,1:fac]))
        len[roun,loo] <- sum(abs(fa[loo,1:fac]-fa[roun,1:fac]))
          }
        }
    #计算最小距离
    maxlen <- which(len==len[which.max(len)],arr.ind=T)
    dmax <- len[maxlen[1,1],maxlen[1,2]]
    for (numb in 1:num){
      len[numb,numb] <- dmax
    }
    minlen <- which(len==len[which.min(len)],arr.ind=T)
    #显示最大化后的最小距离
    lenmin <- len[minlen[1,1],minlen[1,2]]
    lenmatrix2 <- len
    diag(lenmatrix2) <- NA
    lenmatrix2 <- t(matrix(t(lenmatrix2)[which(!is.na(lenmatrix2))],nrow=num-1,ncol=num))
    len_table <- table(lenmatrix2)/2
    lhd <- fa
    # roulist:目标函数每次迭代的结果
    # roumatrix ：相关系数阵
    return(list(roumatrix,rouaim,lenmin,len_table,lhd))
  }
}
