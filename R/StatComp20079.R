#' @title Lagrange interpolation .
#' @description Approximate unknown function with Lagrange interpolation
#' @param x The x value corresponding to the missing value y (numeric)
#' @param X  The x value vector corresponding to the Y vector(numeric)
#' @param Y All known values of y (numeric)
#' @return A vector (y)
#' @examples
#' \dontrun{
#' X <- c(-2.15,-1.00,0.01,1.02,2.03,3.25)
#' Y <- c(17.03,7.24,1.05,2.03,17.06,23.05)
#' f <- function(x)
#'  Lagrange(x,X,Y)
#' curve(f(x),-2.15,3.25)
#' }
#' @export
Lagrange<-function(x,X,Y){
  n<-length(X) 
  lagr<-0 
  if(length(X)!=length(Y))
    stop("input error,the length of input values X and Y should be equal")
  if(n<2)
    stop("the length of X and Y should be bigger than 1")
  for(i in 1:n)
  {
    li<-1
    for(j in 1:n)
    {
      if(i!=j)
        li<-li*(x-X[j])/(X[i]-X[j])
    }
    lagr<-li*Y[i]+lagr
  }
  return(lagr)
}

#' @title FCM
#' @description Fuzzy cluster analysis of data
#' @param dat training set(matrix)
#' @param k Number of clusters(numeric)
#' @return a list (u,v,cluster_id)
#' @importFrom stats rnorm rgamma var 
#' @useDynLib StatComp20079
#' @examples
#' \dontrun{
#' data(iris)
#' dat <- iris[,2:4]
#' result <- FCM(dat,3)
#' plot(iris$Sepal.Width,iris$Petal.Width,col=result$cluster_id,main="Visualization results",pch=19)
#' }
#' @export
FCM <- function(dat, k){
  m <- 2
  iter.max <- 1000
  eps <- 1e-06
  dat <- as.matrix(dat)
  set.seed(6386)
  v0 <- dat[sample(nrow(dat), k),] 
  
  n <- nrow(dat)
  d <- matrix(0, n, k)
  for (j in 1:k) {
    d[, j] = sqrt(rowSums(sweep(dat, 2, v0[j, ], "-")^2))
  }
  
  u <- matrix(NA, n, k)
  for (j in 1:k){
    for (i in 1:n){
      if (any(d[i,] == 0)){
        u[i, ] <- rep(1/k, k)
      }else{
        u[i, j] <- 1/(sum((d[i, j]/d[i, ])^(2/(m - 1))))
      }
    }
  }
  v <- t(u^m) %*% dat/colSums(u^m)
  J <- sum(d^2 * (u^m)) 
  
  iter <- 0
  j.best <- Inf
  while((iter < iter.max) && (abs(j - j.best) > eps)){
    j.best <- J
    for (j in 1:k) {
      d[, j] = sqrt(rowSums(sweep(dat, 2, v[j, ], "-")^2))
    }
    
    for (j in 1:k){
      for (i in 1:n){
        if (any(d[i,] == 0)){
          u[i, ] <- rep(1/k, k)
        }else{
          u[i, j] <- 1/(sum((d[i, j]/d[i, ])^(2/(m - 1))))
        }
      }
    }
    v <- t(u^m) %*% dat/colSums(u^m)
    J <- sum(d^2 * (u^m))
    iter = iter + 1
  }
  cluster_id <-  apply(u, 1, which.max)
  result <- list()
  result$u <- u
  result$v <- v
  result$cluster_id <- cluster_id
  return(result)
}