#' Projection Covariance
#'
#' Calculate the projection covariance of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @return The projection covariance of X and Y
#' @examples X = matrix(rnorm(100*34,1),100,34)
#' @examples Y = matrix(rnorm(100*62,1),100,62)
#' @examples pcov(X,Y)
#' @seealso \code{\link{pcor}} \code{\link{pcor.test}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829–843. https://doi.org/10.1093/biomet/asx043
#' @export
pcov <- function(X,Y) {
  pcov_va(X,Y)
}


#' Projection Correlation
#'
#' Calculate the projection correlation of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @return The projection correlation of X and Y
#' @examples X = matrix(rnorm(100*34,1),100,34)
#' @examples Y = matrix(rnorm(100*62,1),100,62)
#' @examples pcor(X,Y)
#' @seealso \code{\link{pcov}} \code{\link{pcor.test}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829–843. https://doi.org/10.1093/biomet/asx043
#' @export
pcor <- function(X,Y) {
  pcor_va(X,Y)
}


#' Projection Correlatoin Permutation Test
#'
#' Return the test result of projection correlation test. Test whether two vectors are independnet. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @return \code{pcor.value} projection correlation of X and Y
#' @return \code{p.value} the p-value under the null hypothesis two vectors are independent
#' @examples X = matrix(rnorm(10*3,1),10,3)
#' @examples Y = matrix(rnorm(10*6,1),10,6)
#' @examples pcor.test(X,Y)
#' @seealso \code{\link{pcov}} \code{\link{pcor}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829–843. https://doi.org/10.1093/biomet/asx043
#' @export
pcor.test=function(X,Y){

  n <- nrow(X)
  value <- t_va(X,Y)
  t.value <- n*value

  times <-  2000
  pcor.permu <- replicate(times, t_va(X[sample(1:n,n),],Y), simplify = TRUE)
  p.value <- 1-mean(value>pcor.permu)

  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(Y)))
  method <- "Projection Correlation Permutation Test of Independence"
  rval <- list(method = method, data.name = dname, t.value = t.value, p.value = p.value)

  return(rval)
}


