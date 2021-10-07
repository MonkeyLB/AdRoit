#' Initial estimation of the proportions of cell types
#'
#' @param x the reference data matrix (gene x cell type) of estimated gene expressions (mu)
#' @param y a mixture sample to be deconvoluted
#' @param w0 gene weights related to cell type specificity
#' @param w gene weights related to cross sample varibility
#' @param method method to perform the initial estimation of the cell proportions. Available choices: 'nnls'(default), 'mle'
#'
#' @return a vector of unnormalized estimated cell type proportions
#' @export

initBetaEst<-function(x,y,w0,w,method="nnls"){
  res=NULL
  if(method=="nnls"){
    M = nnls::nnls(w0*w*x,  w0*w*y)
    res = M$x
  }else if (method=="mle"){
    theta0 = diff(sort(c(runif(ncol(x) - 1), 0, 1)))
    fn = function(theta) {
      (w0 * w) %*% (y - (x %*% theta))^2
    }
    gn = function(theta) {
      -2 * t(x) %*% (w0 * w * (y - (x %*% theta)))
    }

    op <- options(show.error.messages = FALSE)
    on.exit(op)
    esti = try(optim(theta0, fn, gr = gn,
                     lower = rep(0, ncol(x)),
                     upper = rep(Inf, ncol(x)), method = "L-BFGS-B"),
               silent = T)
    if (class(esti) == "try-error") {
      res = (rep(NA, length(theta0)))
    } else {
      res = (esti$par)
    }

  }
  return(as.numeric(res))

}
