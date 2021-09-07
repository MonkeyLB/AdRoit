#' Estimate the proportions for a single sample which is left out in training data
#'
#' This function estimates the proportions of a single bulk sample while being left out when estimating gene cross-sample variations and cross-platform biases.
#'
#' @param n index of sample (column) to be estimated. Avalibility of multiple bulk sample replicates is assumed 
#' @param bulk.sample a matrix or data.frame with the rows being genes and columns being samples (or spatial spots for spatial transcriptomics data).
#' @param single.ref the reference object built by the function `ref.build()`.
#' @param silent whether to print out messages. Default is FALSE.
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom stats optim
#'
#'
#' @return the cell type proportions for the nth sample.
#' @export

AdRoit.est.loo <- function(n, bulk.sample, single.ref, per.sample.adapt=FALSE,silent = FALSE){
  
  i <- NULL
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  registerDoParallel(no_cores)
  
  genes = rownames(single.ref[[1]])
  
  x = single.ref[[1]]
  z = single.ref[[2]]
  w0 = single.ref[[3]]
  nb = ncol(bulk.sample)
  if(is.null(nb)){
    bulk.sample=cbind(bulk.sample)
    nb=1
  }
  
  
  if (nb < 3) {
    if (is.null(single.ref[[4]])){
      stop(print("no enough bulk replicates; please first use single cell data to estimate cross sample variability"))
    }
    w.mu = single.ref[[4]][[1]]
    w.sigma = single.ref[[4]][[2]]
  } else {
    tmp = foreach(i = genes, .combine = rbind) %dopar%
      negbin.est(as.integer(bulk.sample[i, -n]))
    rownames(tmp) = genes
    w = 1/(1 + tmp[, 2]/tmp[, 1])
    w[which(is.infinite(w) | is.na(w))] = 0
    if(!per.sample.adapt){
      M = nnls::nnls(w0*x, w0*rowMeans(bulk.sample[genes, -n]))
      ptheta = M$x/sum(M$x)
      msf = abs(rowMeans(bulk.sample[genes, -n])/(x %*% ptheta))
      msf[which(is.infinite(msf) | is.na(msf))] = median(msf, na.rm = T)
    }
  }
  
  if (nb < 3 || per.sample.adapt) {
    M = nnls::nnls(w0*x, w0*(bulk.sample[genes, n]))
    ptheta = M$x/sum(M$x)
    msf = abs((bulk.sample[genes, n])/(x %*% ptheta))
    msf[which(is.infinite(msf) | is.na(msf))] = median(msf, na.rm = T)
  }
  if (nb < 3){
    tmp=cbind(matrix(M$x %*% t(w.mu),ncol=1), matrix(M$x^2 %*% t(w.sigma),ncol=1))
    rownames(tmp)=genes
    w=1/(1 + tmp[, 2]/tmp[, 1])
    w[which(is.infinite(w) | is.na(w))] = 0
  }
  
  ns = ncol(x)
  theta0 = diff(sort(c(runif(ns - 1), 0, 1)))
  
  if (silent == FALSE) {
    message(colnames(bulk.sample)[n])
  }
  
  y = bulk.sample[genes, n]
  fn = function(theta) {
    (w0 * w) %*% (y - log(msf + 1) * (x %*% theta))^2 + sum(theta^2)
  }
  gn = function(theta) {
    -2 * t(x) %*% ((w0 * w) * (y - log(msf + 1) * (x %*% theta))) + 2 * theta
  }
  
  op <- options(show.error.messages = FALSE)
  on.exit(op)
  esti = try(optim(theta0, fn, gr = gn, lower = rep(0, ns),
                   upper = rep(Inf, ns), method = "L-BFGS-B"),
             silent = T)
  if (class(esti) == "try-error") {
    res = rep(NA, length(theta0))
  } else {
    res = esti$par
  }
  
  res = as.matrix(data.frame(thetas = res/sum(res)))
  
  colnames(res) = colnames(bulk.sample)[n]
  rownames(res) = colnames(x)
  options(show.error.messages = T)
  
  ## stop nodes ##
  stopImplicitCluster()
  
  return(res)
}


