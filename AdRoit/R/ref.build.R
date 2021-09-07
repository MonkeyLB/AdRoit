#' Build the single cell reference
#'
#' This function takes a vector of gene counts, fits the negative binomial distribution, and estimate the size and mean parameters.
#'
#' @param counts a matrix/data.frame/sparse matrix (dgCMatrix) of single cell counts with the rows being genes and columns being cells.
#' @param annotations a vector annotates the cell types of cells.
#' @param genes the set of genes selected and to be used for deconvolution.
#' @param samples a vector annotates which sample ID the cells come from. Default is `NULL`.
#' @param normalize specify how to normalize the single cell raw counts. It should be chosen from "None", "Total" and "Median". Default is "None".
#' @param scale.factor the total number of reads to be normalized to if `normalize = "Total"`. Default is 1e+5.
#' @param multi.sample.bulk specify whether the bulk to be deconvoluted has multiple samples. Default is TRUE.
#' @param multi.sample.single specify whether the single cell reference has multiple samples. Default is TRUE.
#' @param nbootsids specify the number of samples simulated via bootstapping using the single cell reference data. Default is 5.
#' @param minbootsize specify the minimum number of cells that are used to synthesize one bulk sample for each cell type and each sample. Default is 50.
#' @param silent whether to print out messages. Default is FALSE.
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @return a list of elements including `mus`, `sizes`, `cell.specificity.w` and `cross.sample.w`.
#' \itemize{
#'   \item{mus    }{estimated means for selected genes.}
#'    \item{sizes    }{estimated size parameters for selected genes.}
#'    \item{cell.specificity.w    }{cell-type specificity weights.}
#'    \item{cross.sample.w    }{cross-sample varaibility weights.}
#'    }
#' @export

ref.build <- function(counts,
                       annotations,
                       genes,
                       samples = NULL,
                       normalize = "None",
                       scale.factor = 1e+05,
                       multi.sample.bulk = TRUE,
                       multi.sample.single = TRUE,
                       nbootsids=5, 
                       minbootsize=50,
                       silent = FALSE){
    
    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)
    negbin.par <- list()
    counts=round(counts)
    cell.types = unique(sort(as.vector(annotations)))
    for (c in cell.types) {
      cells.id = colnames(counts)[which(annotations == c)]
      if (normalize == "Total") {
        temp = round(scale.factor * (sweep(counts[, cells.id], 
                                           2, colSums(counts[, cells.id]), `/`)))
      } else if (normalize == "Median") {
        medcnt = median(colSums(counts[, cells.id]))
        temp = round(medcnt * (sweep(counts[, cells.id], 
                                     2, colSums(counts[, cells.id]), `/`)))
      } else if (normalize == "None") {
        temp = counts[, cells.id]
      } else {
        message("normalize has to be one of Total, Median and None")
      }
      
      if (silent == FALSE) {
        message(paste("Estimate means and dispersions for cell type", 
                      c, sep = ": "))
      }
      dim(temp)
      negbin.par[[c]] <- foreach(i = genes, .combine = rbind) %dopar% 
        negbin.est(temp[i, ])
      rownames(negbin.par[[c]]) = genes
      rm(temp)
    }
    mus <- sizes <- NULL
    for (i in 1:length(negbin.par)) {
      mus = cbind(mus, negbin.par[[i]][, 2])
      sizes = cbind(sizes, negbin.par[[i]][, 1])
    }
    vars = mus + mus^2/sizes
    colnames(mus) = colnames(sizes) = cell.types
    midx = apply(mus, 1, which.max)
    w0 <- NULL
    for (i in 1:length(midx)) {
      w0 = c(w0, mus[i, midx[i]]/vars[i, midx[i]])
    }
    w0[which(is.na(w0) | is.infinite(w0))] = 0
    if (multi.sample.bulk == TRUE) {
      ref.est = list(mus = mus, lambda = sizes, cell.specificity.w = w0)
    } else {
      bootcounts=counts
      bootannots=as.vector(annotations)
      if (multi.sample.single == TRUE) {
        if (is.null(samples)) {
          stop("Please provide sample ID for each cell.")
        }
        bootsids = samples
      } else {
        bootcounts=bootannots=bootsids=c()
        for (ct in cell.types){
          selcells= which(annotations==ct)
          if(length(selcells)>=nbootsids*minbootsize){ 
            tmpbootcounts=counts[,selcells,drop=F]
            tmpbootannots=as.vector(annotations[selcells,drop=F])
            tmpbootsids=paste0("sample_",mapper[as.character(sample(1+((1:length(selcells))-1) %% nbootsids))])
          }else{
            tmpbootcounts=tmpbootannots=tmpbootsids=c()
            nrounds=ceiling(nbootsids*minbootsize/length(selcells))
            
            for (rd in 1:nrounds){
              mapper=1:nbootsids; names(mapper)=sample(1:nbootsids)
              addcounts=counts[,selcells,drop=F]
              addannots=as.vector(annotations[selcells,drop=F])
              addsids=paste0("sample_",mapper[as.character(sample(1+((1:length(selcells))-1) %% nbootsids))])
              colnames(addcounts)=names(addannots)=names(addsids)=paste0(colnames(counts[,selcells,drop=F]),"_",rd)
              if(rd<nrounds){
                tmpbootcounts=cbind(tmpbootcounts, addcounts)
                tmpbootannots=c(tmpbootannots,addannots)
                tmpbootsids=c(tmpbootsids, addsids)
              }else{
                rest=minbootsize-table(tmpbootsids)
                whichrest=setdiff(unlist(sapply(names(rest),function(x){as.character(names(addsids)[which(addsids==x)[1:max(0,rest[x])]])})),NA)
                tmpbootcounts=cbind(tmpbootcounts, addcounts[,whichrest,drop=F])
                tmpbootannots=c(tmpbootannots,addannots[whichrest,drop=F])
                tmpbootsids=c(tmpbootsids, addsids[whichrest,drop=F])
              }
            }
          }
          bootcounts=cbind(bootcounts, tmpbootcounts)
          bootannots=c(bootannots, tmpbootannots)
          bootsids=c(bootsids, tmpbootsids)
        }
      }
      w.mu=w.sigma=c()
      single.pool=synthesize.celltype(bootcounts, bootannots, SampleID = bootsids)
      for(eachcelltype in names(single.pool)){
        smv = foreach(i = genes, .combine = rbind) %dopar% negbin.est(single.pool[[eachcelltype]][["simbulk"]][i, 
                                                                                                               ])
        rownames(smv) = genes
        w.mu=cbind(w.mu, smv[, 1])
        w.sigma=cbind(w.sigma, smv[, 2])
      }
      colnames(w.mu)=names(single.pool)
      colnames(w.sigma)=names(single.pool)
      ref.est = list(mus = mus, lambda = sizes, cell.specificity.w = w0, 
                     cross.sample.w = list("w.mu"=w.mu,"w.sigma"=w.sigma))
    }
    stopImplicitCluster()
    return(ref.est)
}
