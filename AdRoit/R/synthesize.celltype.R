#' Synthesize celltype by pooling single cells per cell type per sample
#'
#' @param counts single cell UMI count matrix or sparse matrix (dgCMatrix), with rownames being genes and colnames being cell names.
#' @param annotations cell type annotation, a vector of length equal to the number of columns in counts.
#' @param SampleID sample identifications, a vector of length equal to the number of columns in counts.
#'
#' @return simulated cell type specific profile per sample
#' @export

synthesize.celltype <- function(counts, annotations, SampleID){
  
  clusters = unique(sort(annotations))
  s.SCcounts = counts
  s.annotations = annotations
  s.SampleID = SampleID
  rs=list()
  for(each in clusters){
    simbulk <- NULL
    selcells=which(s.annotations==each)
    for (i in unique(s.SampleID[selcells])){
      idx = which(s.SampleID[selcells] == i)
      if(length(idx)==1){
        simbulk <- cbind(simbulk, (s.SCcounts[,selcells[idx]]))
      }else{
        simbulk <- cbind(simbulk, rowSums(s.SCcounts[,selcells[idx]]))
      }
    }
    count.table = table(s.SampleID[selcells], s.annotations[selcells]) %>% as.data.frame.matrix() %>%
      t()
    rownames(simbulk) = rownames(counts)
    colnames(simbulk) = unique(s.SampleID[selcells])
    rs[[each]]=list("simbulk"=simbulk, "count"=count.table)
  }
  
  return(rs)
}

