#' Imports ImageJ ROIs and control cell ID
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with corrected cell IDs if needed.
#' @param folder .zip folder containg .roi objects
#' @param full.names a logical value.  If TRUE, the directory path is prepended to the file names to give a relative file path. If FALSE, the file names (rather than paths) are returned.Defaul TRUE.
#' @return a list object of rois with corrected IDs adjusted to previous time stack, if needed.
#' @export
#' @examples
#' folder<-'E:/BudFinder/STARDIST/cell_outline_train/Cells_2D/results/results'
#' yeastMovie<-trackCell(folder)

trackCell <- function(folder){
  zipfiles<-dir(folder, full.names = T)
  zipfiles<-zipfiles[which(tools::file_ext(zipfiles) == 'zip')]
  timesSeries<-list()
  roi <- read.ijzip(zipfiles[1])
  timesSeries[[1]]<-get.buds(roi)
  for(i in seq_along(zipfiles)[-1]){
    roi <- read.ijzip(zipfiles[i])
    
    
    cellShape.tmp<-get.buds(roi)
    

    centroid<-lapply(seq_along(cellShape.tmp), function(x){
      if(length(cellShape.tmp[[x]])>1){
        apply(cellShape.tmp[[x]][,1:2], 2, mean)
      }else{
        return(NA)
      }
    })
    
    
    position.correct <- lapply(seq_along(timesSeries[[i-1]]), function(x){
      if(length(timesSeries[[i-1]][[x]])>1){
        point.in.polygon(centroid[[x]][1], centroid[[x]][2], timesSeries[[i-1]][[x]]$X,  timesSeries[[i-1]][[x]]$Y)
      }else{
        return(NA)
      }
    })
    
    has.to.change<-which(position.correct == 0)
    timesSeries[[i]] <- cellShape.tmp
    if(length(has.to.change)>0){
  
      centroid.original<-lapply(has.to.change, function(x){
        if(length(timesSeries[[i-1]][[x]])>1){
          apply(timesSeries[[i-1]][[x]][,1:2], 2, mean)
        }else{
          return(NA)
        }
      })
      
      order.change<-integer()
      for(k in seq_along(centroid.original)){
      match.cell <- lapply(seq_along(cellShape.tmp), function(x){
        if(length(cellShape.tmp[[x]])>1){
          point.in.polygon(centroid.original[[k]][1], centroid.original[[k]][2], cellShape.tmp[[x]]$X,  cellShape.tmp[[x]]$Y)
        }else{
          return(NA)
        }
      })
    
      order.change <- c(order.change, which(unlist(match.cell) == 1) )
      }
    
      for(l in seq_along(order.change))
        timesSeries[[i]][[has.to.change[l]]] <- cellShape.tmp[[order.change[l]]]
    }
  }
 
  return(timesSeries)
}
