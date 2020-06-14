#' Imports ImageJ ROIs and adjusts the cell ID if necessary. 
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with corrected cell IDs if needed.
#' @param folder .zip folder containg .roi objects of yeast cells
#' @param full.names a logical value. If TRUE, the directory path is prepended to the file names to give a relative file path. If FALSE, the file names (rather than paths) are returned.Defaul TRUE.
#' @param THRESHOLD an integer value specifying the area in pixels needed for an overlap to count as mother/budded relation. Default value is 100.
#' @return a list object of rois with corrected IDs adjusted to previous time stack, if needed.
#' @export
#' @examples
#' folder<-system.file('data', package='yeast')
#' yeastMovie<-trackCell(folder)

#reads in all roi folders and combine them to one object
trackCell <- function(folder, THRESHOLD = 1){
  zipfiles<-dir(folder, full.names = T)
  zipfiles<-zipfiles[which(tools::file_ext(zipfiles) == 'zip')]
  timesSeries<-list()
  roi <- RImageJROI::read.ijzip(zipfiles[1])
  timesSeries[[1]]<-get.buds(roi, THRESHOLD = THRESHOLD)
  for(i in seq_along(zipfiles)[-1]){
    roi <- RImageJROI::read.ijzip(zipfiles[i])
    
    #saves the rois into a list object
    cellShape.tmp<-get.buds(roi, THRESHOLD = THRESHOLD)
    
    #identifies the centroid of each cell. Why don't we calculate centroid on [[1]]?
    centroid<-lapply(seq_along(cellShape.tmp), function(x){
      if(length(cellShape.tmp[[x]])>1){
        apply(cellShape.tmp[[x]][,1:2], 2, mean)
      }else{
        return(NA)
      }
    })
    
    #checks if the cell centroid with a specific ID is within the polygon of the previous time stack, then they will be assigned 1. If not they are assigned 0.
    position.correct <- lapply(seq_along(timesSeries[[i-1]]), function(x){
      if(!is.na(timesSeries[[i-1]][[x]]) & length(centroid)>= length(timesSeries[[i-1]]) ){
        point.in.polygon(centroid[[x]][1], centroid[[x]][2], timesSeries[[i-1]][[x]]$X,  timesSeries[[i-1]][[x]]$Y)
      }else{
        return(NA)
      }
    })
    
    # takes the cells that were identified as not correct
    # if there are cells that need to change ID it will compute the original centroid of them 
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
      
      # Checks if the original centroid of the cells is found withing the cells of the previous time stack            
      order.change<-integer()
      for(k in seq_along(centroid.original)){
        match.cell <- lapply(seq_along(cellShape.tmp), function(x){
          if(length(cellShape.tmp[[x]])>1){
            point.in.polygon(centroid.original[[k]][1], centroid.original[[k]][2], cellShape.tmp[[x]]$X,  cellShape.tmp[[x]]$Y)
          }else{
            return(NA)
          }
        })
        
        
        # Changes the ID of cells that have the wrong ID   
        order.change <- c(order.change, which(unlist(match.cell) == 1)[1] )
      }
      
      for(l in seq_along(order.change))
        timesSeries[[i]][[has.to.change[l]]] <- cellShape.tmp[[order.change[l]]]
    }
  }
  
  return(timesSeries)
}
