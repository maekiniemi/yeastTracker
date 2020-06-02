#' Makes movie out of single cell time stacks
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

install.FFmpeg()

makeMovie <- function(timesSeries){
  
  
}