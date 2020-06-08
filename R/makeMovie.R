#' Makes movie out of single cell time stacks
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with corrected cell IDs if needed.
#' @param folders zip folders containing imageJ rois of yeast cells 
#' @param cll id of cell of interest
#' @return a list object of rois with corrected IDs adjusted to previous time stack, if needed.
#' @export
#' @examples
#' folders<-system.file('data/my.zipfiles', package='yeast')
#' yeastMovie<-trackCell(folder)
#' makeMovie(yeastMovie, cll = 7)

#use the example file daniel made "Untitled1". Use the way we plot in the clustering function to save
# as png files. Then use ffmpeg on that folder. 
#ta celler av intresse
#och identifiera på roi ID 
#spara varje polygon plot som png enligt Daniels exempel
#använd ffmpeg enligt daniels exempel

makeMovie <- function(timesSeries, cll=6){
  
  newfolder<-'mymovie'
  olddir<-getwd()
  dir.create(newfolder)
  setwd(newfolder)
  
for (r in seq_along(timesSeries)){
  x.p <- list()
  y.p <- list()
  x.p <- timesSeries[[r]][[cll]][,4]
  y.p <- timesSeries[[r]][[cll]][,5]
  filename<-paste0('pic', formatC(r, digits = 3, flag='0'), '.png')
  png(filename = filename, width=1920, height=1080)
  plot(x.p, y.p, type = "n", asp=1, axes=F, xlab='', ylab='')
  polygon(x.p, y.p)
  dev.off()

}  
  setwd(olddir)
  
}
  
  
  
  
  
 