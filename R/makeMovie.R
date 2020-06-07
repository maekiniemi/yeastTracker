#' Makes movie out of single cell time stacks
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with corrected cell IDs if needed.
#' @param timesSeries 
#' @return a list object of rois with corrected IDs adjusted to previous time stack, if needed.
#' @export
#' @examples
#' folder<-system.file('data/my.zipfiles', package='yeast')
#' yeastMovie<-trackCell(folder)

#use the example file daniel made "Untitled1". Use the way we plot in the clustering function to save
# as png files. Then use ffmpeg on that folder. 

makeMovie <- function(timesSeries){
  
  newfolder<-'mymovie'
  dir.create(newfolder)
  setwd(newfolder)
  
  for (i in (timesSeries[[i]])){
    
xycoord<-list()

    
    xAvg<-numeric()
    yAvg<-numeric()
    for(l in which(clusters == k)){
      x <- dimred[l,1:100]
      y <- dimred[l,101:200]
      polygon(x , y, border=rgb(0,0,0, 0.1) )
      xAvg<-cbind(xAvg, x)
      yAvg<-cbind(yAvg, y)
    }
    
  }
  
  for(q in 1:nrow(data)){
    filename<-paste0('pic', formatC(q, digits = 3, flag='0'), '.png')
    png(filename = filename, width=1920, height=1080)
    plot(c(-1, 1), c(-1, 1), type = "n", asp=1, axes=F, xlab='', ylab='')
    points(data[q,], pch=16, cex=3)
    dev.off()
  }
  
  system('ffmpeg -r 24 -f image2 -s 1920x1080 -i ./pic%04d.png -vcodec libx264 -pix_fmt yuv420p ../output.mp4')
  
}