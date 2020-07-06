#' Makes movie out of single cell time stacks
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with corrected cell IDs if needed.
#' @param folders zip folders containing imageJ rois of yeast cells 
#' @param cll id of cell of interest
#' @return a list object of rois with corrected IDs adjusted to previous time stack, if needed.
#' @export
#' @examples
#' folder<-system.file('data/res', package='yeast')
#' yeastMovie<-trackCell(folder, THRESHOLD = 1)
#' makeMovie(yeastMovie, cll = 6)

#use the example file daniel made "Untitled1". Use the way we plot in the clustering function to save
# as png files. Then use ffmpeg on that folder. 
#ta celler av intresse
#och identifiera på roi ID 
#spara varje polygon plot som png enligt Daniels exempel
#använd ffmpeg enligt daniels exempel i Anteckningar

makeMovie <- function(timesSeries, cll=6){
  
  newfolder<-'mymovie'
  olddir<-getwd()
  dir.create(newfolder)
  setwd(newfolder)
  q<-0
  
for (r in seq_along(timesSeries)){
  cat(paste('Rendering:', q, '\r'))
  if(!is.na(timesSeries[[r]][[cll]])){
    q<-q+1
    x.p <- list()
    y.p <- list()
    x.p <- timesSeries[[r]][[cll]][,1] #4
    y.p <- timesSeries[[r]][[cll]][,2] #5
    filename<-paste0('pic', formatC(q, digits = 3, flag='0'), '.png')
    png(filename = filename, width=1920, height=1080)
    plot(do.call("rbind", timesSeries[[r]])[,1:2], pch=16, cex=1, col=do.call("rbind", timesSeries[[r]])[,6], xlim=c(0,400), asp=1, main=r, ylim=c(0,450))
    #plot(x.p, y.p, type = "n", asp=1, axes=F, xlab='', ylab='')
    polygon(x.p, y.p, lwd=2, lty=3, col='pink')
    x<-tapply(do.call("rbind", timesSeries[[r]])[,1], do.call("rbind", timesSeries[[r]])$roi, mean)
    y<-tapply(do.call("rbind", timesSeries[[r]])[,2], do.call("rbind", timesSeries[[r]])$roi, mean)
    text(x,y, names(x), col=as.integer(names(x)))
    dev.off()
  }
}  
  setwd(olddir)
  
}

system('ffmpeg -r 24 -f image2 -s 1920x1080 -i ./mymovie/pic%04d.png -vcodec libx264 -pix_fmt yuv420p ./RNAmovie.mp4')    
