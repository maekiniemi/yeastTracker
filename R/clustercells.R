#' Cluster cells based on shape
#'
#' This function loads the result of the get.buds function and 
#' It returns a list object with rois normalized into coordinate system as well as mother/bud relations.
#' @param yeastCells list object with the output from the get.buds function
#' @param K number of clusters. Default value = 3. 
#' @return
#' @export
#' @examples
#' filename<-system.file('data/YET629_02_w1L488nm-L561nm_sequence-10000.zip', package='yeast')
#' roi <- read.ijzip(filename)
#' yeastCells<-get.buds(roi)
#' cluster.cells(yeastCells)

cluster.cells<-function(yeastCells, K = 3){
  dimred<-1:200

#takes the normalized x and y coordinates and assigns it to an data frame object called contour. 
#then it combines the x and y coordinates into a matrix and closes it.

  for(j in which( !is.na(yeastCells) )){
    contour <- yeastCells[[j]][,4:5]
    contour<-cbind(contour$X1,contour$Y1)
    contour<-rbind(contour, contour[1,])


#code for giving each cell equal amount of points describing the contour
# is makes contour into an object of class SpatialLines, called river
# then it sample point locations of the river object and makes 100 points with x and y coordinates for the polygon/cell so that all cells have the same number of points
# then it makes the river object into a data frame  
    
# Identifies the center coordinates of the cell/mother.      
     river<-SpatialLines(list(Lines(Line(contour), ID="a")))
    river<-spsample(river, n = 100, type = "regular")

    river<-data.frame(river)


    center<-which.min( sqrt( (river[,2] - max(river[,2]) )^2 + (river[,1] - 0)^2 ) )

#Skriver om ordningen på coordinaterna/punkterna i polygonen
    river <- river[c( (center+1):nrow(river), 1:center),]
    dimred<-rbind(dimred, c(river[,1], river[,2] ) )


  }

  #CLUSTERING

#removes first row
  dimred<-dimred[-1,]

# dist computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix
  dist <- dist(dimred , diag=TRUE)

  # Hierarchical Clustering with hclust
  hc <- hclust(dist)

  clusters<-cutree(hc, k = K)
  par(mfrow=c(1, max(clusters, na.rm=T)+1  ))
  # Plot the result
  plot(hc, main='clustering')

  
#   
  for(k in unique(na.omit(clusters) ) ){
    plot(0, 0, main=paste('Cluster:', k), type='n', axes=F, xlab='', ylab='', xlim=range(dimred), ylim=range(dimred), asp=1)
    xAvg<-numeric()
    yAvg<-numeric()
    for(l in which(clusters == k)){
      x <- dimred[l,1:100]
      y <- dimred[l,101:200]
      polygon(x , y, border=rgb(0,0,0, 0.1) )
      xAvg<-cbind(xAvg, x)
      yAvg<-cbind(yAvg, y)
    }
    xAvg<-apply(xAvg, 1, mean)
    yAvg<-apply(yAvg, 1, mean)
    polygon(xAvg , yAvg, border=rgb(1,0,0, 1), lwd=1 )
  }

  return(clusters)

}

