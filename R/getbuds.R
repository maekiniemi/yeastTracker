#' Compute area of a polygon
#'
#' This  function computes the area of a polygon X.
#' 
#' @param X a polygon dataframe of x and y coordinates
#' @return area of the polygon x.
#' @export
#' @examples
#' filename<-system.file('data/YET629_02_w1L488nm-L561nm_sequence-10000.zip', package='yeast')
#' roi <- read.ijzip(filename)
#' 
#' cordi <- roi$`001_001`$coords
#' X <-  data.frame(X=cordi[,1], Y=cordi[,2])
#' area<-function(X){
#'  X<-cbind(X$X, X$Y)
#'  X<-rbind(X,X[1,])
#'  x<-X[,1]; y<-X[,2]; lx<-length(x)
#'  -sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2
#' }
#' 
#' area(X)

area<-function(X){
  X<-cbind(X$X, X$Y)
  X<-rbind(X,X[1,])
  x<-X[,1]; y<-X[,2]; lx<-length(x)
  -sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2
}


#' Find intersection of two objects
#'
#' This  function identifies the positions of intersection of two matrix objects.
#' 
#' @param a first matrix object
#' @param b second matrix object
#' @return the coordinates of intersection
#' @export


getIntersection<-function(a, b){
  # the intersection [(x1,y1), (x2, y2)]
  #   it might be a line or a single point. If it is a point,
  #   then x1 = x2 and y1 = y2.  */


  if ( isTRUE(all.equal(a[1,1], a[2,1]) )) {
    # Case (A)
    # As a is a perfect vertical line, it cannot be represented
    # nicely in a mathematical way. But we directly know that
    #
    x1 <- a[1,1]
    x2 <- x1
    if ( isTRUE(all.equal(b[1,1] , b[2,1]) ))  {
      # Case (AA): all x are the same!
      # Normalize
      if(a[1,2] > a[1,2]) {
        tmp<-a[1,] ; a[1,] <- a[2,]; a[2,]<-tmp; #a={"first": a[2,, "second": a[1,};
      }
      if(b[1,2] > b[2,2]) {
        tmp<-b[1,] ; b[1,] <- b[2,]; b[2,]<-tmp; #b = {"first": b[2,, "second": b[1,};
      }
      if(a[1,2] > b[1,2]) {
        tmp<-a
        a<-b
        b<-tmp
      }

      # Now we know that the y-value of a[1, is the
      # lowest of all 4 y values
      # this means, we are either in case (AAA):
      #   a: x--------------x
      #   b:    x---------------x
      # or in case (AAB)
      #   a: x--------------x
      #   b:    x-------x
      # in both cases:
      # get the relavant y intervall
      y1 = b[1,2];
      y2 = min(c(a[2,2], b[2,2]) );
    } else {
      # Case (AB)
      # we can mathematically represent line b as
      #     y = m*x + t <=> t = y - m*x
      # m = (y1-y2)/(x1-x2)
      m = (b[1,2] - b[2,2])/
        (b[1,1] - b[2,1]);
      t = b[1,2] - m*b[1,1];
      y1 = m*x1 + t;
      y2 = y1
    }
  } else{ if( isTRUE(all.equal(b[1,1],b[2,1])  )) {
    # Case (B)
    # essentially the same as Case (AB), but with
    # a and b switched
    x1 = b[1,1];
    x2 = x1;

    tmp = a;
    a = b;
    b = tmp;

    m = (b[1,2] - b[2,2])/
      (b[1,1] - b[2,1]);
    t = b[1,2] - m*b[1,1];
    y1 = m*x1 + t;
    y2 = y1
  } else {
    # Case (C)
    # Both lines can be represented mathematically
    ma = (a[1,2] - a[2,2])/
      (a[1,1] - a[2,1]);
    mb = (b[1,2] - b[2,2])/
      (b[1,1] - b[2,1]);
    ta = a[1,2] - ma*a[1,1];
    tb = b[1,2] - mb*b[1,1];
    if ( isTRUE(all.equal(ma , mb) )) {
      # Case (CA)
      # both lines are in parallel. As we know that they
      # intersect, the intersection could be a line
      # when we rotated this, it would be the same situation
      # as in case (AA)

      # Normalize
      if(a[1,1] > a[2,1]) {
        tmp<-a[1,] ; a[1,] <- a[2,]; a[2,]<-tmp;# a = {"first": a["second"], "second": a["first"]};
      }
      if(b[1,1] > b[2,1]) {
        tmp<-b[1,] ; b[1,] <- b[2,]; b[2,]<-tmp; #b = {"first": b["second"], "second": b["first"]};
      }
      if(a[1,1] > b[1,1]) {
        tmp = a;
        a = b;
        b = tmp;
      }

      # get the relavant x intervall
      x1 = b[1,1];
      x2 = min(a[2,1], b[2,1]);
      y1 = ma*x1+ta;
      y2 = ma*x2+ta;
    } else {
      # Case (CB): only a point as intersection:
      # y = ma*x+ta
      # y = mb*x+tb
      # ma*x + ta = mb*x + tb
      # (ma-mb)*x = tb - ta
      # x = (tb - ta)/(ma-mb)
      x1 = (tb-ta)/(ma-mb);
      y1 = ma*x1+ta;
      x2 = x1;
      y2 = y1;
    }
  }
  }

  if( (x1==x2)&&(y1==y2) ){
    return(c(x1,y1))
  }else{
    return(c(x1,y1,x2,y2))
  }
}



#' getPrincipalAxes
#'
#' This function loads the contour of an object and performs a principal component analysis on it.
#' It returns a list object with principal component 1 and principal component 2  as well as their intersection coordinates.
#' @param contour two dimensional array where the first column constist of x values and the second column constist of y values.
#' @param plot a boolean value specifying whether to plot the ROI including PC1 and PC2. Defauls value is FALSE.
#' @return Returns one matrix for each principal component where the first row represents x,y coordinates for the end of the vector with lowest x value. The second row represents x,y coordinates for the end of the PC2 vector with highest x value. It also returns a matrix with the x,y coordinates for hte PC1 and PC2 intersection, which represents the centroid of the polygon.
#' @export
#' @examples
#' filename<-system.file('data/YET629_02_w1L488nm-L561nm_sequence-10000.zip', package='yeast')
#' roi <- read.ijzip(filename)
#' PCAresult<-getPrincipalAxes(roi$`001_005`$coords, plot=T)

getPrincipalAxes<-function(contour, plot=F){
  load<-princomp(contour)$loadings

  slope <- load[2, ]/load[1, ]
  mn <- apply(contour, 2, mean)
  intcpt <- mn[2] - (slope * mn[1])
  height<-c(min(contour[,2]), max(contour[,2]))
  width<-c(min(contour[,1]), max(contour[,1]))
  y1<-height[1]-diff(height)*0.04
  y2<-height[2]+diff(height)*0.04
  x1<-width[1]-diff(width)*0.04
  x2<-width[2]+diff(width)*0.04
  #first components (usually vertical) get end points out of the contour
  PC_1_x<-( (c(y1, y2)-intcpt[1])/slope[1] )
  #second component (usually horizontal) get end points out of the contour
  PC_2_y <-(intcpt[2]+slope[2]*c(x1, x2) )


  intersect<-getIntersection(matrix(c(x1, x2, PC_2_y[1], PC_2_y[2]), ncol=2), matrix(c(PC_1_x[1], PC_1_x[2], y1, y2), ncol=2))

  if(plot){
    plot(contour, axes=F, type='l', ylab='', xlab='', asp=1)
    polygon(contour, col=gray(0.9))

    #first principal components
    arrows(PC_1_x[1], y1, PC_1_x[2], y2, length=0.15, code=3, lwd=2, col='darkred')
    #second principal components
    arrows(x1, PC_2_y[1], x2, PC_2_y[2], length=0.15, code=3, lwd=2, col='darkred')

    points(intersect[1], intersect[2], pch=21, bg='red', cex=1.5)
  }

  principalaxes<-list( PC2=round(matrix(c(x1, x2, PC_2_y[1], PC_2_y[2]), ncol=2), 3),  PC1= round(matrix(c(PC_1_x[1], PC_1_x[2], y1, y2), ncol=2), 3), intersect=round(intersect,3) )

  return(principalaxes)
}


#' euclidean distance
#' what is does
#' 
#' @param x .roi object from read.ijzip().
#' @return
#' @export
#' @examples


edist<-function(coord){
  sqrt(diff(coord$X)^2 + diff(coord$Y)^2 )
}


#' Import and normalize ImageJ ROIs to cell coordinates
#'
#' This function loads ImageJ roi zip files.
#' It returns a list object with rois normalized into coordinate system as well as mother/bud relations.
#' @param x .roi object from read.ijzip().
#' @param THRESHOLD an integer value specifying the area in pixels needed for an overlap to count as mother/budded relation. Default value is 100.
#' @param borderCol color value for overlap boundary.
#' @param fillcol color value for overlap filling
#' @param add plot in the same coordinate system. default FALSE
#' @param xlab x axis label. Default empty.
#' @param ylab y axis label. Default empty.
#' @param main header. Default "get buds"
#' @param asp set the aspect ratio of the plot. Default 1. 
#' @return
#' @export
#' @examples
#' filename<-system.file('data/YET629_02_w1L488nm-L561nm_sequence-10324.zip', package='yeast')
#' roi <- read.ijzip(filename)
#' yeastCells<-get.buds(roi)



get.buds<-function(roi, THRESHOLD = 1, borderCol = rgb(0,0,0,0.2), fillCol = NA, add=FALSE, xlab = "", ylab = "", main = "get buds", asp = 1, ...) {

  ## Base plot
  if (!add) {
    par(mfrow=c(1,4))
    plot(NA, NA, xlim=range(unlist(lapply(roi, function(i) i$xrange)), na.rm = TRUE), ylim= rev(range(unlist(lapply(roi, function(i) i$yrange)), na.rm = TRUE)) , axes = FALSE, xlab = xlab, ylab = ylab, main = main, asp = asp)
  }


#calculates the mean of each column (the Y positions) of an roi and assign the value to the centroid object. 
# so it basically calculates the centroid of the polygon
  lapply(roi, function(i) {tmp <- i
  class(tmp) <- "ijroi"
  plot(tmp, add = TRUE, ...)
  centroid<-apply(tmp$coord, 2, mean)
  text(centroid[1], centroid[2], labels = tmp$name, cex=0.5)
  })

#goes through all rois and adds the coordinates to a list named p.
#creates an object named cellshape that creates a dataframe of the first and second column of p and assign
#them as X and Y, respectively. These are the X and Y coordinates for the cell outlines.
#it also gives them and id, 1. 
#To identify the position of the centroid of the cells, we then take the X and Y coordinates of the cellshape 
#data frame columns and calculates the mean value. The output is the x and y coordinates of the centroid for
#each cell. 
# To normalize the centroid positions we then take the x and y column coordinates and subtract the centroid value.This
# is to get the cells centered in the same coordinate system. 
# The nomalized values are named X1 and Y1.
# Then a new object called cellshape is created where we combine the columns of cellshape and normalizespos into one
# single dataframe. 
# Then we add another column to the data frame which gives each cell a unique number based on the roi. 
  
  allPolygons <- seq_along(roi)
  cellShape <- list()
  for(i in seq_along(roi)){
    p <- roi[[i]]$coords

    #cellShape[[i]] <- data.frame(X=p[,1], Y=p[,2], id = 1)
    cellShape[[i]] <- data.frame(X=p[,1], Y=p[,2], id = 1)

    centroid <- apply(cellShape[[i]][,1:2], 2, mean)
    normalizedPos <- sweep(cellShape[[i]][,1:2], 2, centroid )
    names(normalizedPos) <- c('X1', 'Y1')
    cellShape[[i]] <- cbind(cellShape[[i]], normalizedPos)
    cellShape[[i]]$roi <- i
    #double brackets subsets a list

    
# q.overlap contains a repeat of 0 for the length of allpolygons, which is 11.
# the for loop then iterates over allpolygons except the current one [-i] and assigns the coordinates to qi. 
# Then it checks what coordinates in qi that overlaps with p.
# Then it will check further if the q.overlap is larger than 0 coordinate points. 
# This is done by taking all the coordinates that passes this requirement and assign them to qi. 
# Then two polygon data frames are created from the coordinates in p and qi respectively.
# PID and POS are from package PBS mapping and assigns the polygon an ID and position. 
    
    q.overlap <- rep(0, length(allPolygons))
    for(q in allPolygons[-i]){
      qi <- roi[[q]]$coords
      q.overlap[allPolygons == q] <- sum( point.in.polygon(qi[,1], qi[,2], p[,1], p[,2]) )
    }
    check.further<- which(q.overlap > 0)
    if(length(check.further)>0){
      for(q in check.further){
        qi <- roi[[q]]$coords

        p1 <- data.frame(PID=rep(1, nrow(p)), POS=1:nrow(p), X=p[,1], Y=p[,2]) 
        p2 <- data.frame(PID=rep(2, nrow(qi)), POS=1:nrow(qi), X=qi[,1], Y=qi[,2])

# Then the polygons are joined and assigned to overlap. If the area of the overlap is larger than the set threshold, 
# default > 100, a new polygon is created from the overlapping coordinates with pink filing and red border.
# Then it checks the absolute value of the x and y coordinates of the p1 and p2 polygons. The largest one will be 
# assigned as the mother and the smaller as bud. Then the two cells are combines to one single polygon and plots the 
# cell, including the overlap polygon. 
        
        overlap <- PBSmapping::joinPolys(p1, p2)
        if( area(overlap) > THRESHOLD ){
          polygon(overlap$X, overlap$Y, col='#FF000020', border='red', lwd=2)
          if( abs(area(p1)) >  abs(area(p2)) ){
            #find mother cell
            #merge into larger polygon
            p1 <- data.frame(X = p1$X, Y = p1$Y, id = 1)#[-which(paste(p1$X, p1$Y) %in% paste(overlap$X, overlap$Y))]
            p2 <- data.frame(X = p2$X, Y = p2$Y, id = 2)#[-which(paste(p2$X, p2$Y) %in% paste(overlap$X, overlap$Y))]

            buddingunion<-rbind(p1, p2)
            buddingunion<-buddingunion[-which(paste(buddingunion$X, buddingunion$Y) %in% paste(overlap$X, overlap$Y)),]
            polygon( buddingunion)
            

            e1<-which.max( edist( buddingunion[buddingunion$id==1,] ) )
            e2<-which.max( edist( buddingunion[buddingunion$id==2,] ) ) + sum(buddingunion$id==1)
            e1<-c(e1, e1+1)
            e2<-c(e2, e2+1)

            distanceBetweenE1<- sqrt( apply(sweep(as.matrix( buddingunion[e2,1:2], ncol=2 ),  2, as.matrix(buddingunion[e1[1],1:2], ncol=2) )^2, 1, sum) )

            indMax<-which.max(distanceBetweenE1)
            indMin<-which.min(distanceBetweenE1)

            #plot(buddingunion[,1:2], type='p')
            #text(buddingunion[,1], buddingunion[,2], 1:nrow(buddingunion))

            index1 <- c(e1[1]:1, sum(buddingunion$id == 1):e1[2] )
            index2<-c(e2[indMax]:min(which(buddingunion$id == 2)), nrow(buddingunion):e2[indMin]  )
            index<-c(index1, index2)
            buddingunion<-buddingunion[index,]

            plot(buddingunion[,1:2], type='n', asp=1)
            polygon(buddingunion[buddingunion$id==1,1:2], col='#ff000020', border=NA)
            polygon(buddingunion[buddingunion$id==2,1:2], col='#0000ff20', border=NA)
            polygon(buddingunion[,1:2], border='black')



            cellShape[[i]] <- buddingunion

# rotate around PCA. 
            
            PC<-getPrincipalAxes(buddingunion[,1:2], plot = T)
            angle <-   atan( PC$PC1[2,1]/PC$PC1[2,2] )
            M <- matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2 )
            buddingunion[,1:2] <- as.matrix(buddingunion[,1:2], ncol=2) %*% M
            


# calculate the new centroid
            centroid <- apply(buddingunion[buddingunion$id==1,1:2], 2, mean)
            normalizedPos <- sweep(buddingunion[,1:2], 2, centroid )
            names(normalizedPos) <- c('X1', 'Y1')
            
            #if bud is negative on Y then reflect Y values so bud always faces up
            centroidBud<-apply(normalizedPos[buddingunion$id==2, ], 2, mean)
            if(centroidBud[2] < 0)
              normalizedPos$Y1 <- (-normalizedPos$Y1)
            
            cellShape[[i]] <- cbind(cellShape[[i]], normalizedPos)
            cellShape[[i]]$roi <- i

            plot(cellShape[[i]]$X1, cellShape[[i]]$Y1, asp=1, type='l', ylab='', xlab='', axes=F)
            polygon(cellShape[[i]]$X1, cellShape[[i]]$Y1)

            print(paste("Polygon", i, "is mother to", q))
          }else{
            cellShape[[i]] <- NA
            print(paste("Polygon", q, "is mother to polygon", i))
          }

        }


      }
    }
  }

  axis(1)
  axis(2, las = 2)



  return(cellShape)

}

