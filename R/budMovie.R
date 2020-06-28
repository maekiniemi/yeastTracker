
library(RImageJROI)
library(yeast)

folder<-system.file('data/results', package='yeast')
timesSeries<-trackCell(folder, THRESHOLD = 1)

files<-dir("/Users/Anna/Desktop/RNAresult/", full.names=T)
files<-files[which(tools::file_ext(files) == 'zip')]


centroid <- data.frame(x = numeric(), y = numeric(), time = integer())
MS2<-list()
for(i in seq_along(files)){
  roi <- read.ijzip(files[i])
  MS2[[i]]<-lapply(roi, function(x){
    X<-smooth.spline(x$coords[,1], df=8)$y;
    Y<-smooth.spline(x$coords[,2], df=8)$y; 
    return(list(x=X, y=Y))})
  for(j in seq_along(roi) ){
    tmp <- apply(roi[[j]]$coords, 2, mean)
    centroid.tmp <- data.frame(x = tmp[1], y = tmp[2], time = i)
    centroid<-rbind(centroid, centroid.tmp )
  }
}

color<-rgb(0,0.7,0,1)


newfolder<-'mymovie'
olddir<-getwd()
dir.create(newfolder)
setwd(newfolder)

sizeBud<-numeric()
for(i in seq(1, length(files)) ){
  if(!all(is.na(timesSeries[[i]][[6]])) ){  
    cell.of.interest<-list(X = timesSeries[[i]][[6]]$X, Y= timesSeries[[i]][[6]]$Y)
  }
  cell.of.interest$X <- smooth.spline(cell.of.interest$X, df=length(cell.of.interest$X)/2 )$y
  cell.of.interest$Y <- smooth.spline(cell.of.interest$Y, df=length(cell.of.interest$Y)/2 )$y
  
  sizeBud<-c(sizeBud, abs(area(cell.of.interest)) )
  
}
plot(sizeBud, type='o', cex=0.5, pch=16, xlim=c(0,600),   )
smoothedBud<-smooth.spline(sizeBud, df=10)
lines(smoothedBud, col='red')
lines(old, col='blue')

plot(sizeBud-smoothedBud$y, type='l')
abline(h=c(300, -400))

remove<-which( (sizeBud-smoothedBud$y) > 300 | (sizeBud-smoothedBud$y) < (-405) )
smoothed<-smooth.spline(sizeBud, df=100)

cell.of.interest<-list(X = numeric(), Y = numeric())
sizeBud<-numeric()

bud<-data.frame()
RNAinBud<-0
for(i in seq(1, length(files)) ){

  if(i == 1){
    cell.of.interest<-list(X = timesSeries[[i]][[6]]$X, Y= timesSeries[[i]][[6]]$Y)
    cell.of.interest$X <- smooth.spline(cell.of.interest$X, df=length(cell.of.interest$X)/2 )$y
    cell.of.interest$Y <- smooth.spline(cell.of.interest$Y, df=length(cell.of.interest$Y)/2 )$y
  }else{
  if(!all(is.na(timesSeries[[i]][[6]])) ){  
    cell.of.interestTMP<-list(X = timesSeries[[i]][[6]]$X, Y= timesSeries[[i]][[6]]$Y)
    if( !(i %in% remove) ){
      cell.of.interest<-cell.of.interestTMP
      cell.of.interest$X <- smooth.spline(cell.of.interest$X, df=length(cell.of.interest$X)/2 )$y
      cell.of.interest$Y <- smooth.spline(cell.of.interest$Y, df=length(cell.of.interest$Y)/2 )$y
      
      if(length(unique(timesSeries[[i]][[6]]$id)) >1){
        bud<-timesSeries[[i]][[6]][ timesSeries[[i]][[6]]$id == 2,1:2]
      }
      
    }
  }
  }

  rna.of.interest<-list(x = centroid[centroid$time == i,1], y = centroid[centroid$time == i,2])
  inside<-point.in.polygon(rna.of.interest$x, rna.of.interest$y, cell.of.interest$X, cell.of.interest$Y) == 1
  rna.of.interest<-data.frame(rna.of.interest)
  
  filename<-paste0('pic', formatC(i, digits = 3, flag='0'), '.png')
  
  png(filename = filename, width=1920, height=1080)
  par(mfrow=c(1,2), mar=c(0,0,0,0))
  plot(rna.of.interest[inside,], pch=21, bg=color, main=i, axes=F, xlab='', ylab='', cex=3, xlim=c(150, 256), ylim=c(280, 410), asp=1, type='n')
  img <- readPNG(images[i])
  rasterImage(img, 118, 274, 118+136, 274+136)
  
  lapply(seq_along(MS2[[i]])[inside], function(x)polygon(MS2[[i]][[x]], border=color, lwd=3))
  
  names(cell.of.interest)<-c("x", "y")
  polygon(cell.of.interest, lwd=6, border=rgb(1,0.65,0,0.5), lty=2) 
  names(cell.of.interest)<-c("X", "Y")
  #sizeBud<-c(sizeBud, abs(area(cell.of.interest)) )
  par(mar=c(4,4,4,8))
  plot(smoothed$x[1:i], smoothed$y[1:i], type='l', cex.axis=4, axes=F, ylab='', xlab='', xlim=c(0, length(files) ), lwd=4, ylim=c(1200, 4000))
  #mtext(side=4, line=3.5, "Cell area", cex=4)
  axis(4, at=c(2200), labels=c("Cell area"), cex.axis=4, lty = 0, col.axis=c('black'))
  axis(4, at=c(3200), labels=c("Clb2 in bud"), cex.axis=4, lty = 0, col.axis=color)
  
  points(i, smoothed$y[i], pch=21, bg="orange", cex=3.5)
  
  if(nrow(bud)>0){
    RNAinBud<- c(RNAinBud, sum( point.in.polygon(rna.of.interest$x, rna.of.interest$y, bud$X, bud$Y) ) )
  }else{
    RNAinBud<- c(RNAinBud, RNAinBud[length(RNAinBud)])
  }
  lines(1:length(RNAinBud), 150*RNAinBud+3000, lwd=4, col=color, type='S')
  points(i, 150*RNAinBud[length(RNAinBud)]+3000, pch=21, bg="green3", cex=3.5)
  
  
  dev.off()
  
}
setwd(olddir)


library(png)

images<-dir("/Users/Anna/Documents/images_MS2_bud", full.names=T)

frameId<-tools::file_path_sans_ext(sapply(strsplit(basename(images), "_t"), "[", 2))
frameId<-as.integer(frameId)
images<-images[order(frameId)]


#Replace the directory and file information with your info
img <- readPNG(dir("/Users/Anna/Documents/images_MS2_bud", full.names=T)[213])

#Set up the plot area
plot(c(0,512), c(0,512), type='n', main="Plotting Over an Image", xlab="x", ylab="y")

#Get the plot information so the image will fill the plot box, and draw it
lim <- par()


plot(do.call("rbind", timesSeries[[125]])[do.call("rbind", timesSeries[[125]])$roi==6,1:2] , type='l')
lapply(roi, function(x){X<-smooth.spline(x$coords[,1], df=8)$y;
Y<-smooth.spline(x$coords[,2], df=8)$y;
polygon(X,Y)})

frame29<-do.call("rbind", yeastMovie[[29]])
frame30<-do.call("rbind", yeastMovie[[30]])

plot(frame29$X, frame29$Y, pch=16, cex=0.2)
plot(frame30$X, frame30$Y, pch=16, cex=0.2)
