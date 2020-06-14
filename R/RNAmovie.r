files<-dir('E:/BudFinder/STARDIST/cell_outline_train/RNA_3D/RNApredictResult-20200614T162351Z-001/RNApredictResult', full.names=T)
files<-files[which(tools::file_ext(files) == 'zip')]

library(RImageJROI)


centroid <- data.frame(x = numeric(), y = numeric(), time = integer())
for(i in seq_along(files)){
  roi <- read.ijzip(files[i])
  for(j in seq_along(roi) ){
    tmp <- apply(roi[[j]]$coords, 2, mean)
    centroid.tmp <- data.frame(x = tmp[1], y = tmp[2], time = i)
    centroid<-rbind(centroid, centroid.tmp )
  }
}

plot(centroid[,1], centroid[,2], pch=16, cex=0.25, xlim=c(0, 512), ylim=c(512,0 ), asp=1)
Sys.sleep(time)

