

install.packages('Rtsne')
library(Rtsne)

tsne<-Rtsne(dimred)
install.packages('dbscan')
library("dbscan")


set.seed(1040)
tsne<-Rtsne(dimred, perplexity = 30)
cl <- hdbscan(tsne$Y, minPts = 30)

install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

cluster.color <- brewer.pal(length(unique(cl$cluster)), "Paired")


poly<-structure(list(c(-9.99483202391597, -12.788103044381, -10.925922364071,
                       -5.33938032314089, -3.94274481290836, 0.247161717789205, 9.32529253430062,
                       15.610152330347, 15.610152330347, 12.5841087248432, 8.16142960910685,
                       5.83370375871931, 6.29924892879682, 8.16142960910685, 11.1874732146107,
                       8.85974736422311, 0.479934302827959, -4.87383515306338), c(20.5983873209301,
                                                                                  17.1067985453488, 12.218574259535, 13.6152097697675, 18.9689792256589,
                                                                                  25.9521567768215, 28.279882627209, 32.9353343279841, 37.3580134437204,
                                                                                  42.7117828996118, 43.6428732397668, 41.3151473893792, 38.2891037838754,
                                                                                  36.6596956886042, 33.8664246681391, 30.3748358925578, 26.8832471169765,
                                                                                  24.555521266589)), .Names = c("x", "y"))

cl$cluster[which(point.in.polygon(tsne$Y[,1], tsne$Y[,2], poly$x, poly$y) == 1)] <- NA

dev.copy(pdf,'yeastShape_clustering.pdf')
par(mfrow=c(2,3), mar=c(4,4,1,1))
plot(tsne$Y, pch=21, axes=F, ylab='tSNE 2', xlab='tSNE 2', asp=1, cex=1, bg=cluster.color[cl$cluster+1] )
pos<-apply(tsne$Y[cl$cluster == 5,], 2, mean)
text(pos[1]+3, pos[2], 'cluster 5', pos = 4, col = cluster.color[5+1])
pos<-apply(tsne$Y[cl$cluster == 6,], 2, mean)
text(pos[1], pos[2]+5, 'cluster 6', pos = 3, col = cluster.color[6+1])

par(mar=c(7,0,0,0))
for(k in c(2,3,5,6,7) ){
  plot(0, 0, col.sub = cluster.color[k+1], font.sub = 4,  sub=paste('Cluster:', k), main=''  , type='n', axes=F, xlab='', ylab='', xlim=range(dimred), ylim=range(dimred), asp=1)
  xAvg<-numeric()
  yAvg<-numeric()
  for(l in which(cl$cluster == k)){
    x <- dimred[l,1:100]
    y <- dimred[l,101:200]
    polygon(x , y, border=rgb(0,0,0, 0.1) )
    xAvg<-cbind(xAvg, x)
    yAvg<-cbind(yAvg, y)
  }
  xAvg<-apply(xAvg, 1, mean)
  yAvg<-apply(yAvg, 1, mean)
  polygon(xAvg , yAvg, border=rgb(0,0,0, 1), lwd=3 )
  polygon(xAvg , yAvg, border=rgb(1,0,0, 1), lwd=2 )
}
dev.off()


dimred<-1:200

for(j in which( !is.na(masterYeasts) )){
  contour <- masterYeasts[[j]][,4:5]
  contour<-cbind(contour$X1,contour$Y1)
  contour<-rbind(contour, contour[1,])



  #code for giving each cell equal amount of points describing the contour
  river<-SpatialLines(list(Lines(Line(contour), ID="a")))
  river<-spsample(river, n = 100, type = "regular")

  river<-data.frame(river)


  center<-which.min( sqrt( (river[,2] - max(river[,2]) )^2 + (river[,1] - 0)^2 ) )


  river <- river[c( (center+1):nrow(river), 1:center),]
  dimred<-rbind(dimred, c(river[,1], river[,2] ) )


}

#CLUSTERING


dimred<-dimred[-1,]
