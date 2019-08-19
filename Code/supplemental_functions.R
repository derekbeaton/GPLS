shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {

    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}

svd.line <- function(X,...){
  svd.res <- tolerance.svd(X)
  v1 <- svd.res$d[1] * svd.res$v[,1]
	v2 <- svd.res$d[2] * svd.res$v[,2]
	segments(-v1[1],-v1[2],v1[1],v1[2],...)
}

sp.latentvar_plot <- function(res, axis=1, pch=20, col="mediumorchid4", line.col="grey80", lty=2, lwd=2,
                              main="Latent Variables",
                              pos=3, display_names=T,
                              ...){
  
  these.points <- cbind(res$lx[,axis],res$ly[,axis])
  rownames(these.points) <- rownames(res$lx)
  plot(these.points, type="n", main=main, xlab=paste0("LX ",axis), ylab=paste0("LY ",axis), axes=F)
  abline(h=0,lty=2,lwd=2, col="grey60")
  abline(v=0,lty=2,lwd=2, col="grey60")
  points(these.points, col=col, pch=pch, ...)
  ## will try to employ a "repel" later.
  if(display_names){
    text(these.points,labels=rownames(these.points),pos=pos,col=col)
  }
  
  
}