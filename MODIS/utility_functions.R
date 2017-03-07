# Some R utility fnctions to trim a raster and display a legend on an image

# Trim a raster to a defined bounding box (code taken from the web)
trim2 <- function(x,values=NA,out="matrix"){
  if(!any(out==c("matrix","raster"))) stop("output must be a matrix or raster")
  if(class(x)=="matrix" & out=="raster") stop("if you supply a matrix, you must use out='matrix'")
  if(class(x)=="RasterLayer") {
    if(out=="raster") { cres <- 0.5*res(x); crs <- projection(x); y <- x }
    x <- matrix(as.array(x),nrow=nrow(x),ncol=ncol(x))
  }
  if(class(x)!="matrix") { stop("x must be a matrix or raster")
  } else {
    r.na <- c.na <- c()
    if (is.na(values)) {
      for(i in 1:nrow(x)) r.na <- c(r.na, all(is.na(x[i,])))
      for(i in 1:ncol(x)) c.na <- c(c.na, all(is.na(x[,i])))
    } else {
      for(i in 1:nrow(x)) r.na <- c(r.na, all(x[i,]==values))
      for(i in 1:ncol(x)) c.na <- c(c.na, all(x[,i]==values))
    }
    r1 <- 1 + which(diff(which(r.na))>1)[1]; r2 <- nrow(x) -  which(diff(which(rev(r.na)))>1)[1]
    c1 <- 1 + which(diff(which(c.na))>1)[1]; c2 <- ncol(x) - which(diff(which(rev(c.na)))>1)[1]
    x <- x[r1:r2,c1:c2]
    if(out=="raster") {
      xs <- xFromCol(y,col=c(c1,c2)) + c(-1,1)*cres[1]
      ys <- yFromRow(y,row=c(r2,r1)) + c(-1,1)*cres[2]
      x <- crop(y,extent(xs,ys))
    }
  }
  return(x)
}

##############################################################################
# Display a legend on an image

image.legend <-
  function(x,y, zlim, at.z = NULL, col = heat.colors(12), legnd=NULL,
           lwd = max(3,32/length(col)), bg = NA, bty = "", ...)
    ## * kein y.i -- Benutzer soll rein ueber lwd steuern; sollte reichen.
    ## * legnd koennte interessant sein, falls Text geschrieben werden soll
    ##   (weiss mal wieder nicht, wie man aus legnd legend als option
    ##     macht)
    ## * lwd wird per default in Abh. von col gewaehlt.
  {
    ## Purpose:
    ## Authors: Martin Maechler,   9 Jul 2001
    ##          Martin Schlather, 24 Jul 2001
    
    if (!is.null(legnd) && is.null(at.z))
      stop("at.z must be given if legnd is") ## falls legnd darf at.z
    ##                                nicht automatisch gewaehlt werden
    
    if(!is.numeric(zlim) || zlim[1] > zlim[2])
      stop("`zlim' must be numeric; zlim[1] <= zlim[2]")
    if(is.null(at.z)) {
      ## hier ein Versuch in Abhaengigkeit von n
      ## die Anzahl der labels zu bestimmen:
      n <- min(5, max(1,length(col)/10))
      at.z <- pretty(zlim,n=n,min.n=max(n %/% 3,1))
      
      ## es sieht nicht schoen aus, wenn pretty die letzte oder
      ## erste zahl weit ausserhalb des zlim legt.
      ## heuristisch nur 25%  (oder so) ueberschreitung bzgl
      ## intervalllaenge zulassen:
      tol <- diff(at.z)[1] / 4
      at.z <- at.z[(at.z>=zlim[1]-tol) & (at.z<=zlim[2]+tol)]
    }
    if(!is.numeric(at.z) || is.unsorted(at.z))
      stop("`at.z' must be numeric non-decreasing")
    n.at <- length(at.z)
    nc   <- length(col)
    if(n.at >= nc)
      stop("length(at.z) must be (much) smaller than length(col)")
    dz <- diff(zlim)
    ## The colors must run equidistantly from zlim[1] to zlim[2];
    ## col[i] is for z-interval zlim[1] + [i-1, i) * dz/nc  ; i = 1:nc
    ## i.e., an at.z[] value z0 is color i0 = floor(nc * (z0 - zlim[1])/dz)
    at.i <- floor(nc * (at.z - zlim[1])/dz )
    ## Possibly extend colors by `background' to the left and right
    bgC <- if(is.null(bg)) NA else bg
    if((xtra.l <- 1 - at.i[1]) > 0) {
      at.i <- at.i + xtra.l
      col <- c(rep(bgC, xtra.l), col)
    }
    if((xtra.r <- at.i[n.at] - nc) > 0)
      col <- c(col, rep(bgC, xtra.r))
    lgd <- character(length(col))
    
    ## folgende if-Anweisung ist neu:
    if (is.null(legnd)) lgd[at.i] <-format(at.z, dig = 3)
    else {
      if (length(legnd)!=length(at.z))
        stop("at.z and legnd must have the same length")
      lgd[at.i] <- legnd
    }
    if((V <- R.version)$major <= 1 && V$minor <= 3.0 && V$status == "")
    {
      ## stop-gap fix around the bug that "NA" is not a valid color:
      if(is.na(bgC)) {
        lgd <- lgd[!is.na(col)]
        col <- col[!is.na(col)]
      }
    }
    legend(x,y, legend = rev(lgd), col = rev(col),
           y.i = lwd/16, bty = bty, lwd = lwd, bg = bg, ...)
  }
