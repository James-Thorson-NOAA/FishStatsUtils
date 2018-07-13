#' @export
plotSites <-
function( mapObj, x, y, col, cex=2, ... ){
  # Load package
  require(maptools)

  # Plot map
  if( is.list(mapObj) & length(mapObj)>0 ){
    plot(mapObj, type="l", xlim=mapObj[["bbox"]]['x',], ylim=mapObj[["bbox"]]['y',], ...)
  }else{
    plot(1, type="n", xlim=attr(mapObj,"bbox")['x',], ylim=attr(mapObj,"bbox")['y',], ...)
    if( length(mapObj)>1 ) suppressWarnings( plot(mapObj, add=TRUE) )
  }
  
  #plot site locations
  points(x=x, y=y, pch=20, cex=cex, col=col)
}
