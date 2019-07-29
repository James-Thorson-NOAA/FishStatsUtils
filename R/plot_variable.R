#' @title
#' Plot maps with areal results
#'
#' @description
#' \code{PlotMap_Fn} is a hidden function to plot a map and fill in regions with colors to represent intensity in an areal-interpretion of model results
#'
#' @inheritParams FishStatsUtils::plot_maps
#' @inheritParams sp::CRS
#' @inheritParams par
#'
#' @param Y_gt matrix where values for every column are plotted as a map
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param legend_x two numeric values (between 0 and 1) giving left and right-hand location of color legend
#' @param legend_y two numeric values (between 0 and 1) giving bottom and top location of color legend
#' @param map_list output from \code{FishStatsUtils::make_map_info}
#' @param zlim range for defining bounds of color scale
#' @param add boolean indicating whether to add plot to an existing panel figure, or to define a new panel figure
#' @param outermargintext vector defining text to plot in outer margins of panel figure
#' @param panel_labels vector defining titles to use for each panel; defaults to blank
#'
#' @param ... arguments passed to \code{par}
#'
#' @details
#' This function was necessary to build because \code{mapproj::mapproject} as used in \code{maps::map} has difficulties with both rotations (for most projections) and
#' truncating the cocuntry boundaries within the plotting region (which \code{mapproj::mapproject} appears to do prior to projection,
#' so that the post-projection is often missing boundaries that are within the plotting rectangle).  I use rectangular projections by default, but Lamberts or Albers conformal
#' projections would also be useful for many cases.

#' @export
plot_variable <-
function( Y_gt, map_list, panel_labels, projargs='+proj=longlat', map_resolution="medium",
         file_name="density", working_dir=paste0(getwd(),"/"), Format="png", Res=200, textmargin="", add=FALSE,
         outermargintext=c("Eastings","Northings"), zlim, col, mar=c(0,0,2,0), oma=c(4,4,0,0),
         legend_x=c(0,0.05), legend_y=c(0.05,0.45), mfrow, land_color="grey",
         ...){

  ###################
  # Settings and inputs
  ###################

  # Check for problems and fill in missing stuff
  if( is.vector(Y_gt)){
    Y_gt = matrix(Y_gt, ncol=1)
  }
  if( missing(zlim)){
    zlim = range(Y_gt, na.rm=TRUE)
  }
  if( missing(map_list) || is.null(map_list$MapSizeRatio) ){
    MapSizeRatio = c(3, 3)
  }else{
    MapSizeRatio = map_list$MapSizeRatio
  }
  Y_gt = Y_gt[ map_list$PlotDF[which(map_list$PlotDF[,'Include']>0),'x2i'], , drop=FALSE]
  if( missing(mfrow) ){
    mfrow = ceiling(sqrt(ncol(Y_gt)))
    mfrow = c( mfrow, ceiling(ncol(Y_gt)/mfrow) )
  }
  if( missing(panel_labels) ){
    panel_labels = rep("", ncol(Y_gt))
  }
  if( length(panel_labels) != ncol(Y_gt) ){
    warning( "panel_labels and `ncol(Y_gt)` don't match: Changing panel_labels'")
    panel_labels = 1:ncol(Y_gt)
  }
  if( missing(col)){
    col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
  }
  if( is.function(col)){
    col = col(1000)
  }
  if( all(is.numeric(c(legend_x,legend_y))) ){
    if( any(c(legend_x,legend_y)>1) | any(c(legend_x,legend_y)<0) ){
      stop("Check values for `legend_x` and `legend_y`")
    }
  }
  # Location of extrapolation-grid cells
  loc_g = map_list$PlotDF[which(map_list$PlotDF[,'Include']>0),c('Lon','Lat')]

  # CRS for original and new projections
  CRS_orig = sp::CRS( '+proj=longlat' )
  CRS_proj = sp::CRS( projargs )

  # Data for mapping
  map_data = rnaturalearth::ne_countries(scale=switch(map_resolution, "low"=110, "medium"=50, "high"=10, 50 ))
  map_data = sp::spTransform(map_data, CRSobj=CRS_proj)

  ###################
  # Make panel figure
  ###################

  # Define device
  Par = list( mfrow=mfrow, mar=mar, oma=oma, ...)
  if(Format=="png"){
    png(file=paste0(working_dir,file_name,".png"),
        width=Par$mfrow[2]*MapSizeRatio[2],
        height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
    on.exit( dev.off() )
  }
  if(Format=="jpg"){
    jpeg(file=paste0(working_dir,file_name,".jpg"),
         width=Par$mfrow[2]*MapSizeRatio[2],
         height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
    on.exit( dev.off() )
  }
  if(Format%in%c("tif","tiff")){
    tiff(file=paste0(working_dir,file_name,".tif"),
         width=Par$mfrow[2]*MapSizeRatio[2],
         height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
    on.exit( dev.off() )
  }
  if(add==FALSE) par( Par )

  # Loop across columns (years)
  for( tI in 1:ncol(Y_gt) ){
    # Read extrapolation grid
    Points_orig = sp::SpatialPointsDataFrame( coords=loc_g, data=data.frame("y"=Y_gt[,tI]), proj4string=CRS_orig )

    # Reproject to Lat-Long
    Points_LongLat = sp::spTransform( Points_orig, sp::CRS('+proj=longlat') )

    # Re-project to plotting CRS
    Points_proj = sp::spTransform( Points_orig, CRS_proj )

    # Interpolate to raster
    # library(plotKML)
    cell.size = mean(diff(Points_proj@bbox[1,]),diff(Points_proj@bbox[2,]))/floor(sqrt(nrow(Y_gt)))
    Raster_proj = plotKML::vect2rast( Points_proj, cell.size=cell.size )
    image( Raster_proj, col=col, zlim=zlim )

    # Plot maps using rnaturalearth
    sp::plot( map_data, col=land_color, add=TRUE )

    # Title and box
    title( panel_labels[tI], line=0.1, cex.main=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main), cex=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main) )
    box()
  }

  # Include legend
  if( all(is.numeric(c(legend_x,legend_y))) ){
    xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
    xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
    yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
    yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
    plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(zlim[1],zlim[2],length=4),1), rect.col=col, cex=1, align=c("lt","rb")[2], gradient="y")
  }

  # Margin text
  if(add==FALSE) mtext(side=1, outer=TRUE, outermargintext[1], cex=1.75, line=par()$oma[1]/2)
  if(add==FALSE) mtext(side=2, outer=TRUE, outermargintext[2], cex=1.75, line=par()$oma[2]/2)

  # return stuff as necessary
  return( invisible(list("Par"=Par,"cell.size"=cell.size)) )
}
