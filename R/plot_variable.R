#' @title
#' Plot maps with areal results
#'
#' @description
#' \code{PlotMap_Fn} is a hidden function to plot a map and fill in regions with colors to represent intensity in an areal-interpretion of model results
#'
#' @inheritParams plot_maps
#' @param plot_legend_fig Boolean, whether to plot a separate figure for the heatmap legend or not
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param ... arguments passed to \code{par}
#'
#' @details
#' This function was necessary to build because \code{mapproj::mapproject} as used in \code{maps::map} has difficulties with both rotations (for most projections) and
#' truncating the cocuntry boundaries within the plotting region (which \code{mapproj::mapproject} appears to do prior to projection,
#' so that the post-projection is often missing boundaries that are within the plotting rectangle).  I use rectangular projections by default, but Lamberts or Albers conformal
#' projections would also be useful for many cases.

#' @export
plot_variable <-
function( Y_gt, spatial_list, extrapolation_list, map_info, panel_labels, proj4string='+proj=longlat', map_resolution="medium",
         file_name="density", working_dir=paste0(getwd(),"/"), Format="png", Res=200, textmargin="", add=FALSE, pch=15,
         outermargintext=c("Eastings","Northings"), zlim, col, mar=c(0,0,2,0), oma=c(4,4,0,0),
         Legend=list("use"=FALSE, "x"=c(10,30), "y"=c(10,30)), mfrow, plot_legend_fig=TRUE, land_color="grey",
         ...){

  # Check for problems and fill in missing stuff
  if( is.vector(Y_gt)){
    Y_gt = matrix(Y_gt, ncol=1)
  }
  if( missing(zlim)){
    zlim = range(Y_gt, na.rm=TRUE)
  }
  if( missing(map_info) ){
    MapSizeRatio = c(3, 3)
  }else{
    MapSizeRatio = map_info$MapSizeRatio
  }
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

  # Plotting functions
  if( missing(col)) col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
  if( is.function(col)) col = col(1000)

  # Plot
  Par = list( mfrow=mfrow, mar=mar, oma=oma, ...)
  if(Format=="png"){
    png(file=file.path(working_dir, paste0(file_name,".png")),
        width=Par$mfrow[2]*MapSizeRatio[2],
        height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
  }
  if(Format=="jpg"){
    jpeg(file=paste0(FileName, paste0(file_name,".jpg")),
         width=Par$mfrow[2]*MapSizeRatio[2],
         height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
  }
  if(Format%in%c("tif","tiff")){
    tiff(file=paste0(FileName, paste0(file_name,".tif")),
         width=Par$mfrow[2]*MapSizeRatio[2],
         height=Par$mfrow[1]*MapSizeRatio[1], res=Res, units='in')
  }
    if(add==FALSE) par( Par )          # consider changing to Par=list() input, which overloads defaults a la optim() "control" input
    for( tI in 1:ncol(Y_gt) ){
      # Read extrapolation grid
      if( is.numeric(extrapolation_list$zone) ){
        CRS_orig = CRS( paste0("+proj=utm +units=km +zone=",extrapolation_list$zone-ifelse(extrapolation_list$flip_around_dateline,30,0)) )
      }else{
        CRS_orig = CRS(extrapolation_list$zone)
      }
      Points_orig = SpatialPointsDataFrame( coords=spatial_list$loc_g, data=data.frame("y"=Y_gt[,tI]), proj4string=CRS_orig )

      # Reproject to Lat-Long
      Points_LongLat = spTransform( Points_orig, CRS('+proj=longlat') )

      # Re-project to plotting CRS
      Points_proj = spTransform( Points_orig, CRS_proj )

      # Interpolate to raster
      # library(plotKML)
      cell.size = mean(diff(Points_proj@bbox[1,]),diff(Points_proj@bbox[2,]))/floor(sqrt(nrow(Y_gt)))
      Raster_proj = plotKML::vect2rast( Points_proj, cell.size=cell.size )
      image( Raster_proj, col=col, zlim=zlim )

      # Plot maps using rnaturalearth
      #library( rnaturalearth )
      map_info = rnaturalearth::ne_countries(scale=switch( "na", "low"=110, "medium"=50, "high"=10, 50 ))
      map_info = spTransform(map_info, CRSobj=CRS_proj)
      plot( map_info, col=land_color, add=TRUE )              # #EAEAEA

      # Title and box
      title( panel_labels[tI], line=0.1, cex.main=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main), cex=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main) )
      box()
    }
    # Include legend
    if( Legend$use==TRUE ){
      xl = (100-Legend$x[1])/100*min(coordinates(Points_proj)[,1]) + (Legend$x[1])/100*max(coordinates(Points_proj)[,1])
      xr = (100-Legend$x[2])/100*min(coordinates(Points_proj)[,1]) + (Legend$x[2])/100*max(coordinates(Points_proj)[,1])
      yb = (100-Legend$y[1])/100*min(coordinates(Points_proj)[,2]) + (Legend$y[1])/100*max(coordinates(Points_proj)[,2])
      yt = (100-Legend$y[2])/100*min(coordinates(Points_proj)[,2]) + (Legend$y[2])/100*max(coordinates(Points_proj)[,2])
      plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(zlim[1],zlim[2],length=4),1), rect.col=col, cex=1, align=c("lt","rb")[2], gradient="y")
    }
    # Margin text
    if(add==FALSE) mtext(side=1, outer=TRUE, outermargintext[1], cex=1.75, line=par()$oma[1]/2)
    if(add==FALSE) mtext(side=2, outer=TRUE, outermargintext[2], cex=1.75, line=par()$oma[2]/2)
  if(Format %in% c("png","jpg","tif","tiff")) dev.off()
  return( invisible(list("Par"=Par)) )
}
