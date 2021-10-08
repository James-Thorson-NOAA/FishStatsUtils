#' @title
#' Plot results on a multi-panel figure of maps
#'
#' @description
#' \code{plot_variable} plots a map and fills in regions with colors to represent intensity in an areal-interpretion of model results
#'
#' If maps are generating visual artefacts, please try using argument \code{country} to simplify the polygons used to represent land features.
#'
#' In some sampling designs and/or analyses, it is not appropriate to extrapolate beyond samping stations.  In these cases,
#' the analyst may define the extrapolation-grid to only represent sampling stations, and then use \code{format="points"}
#' to restrict mapping to bullets at those locations.
#'
#' @inheritParams project_coordinates
#' @inheritParams rnaturalearth::ne_countries
#' @inheritParams raster::rasterize
#'
#' @param Y_gt matrix where values for every column are plotted as a map in a panel-figure with one panel per column
#' @param projargs Character passed to \code{\link[sp]{CRS}} and see that documentation for description.
#'        See \url{https://proj.org/operations/projections/index.html} for a list of projections to pass via \code{projargs}.
#'        As a stable default please consider \code{projargs='+proj=natearth +lon_0=0 +units=km'}
#'        where argument \code{+lon_0} allows the user to center eastings on a specified longitude.
#' @param n_cells Integer used to determine the argument \code{cell.size} passed
#'        to \code{\link[plotKML]{vect2rast}} when converting output at extrapolation-grid cells
#'        to a raster prior to plotting mapped outputs using \code{\link{plot_variable}}.  An increased \code{n_cells} results in a decreased
#'        \code{cell.size} and hence higher resolution plots.  The default value for \code{n_cells}
#'        will often need to be modified for a given purpose.
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param legend_x two numeric values (generally between 0 and 1, but slightly lower/higher values
#'        generate colorbars that are just outside the plotting window) giving left and right-hand location of color legend
#' @param legend_y two numeric values (see description for legend_x) giving bottom and top location of color legend
#' @param map_list output from \code{\link{make_map_info}}, which is included in output of \code{\link{plot_results}}
#' @param zlim two numeric values, specifying range for defining bounds of color scale.
#'        If \code{zlim=NULL}, then a constant scale is inferred
#'        from the range of the plotted variable and a color-legend is plotted in the last panel.
#'        If \code{zlim=NA} then a different range is used in each panel from the range of \code{Y_gt[,t]} and a color-legend is plotted in every panel.
#' @param add boolean indicating whether to add plot to an existing panel figure, or to define a new panel figure
#' @param outermargintext vector defining text to plot in outer margins of panel figure
#' @param panel_labels vector defining titles to use for each panel; defaults to blank
#' @param contour_nlevels number of levels used when adding contour lines, passed to \code{\link[graphics]{contour}} as argument \code{nlevels}
#' @param format Character specifying whether to plot maps as a raster (the default), \code{format="raster"}, or as points \code{format="points"}
#' @param cex.points Numeric specifying the size of bullets when \code{format="points"}
#' @param legend_digits number of digits (i.e., value passed to \code{round}) when creating labels for colorbar legend
#'
#' @param ... arguments passed to \code{par}
#'

#' @export
plot_variable <-
function( Y_gt,
          map_list,
          panel_labels,
          projargs = '+proj=longlat',
          map_resolution = "medium",
          file_name = "density",
          working_dir = paste0(getwd(),"/"),
          Format = "png",
          Res = 200,
          add = FALSE,
          outermargintext = c("Eastings","Northings"),
          zlim = NULL,
          col = viridisLite::viridis,
          mar = c(0,0,2,0),
          oma = c(4,4,0,0),
          legend_x = c(0,0.05),
          legend_y = c(0.05,0.45),
          cex.legend = 1,
          mfrow,
          land_color = "grey",
          n_cells = NULL,
          xlim,
          ylim,
          country = NULL,
          contour_nlevels = 0,
          fun = mean,
          format = "raster",
          cex.points = 1,
          legend_digits = 1,
          ...){

  ###################
  # Settings and inputs
  ###################

  # Check for problems and fill in missing stuff
  if( is.vector(Y_gt)){
    Y_gt = matrix(Y_gt, ncol=1)
  }
  if( is.null(zlim)){
    zlim = range(Y_gt, na.rm=TRUE)
  }
  if( missing(map_list) || is.null(map_list$MapSizeRatio) ){
    MapSizeRatio = c(3, 3)
  }else{
    MapSizeRatio = map_list$MapSizeRatio
  }
  if( !("PlotDF" %in% names(map_list)) ) stop("Check input `map_list`")
  Y_gt = Y_gt[ map_list$PlotDF[which(map_list$PlotDF[,'Include']>0),'x2i'], , drop=FALSE]
  if(missing(n_cells) || is.null(n_cells)) n_cells = nrow(Y_gt)
  if( missing(mfrow) ){
    mfrow = ceiling(sqrt(ncol(Y_gt)))
    mfrow = c( mfrow, ceiling(ncol(Y_gt)/mfrow) )
  }
  if( missing(panel_labels) | is.null(panel_labels) ){
    panel_labels = rep("", ncol(Y_gt))
  }
  if( length(panel_labels) != ncol(Y_gt) ){
    warning( "panel_labels and `ncol(Y_gt)` don't match: Changing panel_labels'")
    panel_labels = 1:ncol(Y_gt)
  }
  if( is.null(col)){
    col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
  }
  if( is.function(col)){
    col = col(1000)
  }
  if( !any(is.na(c(legend_x,legend_y))) ){
    if( any(c(legend_x,legend_y) > 1.2) | any(c(legend_x,legend_y) < -0.2) ){
      stop("Check values for `legend_x` and `legend_y`")
    }
  }
  # Location of extrapolation-grid cells
  loc_g = map_list$PlotDF[which(map_list$PlotDF[,'Include']>0),c('Lon','Lat')]

  # CRS for original and new projections
  CRS_orig = sp::CRS( '+proj=longlat' )
  CRS_proj = sp::CRS( projargs )

  # Data for mapping
  #map_data = rnaturalearth::ne_coastline(scale=switch(map_resolution, "low"=110, "medium"=50, "high"=10, 50), continent="america")
  map_data = rnaturalearth::ne_countries(scale=switch(map_resolution, "low"=110, "medium"=50, "high"=10, 50), country=country)
  # Fix warning messages from projecting rnaturalearth object
  # Solution: Recreate SpatialPolygonsDataFrame from output
  map_data = sp::SpatialPolygonsDataFrame( Sr=sp::SpatialPolygons(slot(map_data,"polygons"),proj4string=CRS_orig), data=slot(map_data,"data") )
  # comment(slot(map_data, "proj4string")) =  comment(sp::CRS("+proj=longlat"))
  map_proj = sp::spTransform(map_data, CRSobj=CRS_proj)

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

    # Get Zlim
    Zlim = zlim
    if(is.na(Zlim[1])) Zlim = range(Y_gt[,tI],na.rm=TRUE)
    if(missing(xlim)) xlim = Points_proj@bbox[1,]
    if(missing(ylim)) ylim = Points_proj@bbox[2,]

    # Do plot
    if( format=="raster"  ){
      if( all(is.na(Y_gt[,tI])) ){
        # Empty plot if no data
        plot( x=Points_proj@coords[,1], y=Points_proj@coords[,2], type="n", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, xlab="", ylab="" )
      }else{
        # Interpolate and plot as raster
        cell.size = mean(diff(Points_proj@bbox[1,]),diff(Points_proj@bbox[2,])) / floor(sqrt(n_cells))
        Raster_proj = plotKML::vect2rast( Points_proj, cell.size=cell.size, fun=fun )
        plot( x=Points_proj@coords[,1], y=Points_proj@coords[,2], type="n", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, xlab="", ylab="" )
        image( Raster_proj, col=col, zlim=Zlim, add=TRUE )
        # Add contour lines
        if( contour_nlevels > 0 ){
          contour( Raster_proj, add=TRUE, nlevels=contour_nlevels )
        }
      }
    }else if( format=="points" ){
      # Plot points
      Points_col = col[cut(Y_gt[,tI],breaks=seq(Zlim[1],Zlim[2],length=length(col)),include.lowest=TRUE)]
      points( x=Points_proj@coords[,1], y=Points_proj@coords[,2], col=Points_col, pch=20, cex=cex.points, xaxt="n", yaxt="n" )
    }

    # Plot maps using rnaturalearth
    # Make plot after adding raster or points, to overwrite land
    sp::plot( map_proj, col=land_color, add=TRUE )

    # Title and box
    title( panel_labels[tI], line=0.1, cex.main=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main), cex=ifelse(is.null(Par$cex.main), 1.5, Par$cex.main) )
    box()

    # Include legend
    if( !any(is.na(c(legend_x,legend_y))) & (tI==ncol(Y_gt) | is.na(zlim[1])) ){
      xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
      xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
      yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
      yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
      if( diff(legend_y) > diff(legend_x) ){
        align = c("lt","rb")[2]
        gradient = c("x","y")[2]
      }else{
        align = c("lt","rb")[1]
        gradient = c("x","y")[1]
      }
      plotrix::color.legend(xl=xl, yb=yb, xr=xr, yt=yt, legend=round(seq(Zlim[1],Zlim[2],length=4),legend_digits), rect.col=col, cex=cex.legend, align=align, gradient=gradient)
    }
  }

  # Margin text
  if(add==FALSE) mtext(side=1, outer=TRUE, outermargintext[1], cex=1.75, line=par()$oma[1]/2)
  if(add==FALSE) mtext(side=2, outer=TRUE, outermargintext[2], cex=1.75, line=par()$oma[2]/2)

  # return stuff as necessary
  return( invisible(list("Par"=Par, "cell.size"=cell.size, "n_cells"=n_cells, "xlim"=xlim, "ylim"=ylim)) )
}
