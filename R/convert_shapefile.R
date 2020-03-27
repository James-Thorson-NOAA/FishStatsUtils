

#' Convert shapefile to extrapolation-grid
#'
#' \code{convert_shapefile} reads in a shapefile and creates an extrapolation with a chosen resolution
#'
#' @inheritParams sp::CRS
#' @inheritParams make_extraploation_info
#'
#' @param file_path path for shapefile on harddrive
#' @param make_plots Boolean indicating whether to visualize inputs and outputs as maps

#' @return extrapolation-grid

#' @examples
#' \dontrun{
#'  convert_shapefile( file_path="C:/Users/James.Thorson/Desktop/Work files/AFSC/2020-03 -- Add ICES grids/IBTS grids/BITS/Shapefile.shp", make_plots=TRUE )
#' }

#' @export
convert_shapefile = function( file_path, projargs=NULL, grid_dim_km=c(2,2), make_plots=FALSE, quiet=TRUE ){

  # Read shapefile
  shapefile = rgdal::readOGR( file_path, verbose=!quiet )
  proj_orig = "+proj=longlat +ellps=WGS84 +no_defs"
  shapefile@proj4string = sp::CRS(proj_orig)

  # Infer projargs if missing, and project
  if(is.null(projargs)){
    utm_zone = floor((mean(sp::bbox(shapefile)[1,]) + 180) / 6) + 1
    projargs = paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +datum=WGS84 +units=km +no_defs ")
  }
  shapefile_proj = sp::spTransform(shapefile, CRSobj=sp::CRS(projargs) )

  # Determine bounds for box
  bbox = shapefile_proj@bbox
  bbox[,1] = floor(bbox[,1] / grid_dim_km[1]) * grid_dim_km[1]
  bbox[,2] = ceiling(bbox[,2] / grid_dim_km[2]) * grid_dim_km[2]

  # Make grid
  grid_proj = expand.grid(X=seq(bbox[1,1],bbox[1,2],by=grid_dim_km[1]), Y=seq(bbox[2,1],bbox[2,2],by=grid_dim_km[2]))
  grid_proj = sp::SpatialPointsDataFrame( coords=grid_proj[,c("X","Y")], data=grid_proj[,1:2], proj4string=sp::CRS(projargs) )

  # Restrict points to those within a polygon
  grid_proj@data = sp::over(grid_proj, shapefile_proj)
  grid_proj = subset( grid_proj, !is.na(AreaName) )

  #
  grid_orig = sp::spTransform( grid_proj, CRSobj=sp::CRS(proj_orig) )
  grid_orig = as.data.frame(grid_orig)
  grid_output = data.frame( Lat=grid_orig[,'Y'], Lon=grid_orig[,'X'], Area_km2=prod(grid_dim_km), grid_orig[,'AreaName'] )

  # make plots
  if( make_plots==TRUE ){
    par( mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

    # Plot projected shapefile
    sp::plot( shapefile_proj, main="shapefile in original projection" )
    axis(1); axis(2); box()

    # Plot projected grid
    sp::plot( grid_proj, main="Grid in new projection" )
    axis(1); axis(2); box()

    # plot original shapefile
    sp::plot( shapefile, main="shapefile in new projection" )
    axis(1); axis(2); box()

    # plot original shapefile
    sp::plot( grid_output[,c('Lon','Lat')], main="Grid in original coordinates" )
    axis(1); axis(2); box()
  }

  # Return output
  return( grid_output )
}

