
#' Convert from Lat-Long to UTM
#'
#' \code{project_coordinates} converts from Latitude-Longitude to Universal Transverse Mercator projections for a given location
#'
#' @param Lat vector of latitudes
#' @param Lat vector of latitudes
#' @param Lon vector of longitudes
#' @param zone UTM zone (integer between 1 and 60) or alphanumeric CRS code used by package rgdal to convert latitude-longitude coordinates to projection in kilometers; \code{zone=NA} uses UTM and automatically detects the appropriate zone
#' @param flip_around_dateline boolean specifying whether to flip Lat-Lon locations around the dateline, and then retransform back (only useful if Lat-Lon straddle the dateline)

#' @return A data frame with the following columns
#' \describe{
#'   \item{X}{The UTM eastings for each value of Lon}
#'   \item{Y}{The UTM northings measured from the equator for each Lat}
#' }

#' @export
project_coordinates <-
function( Lon_i, Lat_i, projargs=NA, origargs="+proj=longlat +ellps=WGS84", zone=NA, flip_around_dateline=FALSE, silent=FALSE ){

  # Original projection
  origCRS = sp::CRS(origargs)

  # Default -- Use zone value if supplied for backwards compatibility
  if( !is.na(zone) ){
    if( !is.na(projargs) ){
      stop("Please do not specify both `zone` and `projargs`")
    }
    # Transform around dateline if requested
    if( flip_around_dateline==TRUE ){
      zone = zone + 30
      Lon_i = Lon_i + 180
    }
    # Project to UTM
    lonlat_sp = sp::SpatialPoints( coords=cbind(Lon_i,Lat_i), proj4string=origCRS ) # expects in long-lat format
    projargs = paste0("+proj=utm +datum=WGS84 +units=km",ifelse(mean(Lat_i)<0," +south","")," +zone=",zone)
    if(silent==FALSE) message("For the UTM conversion, used zone ",zone," as specified")
  }

  # New option
  if( is.na(zone) ){
    # Auto-detect the UTM zone to match PBSmapping
    if( is.na(projargs) ){
      lon1 = Lon_i
      lon2 = ifelse( lon1 < 0, lon1 + 360, lon1 )
      # Determine whether to transform around dateline or not
      if( diff(range(lon1)) < diff(range(lon2)) ){
        mean_lon = mean(lon1)
      }else{
        mean_lon = mean(lon2)
      }
      # Detect UTM based on mean of longitudes to match PBSmapping behavior
      zone = 1 + floor( (mean_lon-180)/6 )
      zone = (zone-1) %% 60 + 1
      projargs = paste0("+proj=utm +datum=WGS84 +units=km",ifelse(mean(Lat_i)<0," +south","")," +zone=",zone)
      if(silent==FALSE) message("For the UTM conversion, automatically detected zone ",zone,".")
    }
    lonlat_sp = sp::SpatialPoints( coords=cbind(Lon_i,Lat_i), proj4string=origCRS ) # expects in long-lat format
  }

  # Check for issues
  if( !is.na(zone) && (zone<1 | zone>60) ) stop("Check zone in `project_coordinates(.)`")
  projCRS = sp::CRS( projargs )

  # Conduct projection
  utm = sp::spTransform( lonlat_sp, projCRS )
  utm_i = utm@coords
  colnames(utm_i) = c("X", "Y")
  if( length(grep("+proj=utm",projargs))==1 & length(grep("+units=km",projargs))==1 ){
    colnames(utm_i) = c("E_km", "N_km")
  }
  if( length(grep("+proj=longlat",projargs))==1 ){
    colnames(utm_i) = c("Lon", "Lat")
  }

  # Add attributes such that it includes everything in `PBSmapping::convUL`, i.e., attribute "zone"
  attr(utm_i,"zone") = zone
  attr(utm_i,"origargs") = origargs
  attr(utm_i,"projargs") = projargs
  attr(utm_i,"origCRS") = origCRS
  attr(utm_i,"projCRS") = projCRS

  # Return results
  return( utm_i )
}
