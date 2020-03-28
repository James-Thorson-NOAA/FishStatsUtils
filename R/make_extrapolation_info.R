
#' Build extrapolation grid
#'
#' \code{make_extrapolation_data} builds an object used to determine areas to extrapolation densities to when calculating indices
#'
#' To do area-weighted extrapolation of estimated density for use in calculating abundance indices,
#' it is necessary to have a precise measurement of the footprint for a given survey design.
#' Using VAST, analysts do this by including an "extrapolation grid" where densities are predicted
#' at the location of each grid cell and where each grid cell is associated with a known area within a given survey design.
#' Collaborators have worked with the package author to include the extrapolation-grid for several
#' regions automatically in FishStatsUtils. For new regions an analyst can either (1) detect
#' the grid automatically using \code{Region="Other"}, or (2) input an extrapolation-grid manually
#' using \code{Region="User"}, or supply a GIS shapefile \code{Region="[directory_path/file_name].shp"}.
#' The extrapolation is also used to determine where to drawn pixels when plotting predictions of density.
#' If a user supplies a character-vector with more than one of these, then they are combined to
#' assemble a combined extrapolation-grid.
#'
#' When supplying a shapefile, I recommend using a UTM projection for projargs, which appears to have lower
#' projection errors regarding total area than rnaturalearth.
#'
#' @inheritParams sp::CRS
#' @inheritParams Calc_Kmeans
#'
#' @param Region a character vector, where each element that is matched against potential values to determine the region for the extrapolation grid. Current options are "california_current", "west_coast_hook_and_line", "british_columbia", "eastern_bering_sea", "northern_bering_sea", "bering_sea_slope", "st_matthews_island", "aleutian_islands", "gulf_of_alaska", "northwest_atlantic", "south_africa", "gulf_of_st_lawrence", "new_zealand", "habcam", "gulf_of_mexico", "ATL-IBTS-Q1", "ATL-IBTS-Q4", "BITS", "BTS", "BTS-VIIA", "EVHOE", "IE-IGFS", "NIGFS", "NS_IBTS", "PT-IBTS", "SP-ARSA", "SP-NORTH", "SP-PORC", "stream_network", "user", "other", or the absolute path and file name for a GIS shapefile
#' @param strata.limits an input for determining stratification of indices (see example script)
#' @param zone UTM zone used for projecting Lat-Lon to km distances; use \code{zone=NA} by default to automatically detect UTM zone from the location of extrapolation-grid samples
#' @param flip_around_dateline used applies when using UTM projection, where {flip_around_dateline=TRUE} causes code to convert given latitude on other side of globe (as helpful when data straddle dateline); default value depends upon \code{Region} used
#' @param create_strata_per_region Boolean indicating whether to create a single stratum for all regions listed in \code{Region} (the default), or a combined stratum in addition to a stratum for each individual Region
#' @param observations_LL a matrix with two columns (labeled 'Lat' and 'Lon') giving latitude and longitude for each observation; only used when \code{Region="other"}
#' @param input_grid a matrix with three columns (labeled 'Lat', 'Lon', and 'Area_km2') giving latitude, longitude, and area for each cell of a user-supplied grid; only used when \code{Region="user"}
#' @param grid_dim_km numeric-vector with length two, giving the distance in km between cells in the automatically generated extrapolation grid; only used if \code{Region="other"}
#' @param maximum_distance_from_sample maximum distance that an extrapolation grid cell can be from the nearest sample and still be included in area-weighted extrapolation of density; only used if \code{Region="other"}
#' @param grid_dim_LL same as \code{grid_dim_km} except measured in latitude-longitude coordinates; only used if \code{Region="other"}
#' @param grid_in_UTM Boolean stating whether to automatically generate an extrapolation grid based on sampling locations in km within the UTM projection of within Lat-Lon coordinates; only used if \code{Region="other"}
#' @param strata_to_use strata to include by default for the BC coast extrapolation grid; only used if \code{Region="british_columbia"}
#' @param survey survey to use for New Zealand extrapolation grid; only used if \code{Region="new_zealand"}
#' @param region which coast to use for South Africa extrapolation grid; only used if \code{Region="south_africa"}
#' @param surveyname area of West Coast to include in area-weighted extrapolation for California Current; only used if \code{Region="california_current"}
#' @param max_cells Maximum number of extrapolation-grid cells.  If number of cells in extrapolation-grid is less than this number, then its value is ignored.  Default \code{max_cells=Inf} results in no reduction in number of grid cells from the default extrapolation-grid for a given region.  Using a lower value is particularly useful when \code{fine_scale=TRUE} and using epsilon bias-correction, such that the number of extrapolation-grid cells is often a limiting factor in estimation speed.
#' @param ... other objects passed for individual regions (see example script)

#' @return Tagged list used in other functions
#' \describe{
#'   \item{a_el}{The area associated with each extrapolation grid cell (rows) and strata (columns)}
#'   \item{Data_Extrap}{A data frame describing the extrapolation grid}
#'   \item{zone}{the zone used to convert Lat-Long to UTM by PBSmapping package}
#'   \item{flip_around_dateline}{a boolean stating whether the Lat-Long is flipped around the dateline during conversion to UTM}
#'   \item{Area_km2_x}{the area associated with each row of Data_Extrap, in units square-kilometers}
#' }

#' @export
make_extrapolation_info = function( Region, projargs=NA, zone=NA, strata.limits=data.frame('STRATA'="All_areas"),
  create_strata_per_region=FALSE, max_cells=NULL, input_grid=NULL, observations_LL=NULL, grid_dim_km=c(2,2),
  maximum_distance_from_sample=NULL, grid_in_UTM=TRUE, grid_dim_LL=c(0.1,0.1),
  region=c("south_coast","west_coast"), strata_to_use=c('SOG','WCVI','QCS','HS','WCHG'),
  survey="Chatham_rise", surveyname='propInWCGBTS', flip_around_dateline, nstart=100, ... ){

  # Note: flip_around_dateline must appear in arguments for argument-matching in fit_model
  # However, it requires a different default value for different regions; hence the input format being used.

  # Backwards compatibility for when settings didn't include max_cells, such that settings$max_cells=NULL
  if(is.null(max_cells)) max_cells = Inf

  for( rI in seq_along(Region) ){
    Extrapolation_List = NULL
    if( tolower(Region[rI]) == "california_current" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_WCGBTS_Extrapolation_Data_Fn( strata.limits=strata.limits, surveyname=surveyname, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) %in% c("wcghl","wcghl_domain","west_coast_hook_and_line") ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_WCGHL_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }                      #
    if( tolower(Region[rI]) == "british_columbia" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_BC_Coast_Extrapolation_Data_Fn( strata.limits=strata.limits, strata_to_use=strata_to_use, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "eastern_bering_sea" ){ #
      if(missing(flip_around_dateline)) flip_around_dateline = TRUE
      Extrapolation_List = Prepare_EBS_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "northern_bering_sea" ){ #
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NBS_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "bering_sea_slope" ){ #
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_BSslope_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) %in% c("st_matthews_island","smi") ){ #
      if(missing(flip_around_dateline)) flip_around_dateline = TRUE
      Extrapolation_List = Prepare_SMI_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "aleutian_islands" ){ #
      if(missing(flip_around_dateline)) flip_around_dateline = TRUE
      Extrapolation_List = Prepare_AI_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_alaska" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GOA_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "northwest_atlantic" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NWA_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "south_africa" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_SA_Extrapolation_Data_Fn( strata.limits=strata.limits, region=region, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_st_lawrence" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GSL_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "new_zealand" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NZ_Extrapolation_Data_Fn( strata.limits=strata.limits, survey=survey, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "habcam" ){  #
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_HabCam_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_mexico" ){
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GOM_Extrapolation_Data_Fn( strata.limits=strata.limits, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( toupper(Region[rI]) %in% c("ATL-IBTS-Q1","ATL-IBTS-Q4","BITS","BTS","BTS-VIIA","EVHOE","IE-IGFS","NIGFS","NS_IBTS","PT-IBTS","SP-ARSA","SP-NORTH","SP-PORC") ){
      if( Region[rI]=="SP-ARSA" ) stop("There's some problem with `SP-ARSA` which precludes it's use")
      Conversion = convert_shapefile( file_path=paste0(system.file("region_shapefiles",package="FishStatsUtils"),"/",toupper(Region[rI]),"/Shapefile.shp"),
        projargs=projargs, grid_dim_km=grid_dim_km, ... )
      Extrapolation_List = list( "a_el"=matrix(Conversion$extrapolation_grid[,'Area_km2'],ncol=1), "Data_Extrap"=Conversion$extrapolation_grid,
        "zone"=NA, "projargs"=Conversion$projargs, "flip_around_dateline"=FALSE, "Area_km2_x"=Conversion$extrapolation_grid[,'Area_km2'])
    }
    if( file.exists(Region[rI]) ){
      Conversion = convert_shapefile( file_path=Region[rI], projargs=projargs, grid_dim_km=grid_dim_km, ... )
      Extrapolation_List = list( "a_el"=matrix(Conversion$extrapolation_grid[,'Area_km2'],ncol=1), "Data_Extrap"=Conversion$extrapolation_grid,
        "zone"=NA, "projargs"=Conversion$projargs, "flip_around_dateline"=FALSE, "Area_km2_x"=Conversion$extrapolation_grid[,'Area_km2'])
    }
    if( tolower(Region[rI]) == "stream_network" ){
      if( is.null(input_grid)){
        stop("Because you're using a stream network, please provide 'input_grid' input")
      }
      if( !(all(c("Lat","Lon","Area_km2","child_i") %in% colnames(input_grid))) ){
        stop("'input_grid' must contain columns named 'Lat', 'Lon', 'Area_km2', and 'child_i'")
      }
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_User_Extrapolation_Data_Fn( strata.limits=strata.limits, input_grid=input_grid, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( tolower(Region[rI]) == "user" ){
      if( is.null(input_grid)){
        stop("Because you're using a user-supplied region, please provide 'input_grid' input")
      }
      if( !(all(c("Lat","Lon","Area_km2") %in% colnames(input_grid))) ){
        stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
      }
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_User_Extrapolation_Data_Fn( strata.limits=strata.limits, input_grid=input_grid, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }
    if( is.null(Extrapolation_List) ){
      if( is.null(observations_LL)){
        stop("Because you're using a new Region[rI], please provide 'observations_LL' input")
      }
      if(missing(flip_around_dateline)) flip_around_dateline = FALSE
      Extrapolation_List = Prepare_Other_Extrapolation_Data_Fn( strata.limits=strata.limits, observations_LL=observations_LL,
        grid_dim_km=grid_dim_km, maximum_distance_from_sample=maximum_distance_from_sample,
        grid_in_UTM=grid_in_UTM, grid_dim_LL=grid_dim_LL, projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline, ... )
    }

    # Combine
    if( rI==1 ){
      Return = Extrapolation_List
    }else{
      Return = combine_extrapolation_info( Return, Extrapolation_List, create_strata_per_region=create_strata_per_region )
    }
  }

  # Optionally reduce number of extrapolation-grid cells
  if( max_cells < nrow(Return$Data_Extrap) ){
    message( "# Reducing extrapolation-grid from ",nrow(Return$Data_Extrap)," to ",max_cells," cells for Region(s): ",paste(Region,collapse=", ") )
    # Run K-means only on grid-cells with nonzero area
    loc_orig = Return$Data_Extrap[,c("E_km","N_km")]
      loc_orig = loc_orig[ which(Return$Area_km2_x>0), ]
    Kmeans = Calc_Kmeans( n_x=max_cells, loc_orig=loc_orig, nstart=nstart,
      randomseed=1, iter.max=1000, DirPath=paste0(getwd(),"/"), Save_Results=FALSE )
    Kmeans[["cluster"]] = RANN::nn2( data=Kmeans[["centers"]], query=Return$Data_Extrap[,c("E_km","N_km")], k=1)$nn.idx[,1]
    # Transform Extrapolation_List
    aggregate_vector = function( values_x, index_x, max_index, FUN=sum ){
      tapply( values_x, INDEX=factor(index_x,levels=1:max_index), FUN=FUN )
    }
    # a_el
    a_el = matrix(NA, nrow=max_cells, ncol=ncol(Return$a_el) )
    for(lI in 1:ncol(Return$a_el)){
      a_el[,lI] = aggregate_vector( values_x=Return$a_el[,lI], index_x=Kmeans$cluster, max_index=max_cells )
    }
    # Area_km2_x
    Area_km2_x = aggregate_vector( values_x=Return$Area_km2_x, index_x=Kmeans$cluster, max_index=max_cells )
    # Data_Extrap
    Include = aggregate_vector( values_x=Return$Data_Extrap[,'Include'], index_x=Kmeans$cluster, max_index=max_cells, FUN=function(vec){any(vec>0)} )
    #Area_in_survey_km2 = aggregate_vector( values_x=Return$Data_Extrap[,'Area_in_survey_km2'], index_x=Kmeans$cluster, max_index=max_cells )
    lonlat_g = project_coordinates( X=Kmeans$centers[,"E_km"], Y=Kmeans$centers[,"N_km"], projargs="+proj=longlat +ellps=WGS84", origargs=Return$projargs )
    #Data_Extrap = cbind( "Lon"=lonlat_g[,1], "Lat"=lonlat_g[,2], "Area_in_survey_km2"=Area_in_survey_km2, "Include"=Include, Kmeans$centers )
    Data_Extrap = cbind( "Lon"=lonlat_g[,1], "Lat"=lonlat_g[,2], "Include"=Include, Kmeans$centers )
    # Assemble
    Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=Return$zone, "projargs"=Return$projargs,
      "flip_around_dateline"=Return$flip_around_dateline, "Area_km2_x"=Area_km2_x )
  }

  # Add total across regions if requested
  if( length(Region)>1 & create_strata_per_region==TRUE ){
    Return$a_el = cbind( "Total"=rowSums(Return$a_el), Return$a_el )
  }

  # Return
  class(Return) = "make_extrapolation_info"
  return( Return )
}

#' Plot extrapolation-grid used for spatial inference
#'
#' @inheritParams sp::CRS
#'
#' @title Plot extrapolation-grid
#' @param x Output from \code{\link{make_extrapolation_info}}
#' @param ... Not used
#' @return NULL
#' @method plot make_extrapolation_info
#' @export
plot.make_extrapolation_info <- function(x, cex=0.01, land_color="grey", map_resolution="medium", ...)
{
  par( mfrow=c(1,2), mar=c(3,3,2,0), mgp=c(1.75,0.25,0) )

  # CRS for original and new projections
  CRS_orig = sp::CRS( '+proj=longlat' )
  CRS_proj = sp::CRS( x$projargs )

  # Data for mapping
  map_data = rnaturalearth::ne_countries(scale=switch(map_resolution, "low"=110, "medium"=50, "high"=10, 50 ))

  # Plot #1 -- Latitude
  plot( x$Data_Extrap[which(x$Area_km2_x>0),c('Lon','Lat')], cex=cex, main="Extrapolation (Lat-Lon)", ... )
  map_data_orig = sp::spTransform(map_data, CRSobj=CRS_orig)
  sp::plot( map_data_orig, col=land_color, add=TRUE )

  # Plot #2 -- Projection coordinates
  if( !any(is.na(x$Data_Extrap[,c('E_km','N_km')])) ){
    plot( x$Data_Extrap[which(x$Area_km2_x>0),c('E_km','N_km')], cex=cex, main="Extrapolation (North-East)", ... )
    #map_data_proj = sp::spTransform(map_data, CRSobj=CRS_proj)
    #sp::plot( map_data_proj, col=land_color, add=TRUE )
  }

  invisible(x)
}

#' Extract extrapolation-grid used for spatial inference
#'
#' @title Extract extrapolation-grid
#' @param x Output from \code{\link{make_extrapolation_info}}
#' @param quiet Boolean indicating whether to print to terminal or not
#' @param ... Not used
#' @return Data frame of Latitude and Longitude for each extrapolation-grid cell
#' @method print make_extrapolation_info
#' @export
print.make_extrapolation_info <- function(x, quiet=FALSE, ...)
{
  loc_gz = cbind( x$Data_Extrap[,c('Lon','Lat')], "Area_km2"=x$Area_km2_x )[ which(x$Area_km2_x>0), ]
  rownames(loc_gz) = paste0("Grid_",1:nrow(loc_gz))
  #ans = c( x[c('zone','flip_around_dateline')], list("extrapolation_grid"=loc_gz) )

  if(quiet==FALSE){
    cat("make_extrapolation_info(.) result\n")
    print( summary(loc_gz) )
    if(is.na(x$zone)){
      cat( paste0("\nprojargs: ",x$projargs,"\n") )
    }else{
      cat( paste0("\nUTM zone: ",x$zone,"\n") )
    }
    if( x$flip_around_dateline == TRUE ){
      cat( "Note: Longitude was translated away from dateline (by adding 180, to avoid dateline projection issues in PBSmapping) prior projection using the UTM zone listed above\n" )
    }
  }

  invisible(loc_gz)
}
