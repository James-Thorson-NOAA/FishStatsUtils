
#' Build extrapolation grid
#'
#' \code{make_extrapolation_data} builds an object used to determine areas to extrapolation densities to when calculating indices
#'
#' To do area-weighted extrapolation of estimated density for use in calculating abundance indices, it is necessary to have a precise measurement of the footprint for a given survey design. Using VAST, analysts do this by including an "extrapolation grid" where densities are predicted at the location of each grid cell and where each grid cell is associated with a known area within a given survey design. Collaborators have worked with the package author to include the extrapolation-grid for several regions automatically in FishStatsUtils, but for new regions an analyst must either detect the grid automatically using \code{Region="Other"} or input an extrapolation-grid manually using \code{Region="User"}.  The extrapolation is also used to determine where to drawn pixels when plotting predictions of density.
#'
#' @param Region a character vector, where each element that is matched against potential values to determine the region for the extrapolation grid. Current options are "clifornia_current", "west_coast_hook_and_line", "british_columbia", "eastern_bering_sea", "northern_bering_sea", "bering_sea_slope", "st_matthews_island", "aleutian_islands", "gulf_of_alaska", "northwest_atlantic", "south_africa", "gulf_of_st_lawrence", "new_zealand", "habcam", "gulf_of_mexico", "user", or "other"
#' @param strata.limits an input for determining stratification of indices (see example script)
#' @param zone UTM zone used for projecting Lat-Lon to km distances; use \code{zone=NA} by default to automatically detect UTM zone from the location of extrapolation-grid samples
#' @param observations_LL a matrix with two columns (labeled 'Lat' and 'Lon') giving latitude and longitude for each observation (only used when Region doesn't match known entries)
#' @param input_grid a matrix with three columns (labeled 'Lat', 'Lon', and 'Area_km2') giving latitude, longitude, and area for each cell of a user-supplied grid (only used when \code{Region=="User"})
#' @param grid_dim_km numeric-vector with length two, giving the distance in km between cells in the automatically generated extrapolation grid; only used if \code{Region} doesn't match existing options
#' @param maximum_distance_from_sample maximum distance that an extrapolation grid cell can be from the nearest sample and still be included in area-weighted extrapolation of density; only used if \code{Region} doesn't match existing options
#' @param grid_dim_LL same as \code{grid_dim_km} except measured in latitude-longitude coordinates; only used if \code{Region} doesn't match existing options
#' @param grid_in_UTM Boolean stating whether to automatically generate an extrapolation grid based on sampling locations in km within the UTM projection of within Lat-Lon coordinates; only used if \code{Region} doesn't match existing options
#' @param strata_to_use strata to include by default for the BC coast extrapolation grid; only used if \code{Region="british_columbia"}
#' @param survey survey to use for New Zealand extrapolation grid; only used if \code{Region="new_zealand"}
#' @param region which coast to use for South Africa extrapolation grid; only used if \code{Region="south_africa"}
#' @param surveyname area of West Coast to include in area-weighted extrapolation for California Current; only used if \code{Region="california_current"}
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
make_extrapolation_info = function( Region, zone=NA, strata.limits=data.frame('STRATA'="All_areas"),
  input_grid=NULL, observations_LL=NULL, grid_dim_km=c(2,2), maximum_distance_from_sample=NULL,
  grid_in_UTM=TRUE, grid_dim_LL=c(0.1,0.1), region=c("south_coast","west_coast"),
  strata_to_use=c('SOG','WCVI','QCS','HS','WCHG'), survey="Chatham_rise", surveyname='propInWCGBTS', ... ){

  for( rI in seq_along(Region) ){
    Extrapolation_List = NULL
    if( tolower(Region[rI]) == "california_current" ){
      Extrapolation_List = Prepare_WCGBTS_Extrapolation_Data_Fn( strata.limits=strata.limits, surveyname=surveyname, zone=zone, ... )
    }
    if( tolower(Region[rI]) %in% c("wcghl","wcghl_domain","west_coast_hook_and_line") ){
      Extrapolation_List = Prepare_WCGHL_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }                      #
    if( tolower(Region[rI]) == "british_columbia" ){
      Extrapolation_List = Prepare_BC_Coast_Extrapolation_Data_Fn( strata.limits=strata.limits, strata_to_use=strata_to_use, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "eastern_bering_sea" ){ #
      Extrapolation_List = Prepare_EBS_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "northern_bering_sea" ){ #
      Extrapolation_List = Prepare_NBS_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "bering_sea_slope" ){ #
      Extrapolation_List = Prepare_BSslope_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) %in% c("st_matthews_island","smi") ){ #
      Extrapolation_List = Prepare_SMI_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "aleutian_islands" ){ #
      Extrapolation_List = Prepare_AI_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_alaska" ){
      Extrapolation_List = Prepare_GOA_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "northwest_atlantic" ){
      Extrapolation_List = Prepare_NWA_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "south_africa" ){
      Extrapolation_List = Prepare_SA_Extrapolation_Data_Fn( strata.limits=strata.limits, region=region, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_st_lawrence" ){
      Extrapolation_List = Prepare_GSL_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "new_zealand" ){
      Extrapolation_List = Prepare_NZ_Extrapolation_Data_Fn( strata.limits=strata.limits, survey=survey, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "habcam" ){  #
      Extrapolation_List = Prepare_HabCam_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "gulf_of_mexico" ){
      Extrapolation_List = Prepare_GOM_Extrapolation_Data_Fn( strata.limits=strata.limits, zone=zone, ... )
    }
    if( tolower(Region[rI]) == "user" ){
      if( is.null(input_grid)){
        stop("Because you're using a user-supplied region, please provide 'input_grid' input")
      }
      if( !(all(c("Lat","Lon","Area_km2") %in% colnames(input_grid))) ){
        stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
      }
      Extrapolation_List = Prepare_User_Extrapolation_Data_Fn( strata.limits=strata.limits, input_grid=input_grid, zone=zone, ... )
    }
    if( is.null(Extrapolation_List) ){
      if( is.null(observations_LL)){
        stop("Because you're using a new Region[rI], please provide 'observations_LL' input")
      }
      Extrapolation_List = Prepare_Other_Extrapolation_Data_Fn( strata.limits=strata.limits, observations_LL=observations_LL,
        grid_dim_km=grid_dim_km, maximum_distance_from_sample=maximum_distance_from_sample,
        grid_in_UTM=grid_in_UTM, grid_dim_LL=grid_dim_LL, zone=zone, ... )
    }

    # Combine
    if( rI==1 ){
      Return = Extrapolation_List
    }else{
      Return = combine_extrapolation_info( Return, Extrapolation_List )
    }
  }

  # Return
  return( Return )
}
