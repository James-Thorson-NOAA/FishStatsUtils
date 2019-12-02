#' @export
Prepare_NWA_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = list('All_areas'=1:1e5)
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( northwest_atlantic_grid, package="FishStatsUtils" )
  Data_Extrap <- northwest_atlantic_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']
  
  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  if( length(strata.limits)==1 && strata.limits[1]=="EPU" ){
    # Specify strata by 'stratum_number'
    a_el = matrix(NA, nrow=nrow(Data_Extrap), ncol=length(unique(northwest_atlantic_grid[,'EPU'])), dimnames=list(NULL,unique(northwest_atlantic_grid[,'EPU'])) )
    for(l in 1:ncol(a_el)){
      a_el[,l] = ifelse( Data_Extrap[,'EPU']==unique(northwest_atlantic_grid[,'EPU'])[l], Area_km2_x, 0 )
    }
  }else{
    # Specify strata by 'stratum_number'
    a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
    for(l in 1:ncol(a_el)){
      a_el[,l] = ifelse( Data_Extrap[,'stratum_number'] %in% strata.limits[[l]], Area_km2_x, 0 )
    }
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('E_km','N_km')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}
