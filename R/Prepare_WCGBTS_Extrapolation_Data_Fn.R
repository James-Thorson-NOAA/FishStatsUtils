#' @export
Prepare_WCGBTS_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, surveyname='propInWCGBTS', zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( california_current_grid, package="FishStatsUtils" )
  Data_Extrap <- california_current_grid

  # Survey areas
  Area_km2_x = 4*apply(Data_Extrap[,surveyname,drop=FALSE], MARGIN=1, FUN=min)
  
  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=(-1000)*Data_Extrap[,'Depth_km'], "BEST_LAT_DD"=Data_Extrap[,'Lat'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    #a_el[,l] = apply(Tmp, MARGIN=1, FUN=nwfscDeltaGLM::strata.fn, Strata.df=strata.limits[l,,drop=FALSE])
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)
  Data_Extrap = cbind( Data_Extrap, 'Include'=(Data_Extrap[,'Cowcod']==0 & Data_Extrap[,'Ngdc_m']<(-35)))
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('E_km','N_km')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}
