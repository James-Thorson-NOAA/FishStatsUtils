
#' Combines multiple extrapolation grids
#'
#' \code{combine_extrapolation_info} combines multiple extrapolation grids when combining data from multiple surveysn
#'
#' @param ... a sequence of outputs from \code{FishStatsUtils::make_extrapolation_info}

#' @return Identical output from \code{FishStatsUtils::make_extrapolation_info}, but combined from each input

#' @export
combine_extrapolation_info = function( ... ){

  input_list = list( ... )

  for( lI in 1:length(input_list)){
    if( !all( c("a_el","Data_Extrap","zone","flip_around_dateline","Area_km2_x") %in% names(input_list[[lI]]) )){
      stop( "Check inputs to `combine_extrapolation_grids`" )
    }
  }
  Zone = sapply( input_list, FUN=function(List){List[["zone"]]} )
  Flip = sapply( input_list, FUN=function(List){List[["flip_around_dateline"]]} )
  if( sd(Zone)>0 | sd(Flip)>0 ){
    stop( "Must use same Zone for UTM conversion for all extrapolation grids" )
  }

  # Combine stuff
  a_el = Data_Extrap = Area_km2_x = NULL
  #assign( x="input_list", value=input_list, envir = .GlobalEnv )

  for( lI in 1:length(input_list) ){
    Tmp = input_list[[lI]]$Data_Extrap
    #colnames(Tmp) = ifelse( colnames(Tmp)=="Area_in_survey_km2", "Area_km2", colnames(Tmp) )
    #Data_Extrap = rbind( Data_Extrap, Tmp[,c('E_km','N_km','Lon','Lat','Include','Area_km2')] )
    Data_Extrap = rbind( Data_Extrap, Tmp[,c('E_km','N_km','Lon','Lat','Include')] )
    a_el = rbind( a_el, input_list[[lI]]$a_el )
    Area_km2_x = c( Area_km2_x, input_list[[lI]]$Area_km2_x )
  }

  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=Zone[1], "flip_around_dateline"=Flip[1], "Area_km2_x"=Area_km2_x)
}

