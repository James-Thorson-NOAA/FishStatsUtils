
#' Load an example data set
#'
#' \code{load_example} loads a catch and effort data set from multiple surveys worldwide
#'
#' @param data_set data set to load

#' @export
load_example = function( data_set="EBS_pollock" ){

  # Determine region for a given data set
  region = switch( data_set, "Chatham_rise_hake"="New_Zealand",
                   "WCGBTS_canary"="California_current",
                   "GSL_american_plaice"="Gulf_of_St_Lawrence",
                   "BC_pacific_cod"="British_Columbia",
                   "EBS_pollock"="Eastern_Bering_Sea",
                   "GOA_Pcod"="Gulf_of_Alaska",
                   "GOA_pollock"="Gulf_of_Alaska",
                   "GB_spring_haddock"="Northwest_Atlantic",
                   "GB_fall_haddock"="Northwest_Atlantic",
                   "SAWC_jacopever"="South_Africa",
                   "Aleutian_islands_POP"="Aleutian_Islands",
                   "Other")

  if(data_set=="WCGBTS_canary"){
    data( WCGBTS_Canary_example, package="FishStatsUtils" )
    Year = as.numeric(sapply(WCGBTS_Canary_example[,'PROJECT_CYCLE'], FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]}))
    sampling_data = data.frame( "Catch_KG"=WCGBTS_Canary_example[,'HAUL_WT_KG'], "Year"=Year, "Vessel"=paste(WCGBTS_Canary_example[,"VESSEL"],Year,sep="_"), "AreaSwept_km2"=WCGBTS_Canary_example[,"AREA_SWEPT_HA"]/1e2, "Lat"=WCGBTS_Canary_example[,'BEST_LAT_DD'], "Lon"=WCGBTS_Canary_example[,'BEST_LON_DD'], "Pass"=WCGBTS_Canary_example[,'PASS']-1.5)
  }
  if( data_set %in% c("BC_pacific_cod")){
    data( BC_pacific_cod_example, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=BC_pacific_cod_example[,'PCOD_WEIGHT'], "Year"=BC_pacific_cod_example[,'Year'], "Vessel"="missing", "AreaSwept_km2"=BC_pacific_cod_example[,'TOW.LENGTH..KM.']/100, "Lat"=BC_pacific_cod_example[,'LAT'], "Lon"=BC_pacific_cod_example[,'LON'], "Pass"=0)
  }
  if( data_set %in% c("GSL_american_plaice")){
    data( GSL_american_plaice, package="FishStatsUtils" )
    Print_Message( "GSL_american_plaice" )
    sampling_data = data.frame( "Year"=GSL_american_plaice[,'year'], "Lat"=GSL_american_plaice[,'latitude'], "Lon"=GSL_american_plaice[,'longitude'], "Vessel"="missing", "AreaSwept_km2"=GSL_american_plaice[,'swept'], "Catch_KG"=GSL_american_plaice[,'biomass']*GSL_american_plaice[,'vstd'] )
  }
  if(data_set=="EBS_pollock"){
    data( EBS_pollock_data, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=EBS_pollock_data[,'catch'], "Year"=EBS_pollock_data[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=EBS_pollock_data[,'lat'], "Lon"=EBS_pollock_data[,'long'], "Pass"=0)
  }
  if(data_set=="GOA_Pcod"){
    data( GOA_pacific_cod , package="FishStatsUtils")
    sampling_data = data.frame( "Catch_KG"=GOA_pacific_cod[,'catch'], "Year"=GOA_pacific_cod[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_pacific_cod[,'lat'], "Lon"=GOA_pacific_cod[,'lon'], "Pass"=0)
  }
  if(data_set=="GOA_pollock"){
    data( GOA_walleye_pollock, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=GOA_walleye_pollock[,'catch'], "Year"=GOA_walleye_pollock[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_walleye_pollock[,'lat'], "Lon"=GOA_walleye_pollock[,'lon'], "Pass"=0)
  }
  if(data_set=="Aleutian_islands_POP"){
    data( AI_pacific_ocean_perch, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=AI_pacific_ocean_perch[,'cpue..kg.km.2.'], "Year"=AI_pacific_ocean_perch[,'year'], "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=AI_pacific_ocean_perch[,'start.latitude'], "Lon"=AI_pacific_ocean_perch[,'start.longitude'], "Pass"=0)
  }
  if( data_set=="GB_spring_haddock"){
    data( georges_bank_haddock_spring, package="FishStatsUtils" )
    Print_Message( "GB_haddock" )
    sampling_data = data.frame( "Catch_KG"=georges_bank_haddock_spring[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_spring[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_spring[,'LATITUDE'], "Lon"=georges_bank_haddock_spring[,'LONGITUDE'])
  }
  if( data_set=="GB_fall_haddock"){
    data( georges_bank_haddock_fall, package="FishStatsUtils" )
    Print_Message( "GB_haddock" )
    sampling_data = data.frame( "Catch_KG"=georges_bank_haddock_fall[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_fall[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_fall[,'LATITUDE'], "Lon"=georges_bank_haddock_fall[,'LONGITUDE'])
  }
  if( data_set=="SAWC_jacopever"){
    data( south_africa_westcoast_jacopever, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=south_africa_westcoast_jacopever[,'HELDAC'], "Year"=south_africa_westcoast_jacopever[,'Year'], "Vessel"="missing", "AreaSwept_km2"=south_africa_westcoast_jacopever[,'area_swept_nm2']*1.852^2, "Lat"=south_africa_westcoast_jacopever[,'cen_lat'], "Lon"=south_africa_westcoast_jacopever[,'cen_long'])
  }
  if( data_set %in% c("Iceland_cod")){
    # WARNING:  This data set has not undergone much evaluation for spatio-temporal analysis
    data( iceland_cod, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=iceland_cod[,'Catch_b'], "Year"=iceland_cod[,'year'], "Vessel"=1, "AreaSwept_km2"=iceland_cod[,'towlength'], "Lat"=iceland_cod[,'lat1'], "Lon"=iceland_cod[,'lon1'])
  }
  if( data_set %in% c("Chatham_rise_hake")){
    data( chatham_rise_hake, package="FishStatsUtils" )
    sampling_data = data.frame( "Catch_KG"=chatham_rise_hake[,'Hake_kg_per_km2'], "Year"=chatham_rise_hake[,'Year'], "Vessel"=1, "AreaSwept_km2"=1, "Lat"=chatham_rise_hake[,'Lat'], "Lon"=chatham_rise_hake[,'Lon'])
  }
  sampling_data = na.omit( sampling_data )

  Return = list("sampling_data"=sampling_data, "Region"=region)
  return(Return)
}
