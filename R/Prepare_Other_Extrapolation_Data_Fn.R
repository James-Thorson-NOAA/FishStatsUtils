#' @export
Prepare_Other_Extrapolation_Data_Fn <-
function( strata.limits, observations_LL, projargs=NA, grid_dim_km=c(2,2), maximum_distance_from_sample=NULL,
          zone=NA, grid_in_UTM=TRUE, grid_dim_LL=c(0.1,0.1), flip_around_dateline=FALSE, ... ){

  # Local function
  rename_columns = function( DF, origname=colnames(DF), newname ){
    DF_new = DF
    for(i in 1:length(origname)){
      Match = match( origname[i], colnames(DF_new) )
      if(length(Match)==1) colnames(DF_new)[Match] = newname[i]
    }
    return(DF_new)
  }

  if( grid_in_UTM==TRUE ){
    # Fill in missing inputs
    if( is.null(maximum_distance_from_sample) ) maximum_distance_from_sample = sqrt((grid_dim_km[1]/2)^2+(grid_dim_km[2]/2)^2)

    # Get range
    #TmpUTM = Convert_LL_to_UTM_Fn( Lon=observations_LL[,'Lon'], Lat=observations_LL[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
    TmpUTM = project_coordinates( Lon=observations_LL[,'Lon'], Lat=observations_LL[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
    E_lim = mean(range(TmpUTM[,'E_km'])) + c(-0.6,0.6)*diff(range(TmpUTM[,'E_km']))
    N_lim = mean(range(TmpUTM[,'N_km'])) + c(-0.6,0.6)*diff(range(TmpUTM[,'N_km']))

    # Detect Northern or Southern hemisphere
    #NorthernTF = all( observations_LL[,'Lat']>0 )
    #if( any(observations_LL[,'Lat']>0) & any(observations_LL[,'Lat']<0) ) warning( "PBSmapping doesn't work with observations in both Northern and Southern hemisphere" )

    # Make grid
    E_seq = seq(E_lim[1], E_lim[2], by=grid_dim_km[1])
    N_seq = seq(N_lim[1], N_lim[2], by=grid_dim_km[2])
    if( length(E_seq)*length(N_seq) > 1e5 ){
      warning("Please increase `grid_dim_km` to avoid memory issues caused by creating a large extrapolation-grid")
    }
    if( length(E_seq)*length(N_seq) > 1e7 ){
      stop("Please increase `grid_dim_km` to avoid memory-associated crashes caused by creating a huge extrapolation-grid")
    }
    Data_Extrap = expand.grid( "E_km"=E_seq, "N_km"=N_seq, "Area_km2"=prod(grid_dim_km) )

    # Add LL
    #TmpUTM = rename_columns( Data_Extrap[,c("E_km","N_km")], newname=c("X","Y"))
    #attr(TmpUTM, "projection") = "UTM"
    #attr(TmpUTM, "zone") = attr(TmpUTM,"zone")
    #TmpLL = PBSmapping::convUL(TmpUTM, southern=!NorthernTF )
    #if( flip_around_dateline==TRUE ){
    #  TmpLL[,1] = TmpLL[,1] - 180
    #}
    TmpLL = project_coordinates( Lon=Data_Extrap[,"E_km"], Lat=Data_Extrap[,"N_km"], projargs=attr(TmpUTM,"origargs"),
      origargs=attr(TmpUTM,"projargs"), zone=zone, flip_around_dateline=flip_around_dateline)
    Data_Extrap = cbind( Data_Extrap, rename_columns(TmpLL,newname=c("Lon","Lat")) )

    # Restrict to grid locations near samples
    NN_Extrap = RANN::nn2( query=Data_Extrap[,c("E_km","N_km")], data=TmpUTM[,c("E_km","N_km")], k=1)
    Data_Extrap = cbind( Data_Extrap, "Include"=ifelse(NN_Extrap$nn.dists<maximum_distance_from_sample,1,0))

    # Survey areas
    Area_km2_x = Data_Extrap[,'Area_km2'] * Data_Extrap[,'Include']
  }else{
    # Fill in missing inputs
    if( is.null(maximum_distance_from_sample) ) maximum_distance_from_sample = sqrt((grid_dim_LL[1]/2)^2+(grid_dim_LL[2]/2)^2)

    # Get range
    Lat_lim = mean(range(observations_LL[,'Lat'])) + c(-0.6,0.6)*diff(range(observations_LL[,'Lat']))
    if( flip_around_dateline==FALSE ) Lon_lim = mean(range(observations_LL[,'Lon'])) + c(-0.6,0.6)*diff(range(observations_LL[,'Lon']))
    if( flip_around_dateline==TRUE ){
      Tmp_Lon = 180 + ifelse( observations_LL[,'Lon']>0, observations_LL[,'Lon']-360, observations_LL[,'Lon'])
      Lon_lim = mean(range(Tmp_Lon)) + c(-0.6,0.6)*diff(range(Tmp_Lon))
    }

    # Make grid
    Data_Extrap = expand.grid( "Lon"=seq(Lon_lim[1],Lon_lim[2],by=grid_dim_LL[1]), "Lat"=seq(Lat_lim[1],Lat_lim[2],by=grid_dim_LL[2]), "Area_km2"=110^2*prod(grid_dim_LL) )
    if( flip_around_dateline==TRUE ){
      Data_Extrap[,'Lon'] = ifelse( Data_Extrap[,'Lon']<0, Data_Extrap[,'Lon']+360, Data_Extrap[,'Lon'] ) - 180
    }

    # Add empty UTM
    Data_Extrap = cbind( Data_Extrap, "E_km"=NA, "N_km"=NA, "Area_km2"=prod(grid_dim_km) )
    TmpUTM = NA
    attr(TmpUTM, "zone") = NA

    # Restrict to grid locations near samples
    NN_Extrap = RANN::nn2( query=Data_Extrap[,c("Lat","Lon")], data=observations_LL[,c('Lat','Lon')], k=1)
    Data_Extrap = cbind( Data_Extrap, "Include"=ifelse(NN_Extrap$nn.dists<maximum_distance_from_sample,1,0))

    # Survey areas
    Area_km2_x = Data_Extrap[,'Area_km2'] * Data_Extrap[,'Include']
  }

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(TmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}
