
#' Format habitat covariate matrix
#'
#' \code{format_covariates} formats an array used by \code{Data_Fn} for habitat (a.k.a. density) covariates
#'
#' @param Lat_e, Latitude for covariate sample e
#' @param Lon_e, Longitude for covariate sample e
#' @param t_e, Time (e.g., year) for covariate sample e
#' @param Cov_ep, matrix of covariates
#' @param Spatial_List, Output from \code{FishStatsUtils::Spatial_Information_Fn}, representing spatial location of knots
#' @param FUN, function used to aggregate observations associated with each knot-time combination
#' @param Year_Set, Set of times \code{t_e} used when generating \code{Cov_xtp}
#' @param na.omit, What to do when some knot-time combination has no observation. Options include \code{"error"} which throw an error, or \code{"time-average"} which fills in the average for other years with observations for each knot

#' @return Tagged list of useful output
#' \describe{
#'   \item{Cov_xtp}{3-dimensional array for use in \code{VAST::Data_Fn}}
#'   \item{var_p}{a matrix summarizing the data-level variance, variance among knots, and lost variance when aggregating from data to knots for each covariate}
#' }

#' @export
format_covariates = function( Lat_e, Lon_e, t_e, Cov_ep, Extrapolation_List, Spatial_List, FUN=mean, Year_Set=min(t_e):max(t_e), na.omit="error" ){

  # Knots in UTM: Spatial_List$loc_x
  # Info for projection:  Extrapolation_List[c('zone','flip_around_dateline')]

  # Backwards compatibility
  if( is.null(Spatial_List$fine_scale) ) Spatial_List$fine_scale = FALSE

  # Determine grids to use
  if( Spatial_List$fine_scale==FALSE ){
    loc_g = Spatial_List$loc_x
  }else{
    if(is.null(Spatial_List)) stop("Problem with `Spatial_List`")
    loc_g = Spatial_List$loc_g
  }

  #
  if( is.vector(Cov_ep)) Cov_ep = matrix(Cov_ep,ncol=1)

  # Step 1:  Project Lat_e and Lon_e to same coordinate system as knots
  if( is.numeric(Extrapolation_List$zone) ){
    loc_e = Convert_LL_to_UTM_Fn( Lon=Lon_e, Lat=Lat_e, zone=Extrapolation_List$zone, flip_around_dateline=Extrapolation_List$flip_around_dateline )                                                         #$
    loc_e = cbind( 'E_km'=loc_e[,'X'], 'N_km'=loc_e[,'Y'])
  }else{
    loc_e = Convert_LL_to_EastNorth_Fn( Lon=Lon_e, Lat=Lat_e, crs=Extrapolation_List$zone )
  }

  # Associate each knot with all covariate measurements that are closest to that knot
  if( na.omit %in% c("error","time-average") ){
    # Step 2: Determine nearest knot for each LatLon_e
    NN = RANN::nn2( data=loc_g[,c('E_km','N_km')], query=loc_e[,c('E_km','N_km')], k=1 )$nn.idx[,1]

    # Step 3: Determine average covariate for each knot and year
    Cov_xtp = NULL
    for(pI in 1:ncol(Cov_ep)){
      Cov_xt = tapply( Cov_ep[,pI], INDEX=list(factor(NN,levels=1:nrow(loc_g)), factor(t_e,levels=Year_Set)), FUN=FUN )
      Cov_xtp = abind::abind( Cov_xtp, Cov_xt, along=3 )
    }

    # Step 4: QAQC
    if( any(is.na(Cov_xtp)) ){
      if( na.omit=="error" ){
        stop("No covariate value for some combination of knot and year")
      }
      if( na.omit=="time-average"){
        Cov_xtp = ifelse( is.na(Cov_xtp), aperm(apply(Cov_xtp,MARGIN=c(1,3),FUN=mean, na.rm=TRUE)%o%rep(1,length(Year_Set)),c(1,3,2)), Cov_xtp )
      }
    }

    # Step 5: Calculate variance lost when aggregating to knots
    KnotVar_p = apply( Cov_xtp, MARGIN=3, FUN=function(vec){var(as.vector(vec))} )
    DataVar_p = apply( Cov_ep, MARGIN=2, FUN=var )
    LostVar_p = 1 - KnotVar_p / DataVar_p
    var_p = matrix( c(KnotVar_p,DataVar_p,LostVar_p), nrow=3, byrow=TRUE)
    dimnames(var_p) = list( c("Variance_among_knots","Variance_among_observations","Proportion_of_variance_lost"), colnames(Cov_ep) )
    for(pI in 1:ncol(Cov_ep)){
      message( "Projecting to knots lost ", formatC(100*LostVar_p[pI],format="f",digits=1), "% of variance for variable ", ifelse(is.null(colnames(Cov_ep)[pI]),pI,colnames(Cov_ep)[pI]))
    }
  }

  # Associate each knot with the X closest measurements in a given year
  if( is.numeric(na.omit) ){
    Cov_xtp = array(NA, dim=c(nrow(loc_g),length(Year_Set),ncol(Cov_ep)), dimnames=list(NULL,Year_Set,colnames(Cov_ep)) )

    for(tI in 1:length(Year_Set)){
      locprime_e = loc_e[which(t_e==Year_Set[tI]),c('E_km','N_km'),drop=FALSE]
      if( nrow(locprime_e) == 0 ) stop("No measurements for year ", Year_Set[tI] )
      NN = RANN::nn2( data=locprime_e[,c('E_km','N_km')], query=loc_g[,c('E_km','N_km')], k=na.omit )$nn.idx
      for(xI in 1:dim(Cov_xtp)[1] ){
      for(pI in 1:dim(Cov_xtp)[3] ){
        Cov_xtp[xI,tI,pI] = FUN( Cov_ep[which(t_e==Year_Set[tI]),pI][NN[xI,]] )
      }}
    }
  }

  # Calculate interpolated value
  Return = list( Cov_xtp=Cov_xtp )
  if( na.omit %in% c("error","time-average") ){
    Return[["var_p"]] = var_p
  }
  return( Return )
}

