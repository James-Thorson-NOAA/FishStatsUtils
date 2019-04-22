


#' Explore counter-factual scenario
#'
#' \code{Rerun_Fn} re-builds a model while fixing a subset of parameters at zero
#'
#' @param parhat0 parameter values to use by default (presumably maximum-likelihood estimates)
#' @param turnoff_pars character-vector of parameters to turn off
#' @param loc_x location for each knot used to calculate center-of-gravity for counter-factual runs with colnames c('E_km','N_km','Lon','Lat')
#' @param cov_to_turnoff which covariates to turn off, as indicated by their order in \code{TmbData$X_xtp} (which only matters if "gamma1_ctp" or "gamma2_ctp" is in \code{turnoff_pars})
#' @param calculate_COG Boolean whether to calculate COG for each run
#' @param figname name for figure to plot density in counter-factual scenario
#' @inheritParams Build_TMB_Fn
#' @param MapDetails_List output from \code{FishStatsUtils::MapDetails_Fn}
#' @param year_set set of parameters to include
#' @param c_set set of categories to include
#' @param ... additional arguments passed to \code{VAST::Build_TMB_Fn}

#' @return Tagged list
#' \describe{
#'   \item{Report}{Report output for counter-factual run}
#'   \item{NewBuild_List}{Output from \code{VAST::Build_TMB_Fn} using counter-factual parameters}
#' }

#' @examples
#' \dontrun{
#' # Run without GMRF
#'
#' Rerun_Fn(parhat0 = Obj$env$parList(),
#'          turnoff_pars = c("Epsiloninput1_sft", "Epsiloninput2_sft"),
#'          loc_x = Spatial_List$loc_x,
#'          TmbData = TmbData, Version = "VAST_v4_0_0")
#' }

Rerun_Fn = function( parhat0, turnoff_pars, loc_x, cov_to_turnoff=1:dim(parhat0[["gamma2_ctp"]])[3], calculate_COG=TRUE, figname=NULL,
  Map="generate", MapDetails_List=NULL, year_set=1:ncol(parhat0[["beta1_ct"]]), c_set=1:nrow(parhat0[["beta1_ct"]]), ... ){

  # Local function -- calculate center of gravity
  Calc_COG = function( z_x, B_xt ){
    # Set-up
    if( is.vector(z_x) ) z_x = matrix( z_x, ncol=1 )
    P_xt = B_xt / ( rep(1,nrow(B_xt)) %o% apply(B_xt, MARGIN=2, FUN=sum) )

    # Loop through dimensions
    COG_t = NULL
    for( i in 1:ncol(z_x)){
      COG_t = cbind( COG_t, apply( P_xt, MARGIN=2, FUN=weighted.mean, x=z_x[,i]))
    }
    colnames(COG_t) =colnames(z_x)

    return( COG_t )
  }

  # Local function -- plot density maps
  plot_density = function( Report, MapDetails_List, type, year_set, c_set ){
    if("D_xct" %in% names(Report)) D_xcy = Report$D_xct
    if("D_xcy" %in% names(Report)) D_xcy = Report$D_xcy
    if( type=="year" ){
      Dim = c("Nrow"=ceiling(sqrt(length(year_set)))); Dim = c(Dim,"Ncol"=ceiling(length(year_set)/Dim['Nrow']))
      for( i in 1:(length(c_set)-1) ){
        Mat = log(D_xcy[,i,])
        Zlim = range( log(D_xcy[,i,]) )
        FishStatsUtils:::PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], Mat=Mat, zlim=Zlim, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName="", Year_Set=year_set, Rescale=FALSE, Rotate=MapDetails_List[["Rotate"]], Format="", Res=NA, mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,4,0), Cex=MapDetails_List[["Cex"]], textmargin=textmargin, add=FALSE, pch=20, zone=MapDetails_List[["Zone"]] )
        mtext( side=3, outer=TRUE, text=paste0("Size bounds:",c_set[i]," to ",c_set[i+1]), cex=1.5, line=1 )
      }
    }
    if( type=="category" ){
      Dim = c("Nrow"=ceiling(sqrt(length(c_set)-1))); Dim = c(Dim,"Ncol"=ceiling((length(c_set)-1)/Dim['Nrow']))
      for( i in 1:length(year_set) ){
        Mat = log(D_xcy[,,i])
        Zlim = range( log(D_xcy[,,i]) )
        FishStatsUtils:::PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], Mat=Mat, zlim=Zlim, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName="", Year_Set=c_set[-1], Rescale=FALSE, Rotate=MapDetails_List[["Rotate"]], Format="", Res=NA, mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,4,0), Cex=MapDetails_List[["Cex"]], textmargin=textmargin, add=FALSE, pch=20, zone=MapDetails_List[["Zone"]] )
        mtext( side=3, outer=TRUE, text=paste0("Year:",year_set[i]), cex=1.5, line=1 )
      }
    }    #
  }

  # Check for errors
  if( !all(turnoff_pars %in% c("Epsiloninput1_sft","Epsiloninput2_sft","beta1_ct","beta2_ct","gamma1_ctp","gamma2_ctp")) ) stop("Some parameters in turnoff_pars are not as expected")
  if( !all(turnoff_pars %in% names(parhat0)) ) stop("Some parameters in turnoff_pars are not in parhat0")

  # Eliminate spatio-temporal variation and differences in age-composition
  parhat_mod = parhat0
  if( "Epsiloninput1_sft" %in% turnoff_pars ) parhat_mod[["Epsiloninput1_sft"]][] = 0
  if( "Epsiloninput2_sft" %in% turnoff_pars ) parhat_mod[["Epsiloninput2_sft"]][] = 0
  if( "beta1_ct" %in% turnoff_pars ) parhat_mod[["beta1_ct"]][] = outer( rowMeans(ParHat[["beta1_ct"]]), rep(1,ncol(ParHat[["beta1_ct"]])) )
  if( "beta2_ct" %in% turnoff_pars ) parhat_mod[["beta2_ct"]][] = outer( rowMeans(ParHat[["beta2_ct"]]), rep(1,ncol(ParHat[["beta2_ct"]])) )
  which_cov_to_turnoff = intersect(cov_to_turnoff, 1:dim(parhat_mod[["gamma1_ctp"]])[3])
  if( "gamma1_ctp" %in% turnoff_pars ) parhat_mod[["gamma1_ctp"]][,,which_cov_to_turnoff] = 0
  if( "gamma2_ctp" %in% turnoff_pars ) parhat_mod[["gamma2_ctp"]][,,which_cov_to_turnoff] = 0

  # Rebuild
  NewBuild_List = Build_TMB_Fn( "Parameters"=parhat_mod, Map=Map, ... )
  Obj = NewBuild_List[["Obj"]]
  Report = Obj$report( )
  Return = list("Report"=Report, "NewBuild_List"=NewBuild_List)

  # calculate COG
  if( calculate_COG==TRUE ){
    # Figure out column names to use
    ColNames = c('E_km','N_km')
    if( all(c('Lon','Lat') %in% colnames(loc_x)) ) ColNames = c(ColNames,c('Lon','Lat'))
    # Total
    B_xt = apply(Report$Index_xcyl[,,,1,drop=FALSE], MARGIN=c(1,3), FUN=sum )
    Return[["COG_t"]] = Calc_COG( z_x=loc_x[,ColNames], B_xt=B_xt )
    # By category
    Return[["COG_ct"]] = NULL
    for( cI in 1:dim(Report$Index_xcyl)[2] ){
      Return[["COG_ct"]] = abind::abind( Return[["COG_ct"]], Calc_COG(z_x=loc_x[,ColNames], B_xt=Report$Index_xcyl[,cI,,1]), along=3)
    }
    Return[["COG_ct"]] = aperm(Return[["COG_ct"]], c(3,1,2))
  }

  # Plot
  #if( !is.null(figname) & !is.null(MapDetails_List) ){
  #  message( "Starting plot for counter-factual run" )
  #  ThorsonUtilities::save_fig( file=figname, width=8, height=12, type="pdf", onefile=TRUE )
  #    plot_density( Report=Report, MapDetails_List=MapDetails_List, type="year", year_set=Year_Set, c_set=c_set )
  #  dev.off()
  #}

  # Return
  return( Return )
}

#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams plot_overdispersion
#' @inheritParams VAST::make_data
#' @param covhat estimated covariance used for calculating coherence

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }

calc_coherence = function( Report, Data, covhat=NULL, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  ##################
  # Synchrony
  ##################
  n_z = nrow(yearbounds_zz)

  # Category-specific density ("D"), mean and variance
  varD_xcz = array( NA, dim=c(Data$n_x,Data$n_c,n_z) )

  # Community density ("C") summing across categories
  C_xt = array(0, dim=c(Data$n_x,Data$n_t))
  varC_xz = array(NA, dim=c(Data$n_x,n_z))

  # Area-weighted total biomass ("B") for each category
  B_ct = array(0, dim=c(Data$n_c,Data$n_t))
  B_t = rep(0, Data$n_t)

  # Proportion of total biomass ("P") for each station
  P_xz = array(NA, dim=c(Data$n_x,n_z))

  # Synchrony indices
  phi_xz = array(NA, dim=c(Data$n_x,n_z))
  phi_z = rep(0, n_z)

  # Temporary variables
  temp_xt = array(NA, dim=c(Data$n_x,Data$n_t))
  temp_c = rep(NA, Data$n_c)

  # Derived quantities
  for( tI in 1:Data$n_t ){
    for( cI in 1:Data$n_c ){
      for( xI in 1:Data$n_x ){
        C_xt[xI,tI] = C_xt[xI,tI] + Report$D_xct[xI,cI,tI]
        B_ct[cI,tI] = B_ct[cI,tI] + Report$D_xct[xI,cI,tI] * Data$a_xl[xI,1]
      }
      B_t[tI] = B_t[tI] + B_ct[cI,tI]
    }
  }

  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z){
  for( xI in 1:Data$n_x){
    # Variance for each category
    for( cI in 1:Data$n_c ){
      for( tI in 1:Data$n_t ) temp_xt[xI,tI] = Report$D_xct[xI,cI,tI]
      varD_xcz[xI,cI,zI] = var(temp_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    }
    # Variance for combined across categories
    varC_xz[xI,zI] = var(C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    # Proportion in each category
    P_xz[xI,zI] = Data$a_xl[xI,1] * mean( C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] / B_t[yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] )
    # Synchrony index
    for( cI in 1:Data$n_c ) temp_c[cI] = varD_xcz[xI,cI,zI]
    phi_xz[xI,zI] = varC_xz[xI,zI] / sum(sqrt(temp_c))^2
    phi_z[zI] = phi_z[zI] + (phi_xz[xI,zI] * P_xz[xI,zI])
  }}

  ##################
  # Coherence
  ##################

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_xz"=phi_xz, "phi_z"=phi_z, "psi"=psi, "P_xz"=P_xz, "L_c"=L_c)
  return( Return )
}


#' Calculate synchrony among species or locations
#'
#' \code{calc_synchrony} calculates synchrony
#'
#' @inheritParams plot_overdispersion
#' @inheritParams VAST::make_data

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#' }

calc_synchrony = function( Report, Data, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  # Index lengths
  n_z = nrow(yearbounds_zz)
  n_t = Data$n_t
  n_c = Data$n_c
  n_x = Data$n_x
  n_y = nrow(Data$t_yz)

  # Extract elements
  pow = function(a,b) a^b
  D_xcy = Report$D_xcy
  a_xl = Data$a_xl

  # Density ("D") or area-expanded total biomass ("B") for each category (use B when summing across sites)
  D_xy = matrix( 0, nrow=n_x, ncol=n_y );
  B_cy = matrix( 0, nrow=n_c, ncol=n_y );
  B_y = rep( 0, n_y );
  # Mean
  meanD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  meanD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  meanB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  meanB_z = rep( 0, n_z );
  # Sample variance in category-specific density ("D") and biomass ("B")
  varD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  varD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  varB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  varB_z = rep( 0, n_z );
  varB_xbar_z = rep( 0, n_z );
  varB_cbar_z = rep( 0, n_z );
  maxsdD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  maxsdB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  maxsdB_z = rep( 0, n_z );
  # Proportion of total biomass ("P") for each location or each category
  propB_xz = matrix( 0, nrow=n_x, ncol=n_z );
  propB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  # Synchrony indices
  phi_xz = matrix( 0, nrow=n_x, ncol=n_z );
  phi_cz = matrix( 0, nrow=n_c, ncol=n_z );
  phi_xbar_z = rep( 0, n_z );
  phi_cbar_z = rep( 0, n_z );
  phi_z = rep( 0, n_z );

  # Calculate total biomass for different categories
  for( yI in 1:n_y ){
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        D_xy[xI,yI] = D_xy[xI,yI] + D_xcy[xI,cI,yI];
        B_cy[cI,yI] = B_cy[cI,yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
        B_y[yI] = B_y[yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
      }
    }
  }
  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z ){
    for( xI in 1:n_x ){
      # Variance for biomass in each category, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( cI in 1:n_c ){
        temp_mean = 0;
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xcz[xI,cI,zI] = meanD_xcz[xI,cI,zI] + D_xcy[xI,cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ){
          varD_xcz[xI,cI,zI] = varD_xcz[xI,cI,zI] + pow(D_xcy[xI,cI,yI]-meanD_xcz[xI,cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
        }
      }
      # Variance for combined biomass across categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xz[xI,zI] = meanD_xz[xI,zI] + D_xy[xI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varD_xz[xI,zI] = varD_xz[xI,zI] + pow(D_xy[xI,yI]-meanD_xz[xI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    for( cI in 1:n_c ){
      # Variance for combined biomass across sites, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_cz[cI,zI] = meanB_cz[cI,zI] + B_cy[cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varB_cz[cI,zI] = varB_cz[cI,zI] + pow(B_cy[cI,yI]-meanB_cz[cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    # Variance for combined biomass across sites and categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_z[zI] = meanB_z[zI] + B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
      varB_z[zI] = varB_z[zI] + pow(B_y[yI]-meanB_z[zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
    }
    # Proportion in each site
    for( xI in 1:n_x ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_xz[xI,zI] = propB_xz[xI,zI] + a_xl[xI,1] * D_xy[xI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Proportion in each category
    for( cI in 1:n_c ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_cz[cI,zI] = propB_cz[cI,zI] + B_cy[cI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Species-buffering index (calculate in Density so that areas with zero area are OK)
    for( xI in 1:n_x ){
      for( cI in 1:n_c ){
        maxsdD_xz[xI,zI] = maxsdD_xz[xI,zI] + pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_xz[xI,zI] = varD_xz[xI,zI] / pow( maxsdD_xz[xI,zI], 2);
      varB_xbar_z[zI] = varB_xbar_z[zI] + pow(a_xl[xI,1],2) * varD_xz[xI,zI] * propB_xz[xI,zI];
      phi_xbar_z[zI] = phi_xbar_z[zI] + phi_xz[xI,zI] * propB_xz[xI,zI];
    }
    # Spatial-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_cz[cI,zI] = maxsdB_cz[cI,zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_cz[cI,zI] = varB_cz[cI,zI] / pow( maxsdB_cz[cI,zI], 2);
      varB_cbar_z[zI] = varB_cbar_z[zI] + varB_cz[cI,zI] * propB_cz[cI,zI];
      phi_cbar_z[zI] = phi_cbar_z[zI] + phi_cz[cI,zI] * propB_cz[cI,zI];
    }
    # Spatial and species-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_z[zI] = maxsdB_z[zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
    }
    phi_z[zI] = varB_z[zI] / pow( maxsdB_z[zI], 2);
  }

  # Return stuff
  Return = list( "phi_z"=phi_z, "phi_xz"=phi_xz, "phi_cz"=phi_cz, "phi_xbar_z"=phi_xz, "phi_cbar_z"=phi_cz, "B_y"=B_y, "B_cy"=B_cy, "D_xy"=D_xy, "D_xcy"=D_xcy)
  return( Return )
}


#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams plot_overdispersion
#' @inheritParams VAST::make_data
#' @param covhat estimated covariance used for calculating coherence

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }

Coherence = function( Report, Data, covhat=NULL, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  ##################
  # Synchrony
  ##################
  n_z = nrow(yearbounds_zz)

  # Category-specific density ("D"), mean and variance
  varD_xcz = array( NA, dim=c(Data$n_x,Data$n_c,n_z) )

  # Community density ("C") summing across categories
  C_xt = array(0, dim=c(Data$n_x,Data$n_t))
  varC_xz = array(NA, dim=c(Data$n_x,n_z))

  # Area-weighted total biomass ("B") for each category
  B_ct = array(0, dim=c(Data$n_c,Data$n_t))
  B_t = rep(0, Data$n_t)

  # Proportion of total biomass ("P") for each station
  P_xz = array(NA, dim=c(Data$n_x,n_z))

  # Synchrony indices
  phi_xz = array(NA, dim=c(Data$n_x,n_z))
  phi_z = rep(0, n_z)

  # Temporary variables
  temp_xt = array(NA, dim=c(Data$n_x,Data$n_t))
  temp_c = rep(NA, Data$n_c)

  # Derived quantities
  for( tI in 1:Data$n_t ){
    for( cI in 1:Data$n_c ){
      for( xI in 1:Data$n_x ){
        C_xt[xI,tI] = C_xt[xI,tI] + Report$D_xct[xI,cI,tI]
        B_ct[cI,tI] = B_ct[cI,tI] + Report$D_xct[xI,cI,tI] * Data$a_xl[xI,1]
      }
      B_t[tI] = B_t[tI] + B_ct[cI,tI]
    }
  }

  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z){
  for( xI in 1:Data$n_x){
    # Variance for each category
    for( cI in 1:Data$n_c ){
      for( tI in 1:Data$n_t ) temp_xt[xI,tI] = Report$D_xct[xI,cI,tI]
      varD_xcz[xI,cI,zI] = var(temp_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    }
    # Variance for combined across categories
    varC_xz[xI,zI] = var(C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    # Proportion in each category
    P_xz[xI,zI] = Data$a_xl[xI,1] * mean( C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] / B_t[yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] )
    # Synchrony index
    for( cI in 1:Data$n_c ) temp_c[cI] = varD_xcz[xI,cI,zI]
    phi_xz[xI,zI] = varC_xz[xI,zI] / sum(sqrt(temp_c))^2
    phi_z[zI] = phi_z[zI] + (phi_xz[xI,zI] * P_xz[xI,zI])
  }}

  ##################
  # Coherence
  ##################

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_xz"=phi_xz, "phi_z"=phi_z, "psi"=psi, "P_xz"=P_xz, "L_c"=L_c)
  return( Return )
}


#' K-fold crossvalidation
#'
#' \code{Crossvalidate_Fn} runs a k-fold crossvalidation analysis on a fitted model object
#'
#' @param record_dir, a directory with writing privilege where the runs are stored
#' @param parhat, a tagged list of parameters from the fitted model, e.g., as generated by \code{parhat <- obj$env$parList()} after covergence
#' @param original_data, a tagged list of data for the VAST model
#' @param group_i, a vector of positive integers, indicating the k-fold group for each observation (default=NULL, which generates a new group_i with even probability)
#' @param kfold, the number of crossvalidation batches used (default=10)
#' @param skip_finished boolean specifying whether to rerun (skip_finished==FALSE) or skip (skip_finished==TRUE) previously completed runs (Default=FALSE)
#' @param skip_finished boolean specifying whether to rerun (skip_finished==FALSE) or skip (skip_finished==TRUE) previously completed runs (Default=FALSE)
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param ... Additional arguments to pass to \code{VAST::Build_TMB_Fn}

#' @return Results a matrix with total predictive negative log-likelihood for each crossvalidation partition, and number of crossvalidation samples for that partition

Crossvalidate_Fn = function(record_dir, parhat, original_data, group_i=NULL, kfold=10, newtonsteps=1, skip_finished=FALSE, ... ){
  # Lump observations into groups
  if( is.null(group_i) || length(group_i)!=original_data$n_i ){
    message( "Generating group_i" )
    Group_i = sample( x=1:kfold, size=original_data$n_i, replace=TRUE )
  }else{
    message( "Using input group_i" )
    Group_i = group_i
  }
  save( Group_i, file=paste0(record_dir,"Group_i.RData"))

  # Results
  Results = array(NA, dim=c(kfold,3), dimnames=list(paste0("k=",1:kfold),c("predictive_negloglike","number_of_crossvalidation_samples","predictive_nll_per_sample")) )

  # Loop through
  for(i in 1:kfold){
    # Directory
    CrossvalidationDir = paste0(record_dir,"k=",i,"/")
    dir.create( CrossvalidationDir )

    # Run
    if( "Save.RData"%in%list.files(CrossvalidationDir) & skip_finished==TRUE ){
      load(file=paste0(CrossvalidationDir,"Save.RData"))
    }else{
      # Modify data
      Data = original_data
      Data$PredTF_i = ifelse(Group_i==i,1,0)

      # Build new one
      TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, ...)#, "Random"=NULL)
      #TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, "RunDir"=TmbDir, "Version"=Version, "loc_x"=loc_x, "RhoConfig"=RhoConfig, "TmbDir"=TmbDir, "Use_REML"=Use_REML) #, "Map"=Save$Map, "Random"=NULL)

      # Extract objects
      Obj = TmbList[["Obj"]]
      TmbList[["Upper"]][grep("logkappa",names(TmbList[["Upper"]]))] = Inf

      # Run model
      Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], newtonsteps=newtonsteps )  # , rel.tol=1e-20

      # Reports
      Report = Obj$report( Obj$env$last.par.best )
      ParHat = Obj$env$parList(Opt$par)

      # Save stuff
      Save = list("Opt"=Opt, "Report"=Report, "ParHat"=ParHat, "TmbData"=Data, "Map"=TmbList$Map)
      save(Save, file=paste0(CrossvalidationDir,"Save.RData"))
      capture.output( Opt, file=paste0(CrossvalidationDir,"Opt.txt"))
    }

    # Load results
    Results[i,c("predictive_negloglike","number_of_crossvalidation_samples")] = c( Save$Report$pred_jnll, sum(Save$TmbData$PredTF_i))
    Results[i,"predictive_nll_per_sample"] = Results[i,"predictive_negloglike"] / Results[i,"number_of_crossvalidation_samples"]
  }

  # Return output
  return( Results)
}


plot_eigen = function( Cov, which2plot=1:min(3,ncol(Cov)), names=1:nrow(Cov), digits=2, las=2, add=FALSE, ... ){
  Eigen = eigen(Cov)
  rownames( Eigen$vectors ) = names
  if(add==FALSE) par( mfrow=par()$mfrow, oma=par()$oma, mar=par()$mar, tck=par()$tck )
  for( cI in which2plot ){
    plot( 1, type="n", xlim=c(0.5,ncol(Cov)+0.5), ylim=c(-1,1.2), xlab="", ylab="", xaxt="n", ... )
    for( rI in 1:nrow(Cov)){
      graphics::lines(y=c(0,Eigen$vectors[rI,cI])*sign(Eigen$vectors[which.max(abs(Eigen$vectors[,cI])),cI]), x=rep(rI,2), lwd=5)
    }
    graphics::abline( h=0, lty="dotted" )
    graphics::legend( "top", bty="n", legend=paste0("Eigenvalue = ",formatC(Eigen$values[cI],format="f",digits=2)) )
    axis(1, at=1:ncol(Cov), labels=names, las=las  )
  }
  return( invisible(Eigen) )
}



Plot_Ellipse_Fn = function( filename, Report, yearset=NULL, years2include=NULL, MappingDetails, MapSizeRatio=c("Height(in)"=4,"Width(in)"=4), Xlim, Ylim, ZinUTM=TRUE, zone=NA, ncol_legend=2, Format="png", ... ){
  # Only run if necessary outputs are available
  if( all(c("mean_Z_tl","cov_Z_tl") %in% names(Report)) ){
    # Decide on years
    if( is.null(years2include) ) years2include = 1:nrow(Report$mean_Z_tl)

    # Plot ellipse -- Lat/Lon OR Lat/Depth
    if(Format=="png") png( file=filename, width=MapSizeRatio['Width(in)'], height=MapSizeRatio['Height(in)'], res=200, units="in")
      # Settings
      ColBorder = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red","darkred"))(nrow(Report$mean_Z_tl)) # c("blue","red")# c("darkblue","blue","lightblue","lightgreen","yellow","orange","red","darkred")
      if( ZinUTM==FALSE ){
        maps::map(MappingDetails[[1]], MappingDetails[[2]], ylim=mean(Ylim)+1*c(-0.5,0.5)*diff(Ylim), xlim=mean(Xlim)+1*c(-0.5,0.5)*diff(Xlim), fill=TRUE, ...) # , orientation=c(mean(y.lim),mean(x.lim),15)
      }else{
        FishStatsUtils:::Plot_States_in_UTM_Fn( MappingDetails=MappingDetails, fillcol=NA, xlim=Xlim, ylim=Ylim, zone=zone, ... )
      }
      # Plot ellipsoids
      for(t in years2include){
        Shape = mixtools::ellipse(mu=Report$mean_Z_tl[t,1:2], sigma=Report$cov_Z_tl[t,1:2,1:2], alpha=1-(1-pnorm(1))*2, newplot=FALSE, draw=FALSE) #, type="l")
        polygon( Shape, border=ColBorder[t], lwd=3 )    # , col=Col[t]
      }
      if( !is.null(yearset)) legend( "top", fill=ColBorder[years2include], bty="n", legend=yearset[years2include], ncol=ncol_legend)
      #for(t in 2:length(Col)) arrows( x0=mean_Z_tl[t-1,1], y0=mean_Z_tl[t-1,2], x1=mean_Z_tl[t,1], y1=mean_Z_tl[t,2], length=0.1 )
    if(Format=="png") dev.off()
  }else{
    message("mean_Z_tl or cov_Z_tl not found in Report")
  }
}

plotSites <-
function( mapObj, x, y, col, cex=2, ... ){
  # Load package
  require(maptools)

  # Plot map
  if( is.list(mapObj) & length(mapObj)>0 ){
    plot(mapObj, type="l", xlim=mapObj[["bbox"]]['x',], ylim=mapObj[["bbox"]]['y',], ...)
  }else{
    plot(1, type="n", xlim=attr(mapObj,"bbox")['x',], ylim=attr(mapObj,"bbox")['y',], ...)
    if( length(mapObj)>1 ) suppressWarnings( plot(mapObj, add=TRUE) )
  }

  #plot site locations
  points(x=x, y=y, pch=20, cex=cex, col=col)
}

Timeseries_Fn <-
function(Report, FileName, Year_Set, ControlList=list("Width"=4*3, "Height"=2*3, "Res"=200, "Units"='in')){
  # Extract stuff from report
  D_it = Report$D_xt[NN_Extrap$nn.idx,]
  R1_it = Report$R1_xt[NN_Extrap$nn.idx,]
  R2_it = Report$R2_xt[NN_Extrap$nn.idx,]
  D_xt = Report$D_xt
  # Plot time series
  png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
    par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(3,3,2,0), mgp=c(1.5,0.5,0), tck=-0.02)
    # Entropy
    EntropyFn = function(Vec){ sum(Vec * log(Vec+1e-250)/log(length(Vec)) ) }
    Entropy_t = apply( Report$D_xt, MARGIN=2, FUN=EntropyFn )
    plot( x=Year_Set, y=Entropy_t, type="l", main="Entropy")
    # Variance
    Var_t = apply( Report$D_xt, MARGIN=2, FUN=var )
    plot( x=Year_Set, y=Var_t, type="l", main="Variance", ylim=c(0,max(Var_t)) )
    # CV
    CV_t = apply( Report$D_xt, MARGIN=2, FUN=function(Vec){ sd(Vec)/mean(Vec) } )
    plot( x=Year_Set, y=CV_t, type="l", main="CV", ylim=c(0,max(CV_t)) )
    # Occupancy probability
    Occup_t = colMeans( R1_it )
    plot( x=Year_Set, y=Occup_t, type="l", main="Occup_t", ylim=c(0,1) )
    # Density
    CondDens_t = colMeans( R2_it )
    plot( x=Year_Set, y=CondDens_t, type="l", main="CondDens_t", ylim=c(0,max(CondDens_t)) )
    # Density
    Index_t = colMeans( D_it )
    plot( x=Year_Set, y=Index_t, type="l", main="Index_t", ylim=c(0,max(Index_t)) )
    # correlation between occupancy and abundance
    Cor_t = sapply(1:length(Year_Set), FUN=function(Num){cor(R1_it[,Num],R2_it[,Num],method="spearman")})
    plot( x=Year_Set, y=Cor_t, type="l", main="Cor_t", ylim=c(-1,1) ); abline(h=0, lty="dotted")
  dev.off()
}


#' @title
#' Plot vessel effects
#'
#' @description
#' \code{Vessel_Fn} plots estimated vessel effects for model components
#'
#' @param TmbData Formatted data inputs, from `VAST::Data_Fn(...)`
#' @param Sdreport TMB output from `TMB::sdreport(Obj)`
#' @param FileName_VYplot Full filename (including directory) for center-of-gravity plot
#'
#' @return Return Tagged list of output
#'

Vessel_Fn <-
function( TmbData, Sdreport, FileName_VYplot=NULL ){
  Summary = TMB:::summary.sdreport(Sdreport)
  Return = NULL

  if( !any(c("nu1_vt","nu2_vt") %in% rownames(Summary)) ){
    message( "Not plotting vessel effects because none are present" )
  }else{
    message( "\nPlotting vessel effects..." )
    nu_vt = array(NA, dim=c(TmbData$n_v, TmbData$n_t, 2, 2))
    for(vI in 1:TmbData$n_v){
    for(tI in 1:TmbData$n_t){
      Num = (vI-1)*TmbData$n_t + tI
      nu_vt[vI,tI,1,] = Summary[which(rownames(Summary)=="nu1_vt")[Num],]
      nu_vt[vI,tI,2,] = Summary[which(rownames(Summary)=="nu2_vt")[Num],]
    }}

    # Make plot
    if( !is.null(FileName_VYplot) ) jpeg(FileName_VYplot, width=1.5*TmbData$n_t,height=5,res=200,units="in")
      par(mfrow=c(2,TmbData$n_t), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(0,3,0,0))
      for(eI in 1:2){
      for(tI in 1:TmbData$n_t){
        plot(x=1:TmbData$n_v, y=1:TmbData$n_v, type="n", ylim=range( c(nu_vt[,,eI,1]+nu_vt[,,eI,2],nu_vt[,,eI,1]-nu_vt[,,eI,2]) ), xlab="Vessel", ylab="Effect", main=TmbData$Year_Set[tI])
        if(tI==1) mtext( side=2, outer=FALSE, line=2, text=c("Presence/absence","Positive catch rate")[eI])
        for(vI in 1:TmbData$n_v){
          points( x=vI, y=nu_vt[vI,tI,eI,1])
          lines( x=rep(vI,2), y=c(nu_vt[vI,tI,eI,1]-nu_vt[vI,tI,eI,2],nu_vt[vI,tI,eI,1]+nu_vt[vI,tI,eI,2]))
        }
      }}
    if( !is.null(FileName_VYplot) ) dev.off()

    # Add stuff
    Return = c(Return, list("nu_vt"=nu_vt) )
  }

  # Return stuff
  return( Return )
}
