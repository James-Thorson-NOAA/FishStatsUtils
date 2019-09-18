


#' @title
#' Plot shifts in range edges
#'
#' @description
#' \code{plot_range_edge} plots range edges
#'
#' @inheritParams plot_biomass_index
#' @param Sdreport TMB output from `TMB::sdreport(Obj)`
#' @param Obj Fitted TMB object from package `VAST`, i.e., output from `fit_model(...)$tmb_list$Obj`
#' @param working_dir Directory for plots
#' @param quantiles vector
#'
#'

#' @export
plot_range_edge = function( Sdreport, Obj, Year_Set=NULL, Years2Include=NULL, strata_names=NULL,
  category_names=NULL, working_dir=paste0(getwd(),"/"), quantiles=c(0.05,0.95), n_samples=100,
  interval_width=1, width=NULL, height=NULL, ...){

  # Unpack
  Report = Obj$report()
  TmbData = Obj$env$data

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")
  if( !("jointPrecision" %in% names(Sdreport))) stop("jointPrecision not present in Sdreport; please re-run with `getJointPrecision=TRUE`")
  if( any(quantiles<0) | any(quantiles>1) ) stop("Please provide `quantiles` between zero and one")
  if( all(TmbData$Z_gm==0) ) stop("Please re-run with 'Options['Calculate_Range']=TRUE' to calculate range edges")

  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    stop("Not implemente")
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    stop("Not implemente")
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }

  # Default inputs
  if( is.null(Year_Set)) Year_Set = 1:TmbData$n_t
  if( is.null(Years2Include) ) Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c
  if( is.null(colnames(TmbData$Z_gm)) ){
    m_labels = paste0("axis",1:ncol(TmbData$Z_gm))
  }else{
    m_labels = colnames(TmbData$Z_gm)
  }

  #### Local function
  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- Matrix::Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }

  # Sample from densities
  u_zr = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Sdreport$jointPrecision, n.sims=n_samples)
  D_gcyr = array( NA, dim=c(dim(Report$D_gcy),n_samples) )
  for( rI in 1:n_samples ){
    if( rI%%floor(n_samples/10) == 1 ) message( "Obtaining sample ", rI, " from predictive distribution for density" )
    D_gcyr[,,,rI] = Obj$report( par=u_zr[,rI] )$D_gcy
  }

  # Calculate quantiles from observed and sampled densities D_gcy
  Z_zm = TmbData$Z_gm
  E_zctm = array(NA, dim=c(length(quantiles),dim(Report$D_gcy)[2:3],ncol(Z_zm)) )
  E_zctmr = array(NA, dim=c(length(quantiles),dim(Report$D_gcy)[2:3],ncol(Z_zm),n_samples) )
  prop_zctm = array(NA, dim=c(dim(Report$D_gcy)[1:3],ncol(Z_zm)) )
  prop_zctmr = array(NA, dim=c(dim(Report$D_gcy)[1:3],ncol(Z_zm),n_samples) )
  for( rI in 0:n_samples ){
  for( mI in 1:ncol(TmbData$Z_gm) ){
    order_g = order(TmbData$Z_gm[,mI], decreasing=FALSE)
    Z_zm[,mI] = Z_zm[order_g,mI]
    if(rI==0) prop_zctm[,,,mI] = apply( Report$D_gcy, MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
    if(rI>=0) prop_zctmr[,,,mI,rI] = apply( D_gcyr[,,,rI,drop=FALSE], MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )

    # Calculate edge
    for( zI in 1:dim(E_zctm)[1] ){
    for( cI in 1:dim(E_zctm)[2] ){
    for( tI in 1:dim(E_zctm)[3] ){
      if(rI==0){
        index_tmp = which.min( (prop_zctm[,cI,tI,mI]-quantiles[zI])^2 )
        E_zctm[zI,cI,tI,mI] = TmbData$Z_gm[order_g[index_tmp],mI]
      }
      if(rI>=1){
        index_tmp = which.min( (prop_zctmr[,cI,tI,mI,rI]-quantiles[zI])^2 )
        E_zctmr[zI,cI,tI,mI,rI] = TmbData$Z_gm[order_g[index_tmp],mI]
      }
    }}}
  }}
  SE_zctm = apply( E_zctmr, MARGIN=1:4, FUN=sd )
  Edge_zctm = abind::abind( "Estimate"=E_zctm, "Std. Error"=SE_zctm, along=5 )
  dimnames(Edge_zctm)[[1]] = paste0("quantile_",quantiles)

  # Plot cumulative distribution
  # matplot( x=Z_zm[,mI], y=prop_zcym[,1,,mI], type="l" )

  # Plot edges
  for( mI in 1:dim(E_zctm)[4] ){
    Index_zct = array(Edge_zctm[,,,mI,'Estimate'],dim(Edge_zctm)[1:3])
    sd_Index_zct = array(Edge_zctm[,,,mI,'Std. Error'],dim(Edge_zctm)[1:3])
    plot_index( Index_ctl=aperm(Index_zct,c(2,3,1)),
      sd_Index_ctl=aperm(sd_Index_zct,c(2,3,1)),
      Year_Set=Year_Set, Years2Include=Years2Include, strata_names=quantiles, category_names=category_names,
      DirName=working_dir, PlotName=paste0("RangeEdge_",m_labels[mI],".png"), Yrange=c(NA,NA),
      interval_width=interval_width, width=width, height=height, xlab="Year", ylab=paste0("Quantiles (",m_labels[mI],")") )
  }


  # Return list of stuff
  Return = list( "Year_Set"=Year_Set, "Edge_zctm"=Edge_zctm )
  return( invisible(Return) )
}

