#' @title
#' Plot quantile residuals on map
#'
#' @description
#' \code{plot_quantile_residuals} shows aggregated quantile residuals at a user-specified spatial resolution
#'
#' @inheritParams plot_variable
#'
#' @param x Output from \code{\link{fit_model}}
#' @param ... arguments passed to \code{FisshStatsUtils::plot_variable}
#'
#' @export
plot_quantile_residuals = function( dharmaRes, fit, file_name="quantile_residuals_on_map",
  Year_Set=NULL, Years2Include=NULL, ... ){

  # labels
  if( is.null(Year_Set) ) Year_Set = 1:fit$data_list$n_t
  if( is.null(Years2Include) ) Years2Include = 1:fit$data_list$n_t

  # Calculate DHARMa residuals
  if(missing(dharmaRes)){
    dharmaRes = summary( fit, what="residuals", working_dir=NA )
  }

  # Aggregate quantile residuals
  # See Eq. 1 here: https://www.researchgate.net/publication/259073068_Giants'_shoulders_15_years_later_Lessons_challenges_and_guidelines_in_fisheries_meta-analysis
  aggregate_pvalues = function( pvec, na.rm=TRUE ){
    chisq = -2*sum(log(pvec), na.rm=na.rm)
    p = pchisq( chisq, df=2*length(pvec) )
  }

  # Sanity check that equation holds for p-values
  #pvalues = sapply(dharmaRes$scaledResiduals,FUN=aggregate_pvales)
  #summary(dharmaRes$scaledResiduals - (1-pvalues))

  # Make plot using plot_variable
  PlotDF = cbind( "Lat"=fit$data_frame[,'Lat_i'], "Lon"=fit$data_frame[,'Lon_i'], "x2i"=1:fit$data_list$n_i, "Include"=TRUE)
  Y_gt = matrix(NA, nrow=nrow(fit$data_frame), ncol=fit$data_list$n_t )
  Y_gt[ cbind(1:fit$data_list$n_i,fit$data_list$t_i+1) ] = dharmaRes$scaledResiduals
  Y_gt = Y_gt[,Years2Include,drop=FALSE]
  col_function = colorRampPalette(colors=c("darkblue","lightblue","white","pink","red"))
  plot_variable( Y_gt=Y_gt, map_list=list(PlotDF=PlotDF), file_name=file_name,
    fun=aggregate_pvalues, col=col_function,
    panel_labels=Year_Set[Years2Include], ... )

  return( NULL )
}
