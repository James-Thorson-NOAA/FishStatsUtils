#' @title
#' Plot quantile residuals on map
#'
#' @description
#' \code{plot_quantile_residuals} shows aggregated quantile residuals at a user-specified spatial resolution
#'
#' This function plots quantile residuals on a map, and is useful to check for persistent spatial patterns
#' in those residuals.  For a check for outliers (i.e., the scale of residuals) please see the separate Q-Q plot.
#' See \code{\link{summary.fit_model}} for more details regarding inputs
#'
#' @inheritParams plot_variable
#' @inheritParams plot_maps
#' @param output from \code{summary.fit_model(x,what="residuals")}
#'
#' @param x Output from \code{\link{fit_model}}
#' @param ... arguments passed to \code{\link{plot_variable}}
#'
#' @export
plot_quantile_residuals <-
function( dharmaRes,
          fit,
          file_name = "quantile_residuals_on_map",
          zlim = NULL,
          n_cells = NULL,
          year_labels = fit$year_labels,
          years_to_plot = fit$years_to_plot,
          ... ){

  # labels
  if( is.null(year_labels) ) year_labels = 1:fit$data_list$n_t
  if( is.null(years_to_plot) ) years_to_plot = 1:fit$data_list$n_t

  # Calculate DHARMa residuals
  if(missing(dharmaRes)){
    dharmaRes = summary( fit, what="residuals", working_dir=NA )
  }

  # Aggregate quantile residuals
  # See Eq. 1 here: https://www.researchgate.net/publication/259073068_Giants'_shoulders_15_years_later_Lessons_challenges_and_guidelines_in_fisheries_meta-analysis
  aggregate_pvalues = function( pvec, na.rm=TRUE ){
    chisq = -2*sum(log(pvec), na.rm=na.rm)
    p = 1 - pchisq( chisq, df=2*length(pvec) )
  }

  # Sanity check that equation holds for p-values
  #pvalues = sapply(dharmaRes$scaledResiduals,FUN=aggregate_pvales)
  #summary(dharmaRes$scaledResiduals - (1-pvalues))

  # Make plot using plot_variable
  PlotDF = cbind( "Lat"=fit$data_frame[,'Lat_i'], "Lon"=fit$data_frame[,'Lon_i'], "x2i"=1:fit$data_list$n_i, "Include"=TRUE)
  Y_gt = matrix(NA, nrow=nrow(fit$data_frame), ncol=fit$data_list$n_t )
  Y_gt[ cbind(1:fit$data_list$n_i,fit$data_list$t_i+1) ] = dharmaRes$scaledResiduals
  Y_gt = Y_gt[,years_to_plot,drop=FALSE]
  col_function = colorRampPalette(colors=c("darkblue","lightblue","white","pink","red"))
  plot_variable( Y_gt = Y_gt,
    map_list = list(PlotDF = PlotDF),
    file_name = file_name,
    fun = aggregate_pvalues,
    col = col_function,
    panel_labels = year_labels[years_to_plot],
    zlim = zlim,
    n_cells = n_cells,
    ... )

  return( NULL )
}
