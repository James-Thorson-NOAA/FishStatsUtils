#' @title
#' Plot semivariance for quantile residuals
#'
#' @description
#' \code{plot_residual_semivariance} calculates the spatial and temporal semivariance from standard-normal quantile residuals
#'
#' This function plots takes quantile residuals, converts them to a standard-normal distribution,
#'      fits the two-dimensional semi-variance in space and time, and plots this.
#'      A well-performing model will have a semivariance for standard-normal-transformed
#'      quantile residuals of 1.0 for all space and time lags.
#'
#' A sill (asymptotic variance for large distance or time) less than one implies
#'      underdispersion, and greater than one implies overdispersion, and both of
#'      these should also be visible in the quantile-quantile plot.
#'
#' A nugget (variance for small distance or time) less than one implies that the
#'      model has some residual spatial or temporal pattern in residuals.
#'
#' @export
plot_residual_semivariance <-
function( fit,
          dharma_raster,
          dharmaRes,
          file_name = "quantile_residuals_semivariance",
          working_dir = paste0(getwd(),"/") ){

  # Transform
  tmp = dharma_raster$Raster_proj
  for(i in seq_along(tmp)) raster::values(tmp[[i]]) = qnorm((1-1/dharmaRes$nSim)*raster::values(tmp[[i]]) + 1/dharmaRes$nSim/2)
  tmp = raster::stack(tmp)

  # Set up
  sp = raster::rasterToPoints(tmp, spatial=TRUE)
  time = as.POSIXct( paste0(fit$year_labels,"-01-01") )
  mydata = data.frame(values = unlist(sp@data) ) # , ID=IDs)
  if(length(time)==1){
    endTime = time+1
  }else{
    endTime = spacetime::delta(time)
  }
  stfdf = spacetime::STFDF(sp=sp, time=time, data=mydata, endTime=endTime)
  residual_semivariance = gstat::variogram(values~1, stfdf, width=20, cutoff = 500, tlags=0:min(5,length(time)-1) )

  # Plot
  ThorsonUtilities::save_fig( file=paste0(run_dir,file_name), width=6, height=5 )
    # gstat:::plot.gstatVariogram
    # gstat:::plot.StVariogram
    semivariance_plot = plot( residual_semivariance,
          at = seq(0, max(residual_semivariance$gamma,c(1.2),na.rm=TRUE), length=23) )#, wireframe=TRUE), plot.numbers=TRUE
    plot(semivariance_plot)
  dev.off()

  #
  return( invisible(residual_semivariance) )
}
