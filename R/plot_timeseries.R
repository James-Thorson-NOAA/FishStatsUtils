
#' Plot timeseries estimates
#'
#' \code{plot_timeseries} plots a timeseries with standard errors shown as whiskers or a shaded polygon
#'
#' @param x numeric-vector of x-coordinates to plot
#' @param y numeric-vector of y-coordinates to plot
#' @param ysd numeric-vector of standard errors in y-coordinates, used in normal approximation to interval
#' @param ybounds alternate specification of y-coordinate whiskers; useful to avoid normal approximation
#' @param fn what plotting function to use; default \code{fn=lines} does not create a new plotting window
#' @param bounds_type which type either \code{"whiskers"} or \code{"shading"}
#' @param bounds_args tagged list of additional arguments to pass to \code{lines} or \code{polygon}
#' @param interval_width width of interval in normal approximation; only used when specifying \code{ysd} without \code{ybounds}
#'
#' @export
plot_timeseries = function( x, y, y_sd, ybounds, fn=lines, bounds_type="whiskers",
  bounds_args=list(), interval_width=1, ... ){

  # fill in missing
  if(missing(ybounds)) ybounds = cbind(y-interval_width*y_sd, y+interval_width*y_sd)

  # Plot lines
  fn( y=y, x=x, ... )

  # Plot shading
  if( bounds_type=="whiskers" ){
    for(t in 1:length(y)){
      lines_args = combine_lists( input=bounds_args, default=list(col="black", lty="solid", lwd=1, x=rep(x[t],2), y=ybounds[t,]) )
      do.call( what=lines, args=lines_args )
    }
  }
  if( bounds_type=="shading" ){
    # Combine entries
    polygon_args = combine_lists( input=bounds_args, default=list(col="black", border=NA, lty="solid", lwd=1, x=c(x,rev(x)), y=c(ybounds[,1],rev(ybounds[,2]))) )
    #polygon(, col=col_bounds, border=border, lty=border_lty)
    do.call( what=polygon, args=polygon_args )
  }
}
