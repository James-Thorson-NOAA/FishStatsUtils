
#' @export
plot_timeseries = function( x, y, ybounds, fn=lines, bounds_type="whiskers",
  bounds_args=list(), ... ){

  # Local function -- combine two lists
  combine_lists = function( default, input ){
    output = default
    for( i in seq_along(input) ){
      if( names(input)[i] %in% names(default) ){
        output[[names(input)[i]]] = input[[i]]
      }else{
        output = c( output, input[i] )
      }
    }
    return( output )
  }


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
