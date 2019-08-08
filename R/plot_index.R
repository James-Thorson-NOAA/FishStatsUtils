
#' Plot an index
#'
#' \code{plot_index} takes output from a VAST run and plots a panel figure of time-series estimates
#'
#' @param Index_ctl A matrix or array of time-series estimates for multiple categories \code{c}, years \code{t}, and strata \code{l}
#' @param sd_Index_ctl A matrix or array of variances for each estimate
#' @inheritParams plot_biomass_index
#' @inheritParams plot_lines
#' @param Yrange lower and upper bound for left-hand y-axis, corresponding to input \code{Index_ctl} (use \code{Yrange[1]=NA} and/or \code{Yrange[2]=NA} for using the lower and upper bound of estimate intervals)
#' @param Y2range lower and upper bound for right-hand y-axis, corresponding to input \code{SampleSize_ctz} (see Yrange for more info)
#' @param SampleSize_ctz optional array of sample sizes for each category and year to be plotted on each panel
#' @param plot_args additional arguments to pass to \code{plot}
#' @param plot_lines_args additional arguments to pass to \code{plot_lines}
#' @param ... list of settings to pass to \code{par} when making plot
#'
#' @return NULL
#'
#' @export
plot_index = function( Index_ctl, sd_Index_ctl=array(0,dim(Index_ctl)), Year_Set=NULL, Years2Include=NULL,
  strata_names=NULL, category_names=NULL, scale="uniform",
  plot_legend=NULL, DirName=paste0(getwd(),"/"), PlotName="Index.png",
  interval_width=1, width=NULL, height=NULL, xlab="Year", ylab="Index", bounds_type="whiskers", col=NULL,
  col_bounds=NULL, Yrange=c(0,NA), type="b", plot_lines_args=list(), plot_args=list(),
  SampleSize_ctz=NULL, Y2range=c(0,NA), y2lab="", ... ){

  # Change inputs
  if( length(dim(Index_ctl))==length(dim(sd_Index_ctl)) ){
    if( length(dim(Index_ctl))==2 ){
      Index_ctl = Index_ctl %o% rep(1,1)
      sd_Index_ctl = sd_Index_ctl %o% rep(1,1)
    }
  }else{
    stop("Mismatch in dimensions for `Index_ctl` and `sd_Index_ctl` in `plot_index`")
  }
  n_categories = dim(Index_ctl)[1]
  n_years = dim(Index_ctl)[2]
  n_strata = dim(Index_ctl)[3]
  mfrow = c( ceiling(sqrt(n_categories)), ceiling(n_categories/ceiling(sqrt(n_categories))) )
  if( !is.null(SampleSize_ctz) ){
    if( !all( dim(SampleSize_ctz)[1:2] == dim(Index_ctl)[1:2] ) ){
      stop("Check input `SampleSize_ctz`")
    }
  }

  if( all(is.numeric(Year_Set)) ){
    year_names = Year_Set
  }else{
    year_names = Year_Set
    Year_Set = 1:length(Year_Set)
  }
  Pretty = function(vec){
    Return = pretty(vec)
    Return = Return[which(Return %in% vec)]
    return(Return)
  }

  # Defaults
  if( is.null(col)) col = rainbow(n_strata)
  if( is.null(col_bounds)) col_bounds = rainbow(n_strata)
  if( is.null(width)) width = mfrow[2] * 3
  if( is.null(height)) height = mfrow[1] * 3

  # Fill in missing
  if( is.null(Year_Set) ) Year_Set = 1:n_years
  if( is.null(Years2Include) ) Years2Include = 1:n_years
  if( is.null(strata_names) ) strata_names = 1:n_strata
  if( is.null(category_names) ) category_names = 1:n_categories
  if( is.null(plot_legend)) plot_legend = ifelse(n_strata>1, TRUE, FALSE)

  # Plot
  Par = combine_lists( default=list(mar=c(2,2,1,0),mgp=c(2,0.5,0),tck=-0.02,yaxs="i",oma=c(2,2,0,0),mfrow=mfrow), input=list(...) )
  if(!is.na(PlotName)){
    png( file=paste0(DirName,PlotName), width=width, height=height, res=200, units="in")  # paste0(DirName,ifelse(DirName=="","","/"),PlotName)
    on.exit( dev.off() )
  }
  par( Par )
  for( z1 in 1:n_categories ){
    # Calculate y-axis limits
    if(scale=="uniform") Ylim = range(Index_ctl[z1,,]%o%c(1,1) + sd_Index_ctl[z1,,]%o%c(-interval_width,interval_width)*1.05, na.rm=TRUE)
    if(scale=="log") Ylim = range(Index_ctl[z1,,]%o%c(1,1) * exp(sd_Index_ctl[z1,,]%o%c(-interval_width,interval_width)*1.05), na.rm=TRUE)
    Ylim = ifelse( is.na(Yrange), Ylim, Yrange )
    Xlim = range(Year_Set) + c(-1,1) * diff(range(Year_Set))/20
    # Plot stuff
    plot_inputs = combine_lists( default=list(1, type="n", xlim=Xlim, ylim=Ylim, xlab="", ylab="", main=ifelse(n_categories>1,category_names[z1],""), xaxt="n"),
      input=plot_args )
    do.call( what=plot, args=plot_inputs )
    for(z3 in 1:n_strata){
      if(scale=="uniform") ybounds = Index_ctl[z1,,z3]%o%c(1,1) + sd_Index_ctl[z1,,z3]%o%c(-interval_width,interval_width)
      if(scale=="log") ybounds = Index_ctl[z1,,z3]%o%c(1,1) * exp(sd_Index_ctl[z1,,z3]%o%c(-interval_width,interval_width))
      if( n_strata==1 ) x_offset = 0
      if( n_strata>=2 ) x_offset = seq(-0.1, 0.1, length=n_strata)[z3]
      plot_lines_defaults = list(y=Index_ctl[z1,,z3], x=Year_Set+x_offset, ybounds=ybounds, type=type, col=col[z3], col_bounds=col_bounds[z3], ylim=Ylim, bounds_type=bounds_type)
      plot_lines_inputs = combine_lists( default=plot_lines_defaults, input=plot_lines_args )
      do.call( what=plot_lines, args=plot_lines_inputs )
    }
    # Plot lines for sample size
    if( !is.null(SampleSize_ctz) ){
      Y2lim = c(1, 1.2) * range(SampleSize_ctz[z1,,], na.rm=TRUE)
      Y2lim = ifelse( is.na(Y2range), Y2lim, Y2range )
      Labels = pretty(Y2lim)
      At = Labels / diff(range(Y2lim,na.rm=TRUE)) * diff(Ylim) + Ylim[1]
      axis( side=4, at=At, labels=Labels )
      for( z3 in 1:dim(SampleSize_ctz)[3] ){
        Y = SampleSize_ctz[z1,,z3] / diff(range(Y2lim,na.rm=TRUE)) * diff(Ylim) + Ylim[1]
        lines( x=Year_Set, y=Y, col=col[z3], lwd=3, lty="dotted" )
        #points( x=Year_Set, y=Y, col=col[z3], cex=1.5 )
      }
    }
    if(plot_legend==TRUE){
      legend( "top", bty="n", fill=rainbow(n_strata), legend=as.character(strata_names), ncol=2 )
    }
    axis( 1, at=Pretty(Year_Set), labels=year_names[match(Pretty(Year_Set),Year_Set)] )
  }
  mtext( side=c(1,2,4), text=c(xlab,ylab,y2lab), outer=TRUE, line=c(0,0) )

  return(invisible(plot_lines_inputs))
}
