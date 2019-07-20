
#' Plot an index
#'
#' \code{plot_index} takes output from a VAST run and plots a panel figure of time-series estimates
#'
#' @param Index_ctl A matrix or array of time-series estimates for multiple categories \code{c}, years \code{t}, and strata \code{l}
#' @param sd_Index_ctl A matrix or array of variances for each estimate
#' @inheritParams plot_biomass_index
#' @inheritParams plot_lines
#' @param Ymin lower bound for y-axis (use \code{Ymin=NA} for using the lower range of estimate intervals)
#' @param ... list of settings to pass to \code{par} when making plot
#'
#' @return NULL
#'
#' @export
plot_index = function( Index_ctl, sd_Index_ctl=array(0,dim(Index_ctl)), Year_Set=NULL, Years2Include=NULL,
  strata_names=NULL, category_names=NULL, scale="uniform",
  plot_legend=NULL, DirName=paste0(getwd(),"/"), PlotName="Index.png",
  interval_width=1, width=NULL, height=NULL, xlab="Year", ylab="Index", bounds_type="whiskers", col=NULL,
  col_bounds=NULL, Ymin=0, ... ){

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
  Par = list( mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i", oma=c(2,2,0,0), mfrow=mfrow, ... )
  png( file=paste0(DirName,PlotName), width=width, height=height, res=200, units="in")  # paste0(DirName,ifelse(DirName=="","","/"),PlotName)
    par( Par )
    for( z1 in 1:n_categories ){
      # Calculate y-axis limits
      if(scale=="uniform") Ylim = range(Index_ctl[z1,,]%o%c(1,1) + sd_Index_ctl[z1,,]%o%c(-interval_width,interval_width)*1.05, na.rm=TRUE)
      if(scale=="log") Ylim = range(Index_ctl[z1,,]%o%c(1,1) * exp(sd_Index_ctl[z1,,]%o%c(-interval_width,interval_width)*1.05), na.rm=TRUE)
      if( !is.na(Ymin) ) Ylim[1] = Ymin
      # Plot stuff
      plot(1, type="n", xlim=range(Year_Set), ylim=Ylim, xlab="", ylab="", main=ifelse(n_categories>1,category_names[z1],""), xaxt="n" )
      for(z3 in 1:n_strata){
        if(scale=="uniform") ybounds = Index_ctl[z1,,z3]%o%c(1,1) + sd_Index_ctl[z1,,z3]%o%c(-interval_width,interval_width)
        if(scale=="log") ybounds = Index_ctl[z1,,z3]%o%c(1,1) * exp(sd_Index_ctl[z1,,z3]%o%c(-interval_width,interval_width))
        FishStatsUtils::plot_lines( y=Index_ctl[z1,,z3], x=Year_Set+seq(-0.1,0.1,length=n_strata)[z3], ybounds=ybounds, type="b", col=col[z3], col_bounds=col_bounds[z3], ylim=Ylim, bounds_type=bounds_type)
      }
      if(plot_legend==TRUE) legend( "top", bty="n", fill=rainbow(n_strata), legend=as.character(strata_names), ncol=2 )
      axis( 1, at=Pretty(Year_Set), labels=year_names[match(Pretty(Year_Set),Year_Set)] )
    }
    mtext( side=1:2, text=c(xlab,ylab), outer=TRUE, line=c(0,0) )
  dev.off()
}
