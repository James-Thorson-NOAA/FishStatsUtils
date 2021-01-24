
#' @title
#' Plot location of sampling data
#'
#' @description
#' \code{plot_data} produces diagnostics plots for the spatial distribution of data and knots
#'
#' @inheritParams plot_variable
#' @inheritParams plot_maps
#'
#' @param Extrapolation_List Output from \code{Prepare_Extrapolation_Data_Fn}
#' @param Spatial_List Output from \code{Spatial_Information_Fn}
#' @param Data_Geostat data-frame of data (with columns 'E_km', 'N_km', 'Year', 'Lon', 'Lat' at a minimum)
#' @param PlotDir Directory for plots
#' @param ... addition inputs to \code{plot}
#'

#' @export
plot_data <-
function( Extrapolation_List,
          Spatial_List,
          Data_Geostat = NULL,
          Lat_i = Data_Geostat[,'Lat'],
          Lon_i = Data_Geostat[,'Lon'],
          Year_i = Data_Geostat[,'Year'],
          PlotDir = paste0(getwd(),"/"),
          Plot1_name = "Data_and_knots.png",
          Plot2_name = "Data_by_year.png",
          col = "red",
          cex = 0.01,
          pch = 19,
          year_labels,
          projargs = '+proj=longlat',
          map_resolution = "medium",
          land_color = "grey",
          country = NULL,
          ...){

  # Check for issues
  if( is.null(Lat_i) | is.null(Lon_i) | is.null(Year_i) ){
    stop("Problem with inputs")
  }

  # Override defaults
  if( length(cex) == 1 ){
    cex = rep( cex, length(Year_i) )
  }else{
    if(length(cex)!=length(Year_i)) stop("input `cex` has wrong length")
  }
  if( length(pch) == 1 ){
    pch = rep( pch, length(Year_i) )
  }else{
    if(length(pch)!=length(Year_i)) stop("input `pch` has wrong length")
  }
  if( length(col) == 1 ){
    col = rep( col, length(Year_i) )
  }else{
    if(length(col)!=length(Year_i)) stop("input `col` has wrong length")
  }

  # CRS for original and new projections
  CRS_proj = sp::CRS( projargs )

  # Data for mapping
  map_data = rnaturalearth::ne_countries(scale=switch(map_resolution, "low"=110, "medium"=50, "high"=10, 50), country=country)
  map_data = sp::spTransform(map_data, CRSobj=CRS_proj)

  # Plot data and grid
  if( !missing(Extrapolation_List) & !missing(Spatial_List) ){
    png( file=paste0(PlotDir,Plot1_name), width=6, height=6, res=200, units="in")
      par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(1.75,0.25,0) )
      which_rows = which( Extrapolation_List[["Area_km2_x"]]>0 & rowSums(Extrapolation_List[["a_el"]])>0 )
      plot( Extrapolation_List$Data_Extrap[which_rows,c('Lon','Lat')], cex=0.01, main="Extrapolation (Lat-Lon)" )
      sp::plot( map_data, col=land_color, add=TRUE )
      if( !any(is.na(Extrapolation_List$Data_Extrap[,c('E_km','N_km')])) ){
        plot( Extrapolation_List$Data_Extrap[which_rows,c('E_km','N_km')], cex=0.01, main="Extrapolation (North-East)" )
      }
      plot( Spatial_List$loc_x, col="red", pch=20, main="Knots (North-East)")
    dev.off()
  }

  # Plot data by year
  # Use Data_Geostat, instead of TmbData, because I want raw locations, not locations at knots
  if(missing(year_labels)) year_labels = min(Year_i):max(Year_i)
  if( !any(unique(Year_i) %in% year_labels) ) year_labels = sort(unique(Year_i))
    Nrow = ceiling( sqrt(length(year_labels)) )
    Ncol = ceiling( length(year_labels)/Nrow )
  if(!is.null(Plot2_name)) png( file=paste0(PlotDir,Plot2_name), width=Ncol*2, height=Nrow*2, res=200, units="in")
    par( mfrow=c(Nrow,Ncol), mar=c(0,0,2,0), mgp=c(1.75,0.25,0), oma=c(4,4,0,0) )
    for( t in 1:length(year_labels) ){
      plot( 1, type="n", xlim=range(Lon_i), ylim=range(Lat_i), main=year_labels[t], xaxt="n", yaxt="n" )
      sp::plot( map_data, col=land_color, add=TRUE )
      Which = which( Year_i == year_labels[t] )
      if( length(Which)>0 ){
        points( x=Lon_i[Which], y=Lat_i[Which], cex=cex[Which], col=col[Which], pch=pch[Which], ... )
      }
      if( t>(length(year_labels)-Ncol) ) axis(1)
      if( t%%Ncol == 1 ) axis(2)
      mtext( side=c(1,2), text=c("Longitude","Latitude"), outer=TRUE, line=1)
    }
  if(!is.null(Plot2_name))dev.off()
}
