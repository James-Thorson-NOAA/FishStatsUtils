
#' @title
#' Plot location of sampling data
#'
#' @description
#' \code{plot_data} produces diagnostics plots for the spatial distribution of data and knots
#'
#' @param Extrapolation_List Output from \code{Prepare_Extrapolation_Data_Fn}
#' @param Spatial_List Output from \code{Spatial_Information_Fn}
#' @param Data_Geostat data-frame of data (with columns 'E_km', 'N_km', 'Year', 'Lon', 'Lat' at a minimum)
#' @param PlotDir Directory for plots
#' @param maxpanel The maximum number of rows or columns you want per panel. The default
#' generates 2 x 2 panels of years, where each set is saved as a separate file.
#' @param ... addition inputs to \code{plot}
#'

#' @export
plot_data = function( Extrapolation_List, Spatial_List, Data_Geostat, PlotDir=paste0(getwd(),"/"),
  Plot1_name="Data_and_knots.png", Plot2_name="Data_by_year.png", col=rep("red",nrow(Data_Geostat)), cex=0.01, maxpanel = 2, ...){

  # Override defaults
  if( length(cex) == 1 ){
    cex = rep( cex, nrow(Data_Geostat) )
  }else{
    if(length(cex)!=nrow(Data_Geostat)) stop("input `cex` has wrong length")
  }
  if( length(pch) == 1 ){
    pch = rep( pch, nrow(Data_Geostat) )
  }else{
    if(length(pch)!=nrow(Data_Geostat)) stop("input `pch` has wrong length")
  }
  # Plot data and grid
  png( file=paste0(PlotDir,Plot1_name), width=6, height=6, res=200, units="in")
    par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(1.75,0.25,0) )
    plot( Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('Lon','Lat')], cex=0.01, main="Extrapolation (Lat-Lon)" )
    map( "world", add=TRUE )
    if( !any(is.na(Extrapolation_List$Data_Extrap[,c('E_km','N_km')])) ){
      plot( Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('E_km','N_km')], cex=0.01, main="Extrapolation (North-East)", ylab = "Northing (km)", xlab = "Easting (km)" )
    }
    plot( Spatial_List$loc_x, col="red", pch=20, main="Knots (North-East)")
    if( all(c('E_km','N_km')%in%names(Data_Geostat)) ){
      plot( Data_Geostat[,c('E_km','N_km')], col="blue", pch=20, cex=0.1, main="Data (North-East)", ylab = "Northing (km)", xlab = "Easting (km)")
    }
  dev.off()

  # Plot data by year
  # Use Data_Geostat, instead of TmbData, because I want raw locations, not locations at knots
  if(missing(Year_Set)) Year_Set = unique(Data_Geostat[,'Year'])[order(unique(Data_Geostat[,'Year']))]
    Nrow = ceiling( sqrt(length(Year_Set)) )
    Ncol = ceiling( length(Year_Set)/Nrow )
  npanels <- ifelse(Nrow > maxpanel | Ncol > maxpanel,
    ceiling(length(Year_Set) / (maxpanel * maxpanel)), 1)
  yearsall <- Year_Set
  for (panel in 1:npanels) {
    if (npanels != 1) Plot2_name <- gsub("(.)[0-9]*\\.png",
      paste0("\\1", sprintf(paste0("%0", nchar(npanels), "d"), panel), ".png"), Plot2_name)
    Year_Set <- na.omit(yearsall[1:(maxpanel * maxpanel)])
    Nrow <- ceiling( sqrt(length(Year_Set)) )
    Ncol = ceiling( length(Year_Set)/Nrow )
  png( file=paste0(PlotDir,Plot2_name), width=Ncol*2, height=Nrow*2, res=200, units="in")
    par( mfrow=c(Nrow,Ncol), mar=c(0,0,2,0), mgp=c(1.75,0.45,0), oma=c(3,3,0.1,0.1) )
    for( t in 1:length(Year_Set) ){
      Which = which( Data_Geostat[,'Year'] == Year_Set[t] )
      plot( x=Data_Geostat[Which,'Lon'], y=Data_Geostat[Which,'Lat'], cex=cex[Which], main=Year_Set[t], xlim=range(Data_Geostat[,'Lon']), ylim=range(Data_Geostat[,'Lat']), xaxt="n", yaxt="n", col=col[Which], pch=pch[Which], ... )
      map( "world", add=TRUE, fill=TRUE, col="grey" )
      if( t>(length(Year_Set)-Ncol) ) axis(1, at = axisTicks(range(Data_Geostat[Which,'Lon']), log = FALSE))
      if( t%%Ncol == 1 ) axis(2)
      mtext( side=c(1,2), text=c("Longitude","Latitude"), outer=TRUE, line=1.35)
    }
  dev.off()
  yearsall <- yearsall[-(1:(maxpanel * maxpanel))]
  }
}
