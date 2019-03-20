
#' @title
#' Plot results
#'
#' @description
#' \code{plot_results} plots diagnostics, results, and indices for a given fitted model
#'
#' @inheritParams plot_maps
#'

#' @export
plot_results = function( fit, settings, working_dir=paste0(getwd(),"/"), year_labels=fit$year_labels, years_to_plot=fit$years_to_plot ){

  # plot data
  #plot_data(Extrapolation_List=fit$extrapolation_list, Spatial_List=fit$spatial_list, Data_Geostat=Data_Geostat, PlotDir=working_dir )

  # PLot settings
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  # Plot diagnostic for encounter probability
  message("\n### Making plot of encounter probability")
  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=working_dir)

  # Plot anisotropy
  message("\n### Making plot of anisotropy")
  plot_anisotropy( FileName=paste0(working_dir,"Aniso.png"), Report=fit$Report, TmbData=fit$data_list )

  # Plot index
  message("\n### Making plot of abundance index")
  Index = plot_biomass_index( DirName=working_dir, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=year_labels, Years2Include=years_to_plot, use_biascorr=TRUE )

  # Plot range indices
  message("\n### Making plot of spatial indices")
  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=working_dir, Year_Set=year_labels)

  # Plot densities
  message("\n### Making plot of densities")
  Dens_xt = plot_maps(plot_set=c(3), MappingDetails=map_list[["MappingDetails"]], Report=fit$Report, Sdreport=Opt$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=working_dir, Year_Set=year_labels, Years2Include=years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)

  # Plot quantile-quantile plot
  message("\n### Making Q-Q plot")
  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=working_dir )

  # Pearson residuals
  if( "n_x" %in% names(fit$data_list) ){
    message("\n### Making plot of Pearson residuals")
    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, Report=fit$Report, Q=Q, savedir=working_dir, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=working_dir, Year_Set=year_labels, Years2Include=years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)
  }else{
    message("\n### Skipping plot of Pearson residuals")
  }

  # return
  Return = list( "Q"=Q, "Index"=Index, "Dens_xt"=Dens_xt )
  return( invisible(Return) )
}

