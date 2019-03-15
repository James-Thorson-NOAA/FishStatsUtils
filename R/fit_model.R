
#' @title
#' Fit VAST to data
#'
#' @description
#' \code{fit_model} fits a spatio-temporal model to data
#'
#' @inheritParams VAST::make_data
#'

#' @export
fit_model = function( settings, Lat_i, Lon_i, t_i, c_i, b_i, a_i, v_i, working_dir=paste0(getwd(),"/"),
  getsd=TRUE, newtonsteps=1, X_xtp=NULL, Xconfig_zcp=NULL, X_gtp=NULL, X_itp=NULL, Q_ik=NULL, ... ){

  # Assemble inputs
  data_frame = data.frame( "Lat_i"=Lat_i, "Lon_i"=Lon_i, "a_i"=a_i, "v_i"=v_i, "b_i"=b_i )

  # Save record
  save( settings, file=file.path(working_dir,"Record.RData"))
  capture.output( settings, file=file.path(working_dir,"Record.txt"))

  # Build extrapolation grid
  message("\n### Making extrapolation-grid")
  extrapolation_list = make_extrapolation_info( Region=Region, strata.limits=settings$strata.limits, ... )

  # Build information regarding spatial location and correlation
  message("\n### Making spatial information")
  spatial_list = make_spatial_info( grid_size_km=settings$grid_size_km, n_x=settings$n_x, Method=settings$Method, Lon_i=Lon_i, Lat_i=Lat_i,
    Extrapolation_List=extrapolation_list, DirPath=working_dir, Save_Results=TRUE, fine_scale=settings$fine_scale )

  # Build data
  message("\n### Making data object")
  data_list = make_data("Version"=settings$Version, "FieldConfig"=settings$FieldConfig, "OverdispersionConfig"=settings$OverdispersionConfig,
    "RhoConfig"=settings$RhoConfig, "ObsModel"=settings$ObsModel, "c_i"=c_i, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i,
    "s_i"=spatial_list$knot_i-1, "t_i"=t_i, "spatial_list"=spatial_list, "Options"=settings$Options, "Aniso"=settings$use_anisotropy )
  #return( list("data_list"=data_list, "spatial_list"=spatial_list) )

  # Build object
  message("\n### Making TMB object")
  tmb_list = make_model("TmbData"=data_list, "RunDir"=working_dir, "Version"=settings$Version, "RhoConfig"=settings$RhoConfig,
    "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method)

  # Optimize object
  message("\n### Estimating parameters")
  parameter_estimates = TMBhelper::Optimize( obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
    getsd=getsd, savedir=working_dir, bias.correct=settings$bias.correct, newtonsteps=newtonsteps,
    bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=settings$vars_to_correct) )

  # Extract standard outputs
  Report = tmb_list$Obj$report()
  ParHat = tmb_list$Obj$env$parList( parameter_estimates$par )

  # Build and output
  Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list, "data_list"=data_list,
    "tmb_list"=tmb_list, "parameter_estimates"=parameter_estimates, "Report"=Report, "ParHat"=ParHat)
  return( Return )
}
