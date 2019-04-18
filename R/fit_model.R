
#' Fit VAST to data
#'
#' \code{fit_model} fits a spatio-temporal model to data
#'
#' This function is the user-interface for the functions that determine the extrapolation-grid, define spatial objects, assemble data, build model, and estimate parameters.
#'
#' @param settings Output from \code{make_settings}
#' @inheritParams make_extrapolation_info
#' @inheritParams make_spatial_info
#' @inheritParams VAST::make_data
#' @inheritParams VAST::make_model
#' @inheritParams TMBhelper::Optimize
#' @param extrapolation_args tagged list of optional arguments to pass to \code{FishStatsUtils::make_extrapolation_info}
#' @param optimize_args tagged list of optional arguments to pass to \code{TMBhelper::Optimize}
#' @param model_args tagged list of optional arguments to pass to \code{VAST::make_model}
#' @param run_model Boolean indicating whether to run the model or simply return the inputs and built TMB object
#' @param ... additional parameters to pass to \code{VAST::make_data}
#'
#' @return Returns a tagged list of internal objects, the TMB object, and slot \code{parameter_estimates} containing the MLE estimates
#'
#' @family wrapper functions
#' @seealso \code{?VAST} for general documentation, \code{?make_settings} for generic settings, \code{?fit_model} for model fitting, and \code{?plot_results} for generic plots
#'
#' @examples
#' \dontrun{
#' # Load packages
#' library(TMB)
#' library(VAST)
#'
#' # load data set
#' # see `?load_example` for list of stocks with example data
#' # that are installed automatically with `FishStatsUtils`.
#' example = load_example( data_set="EBS_pollock" )
#'
#' # Make settings
#' settings = make_settings( n_x=50, Region=example$Region, purpose="index",
#'   strata.limits=example$strata.limits )
#'
#' # Run model
#' fit = fit_model( "settings"=settings, "Lat_i"=example$sampling_data[,'Lat'],
#'   "Lon_i"=example$sampling_data[,'Lon'], "t_i"=example$sampling_data[,'Year'],
#'   "c_i"=rep(0,nrow(example$sampling_data)), "b_i"=example$sampling_data[,'Catch_KG'],
#'   "a_i"=example$sampling_data[,'AreaSwept_km2'], "v_i"=example$sampling_data[,'Vessel'] )
#'
#' # Plot results
#' plot_results( settings=settings, fit=fit )
#' }
#'
#' @export
fit_model = function( settings, Lat_i, Lon_i, t_iz, c_iz, b_i, a_i,
  v_i=rep(0,length(b_i)), working_dir=paste0(getwd(),"/"),
  Xconfig_zcp=NULL, X_gtp=NULL, X_itp=NULL, Q_ik=NULL, newtonsteps=1,
  extrapolation_args=list(), optimize_args=list(), model_args=list(), silent=TRUE, run_model=TRUE, ... ){

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

  # Assemble inputs
  data_frame = data.frame( "Lat_i"=Lat_i, "Lon_i"=Lon_i, "a_i"=a_i, "v_i"=v_i, "b_i"=b_i )
  # Decide which years to plot
  year_labels = seq( min(t_iz), max(t_iz) )
  years_to_plot = which( unique(t_iz) %in% sort(unique(t_iz)))

  # Save record
  dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
  save( settings, file=file.path(working_dir,"Record.RData"))
  capture.output( settings, file=file.path(working_dir,"Record.txt"))

  # Build extrapolation grid
  message("\n### Making extrapolation-grid")
  extrapolation_args = combine_lists( input=extrapolation_args, default=list(Region=settings$Region, strata.limits=settings$strata.limits, zone=settings$zone) )
  extrapolation_list = do.call( what=make_extrapolation_info, args=extrapolation_args )

  # Build information regarding spatial location and correlation
  message("\n### Making spatial information")
  spatial_list = make_spatial_info( grid_size_km=settings$grid_size_km, n_x=settings$n_x, Method=settings$Method, Lon_i=Lon_i, Lat_i=Lat_i,
    Extrapolation_List=extrapolation_list, DirPath=working_dir, Save_Results=TRUE, fine_scale=settings$fine_scale )

  # Build data
  message("\n### Making data object") # VAST::
  data_list = VAST::make_data("Version"=settings$Version, "FieldConfig"=settings$FieldConfig, "OverdispersionConfig"=settings$OverdispersionConfig,
    "RhoConfig"=settings$RhoConfig, "ObsModel"=settings$ObsModel, "c_iz"=c_iz, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i,
    "s_i"=spatial_list$knot_i-1, "t_iz"=t_iz, "spatial_list"=spatial_list, "Options"=settings$Options, "Aniso"=settings$use_anisotropy,
    Xconfig_zcp=Xconfig_zcp, X_gtp=X_gtp, X_itp=X_itp, Q_ik=Q_ik, ... )

  # Build object
  message("\n### Making TMB object")
  model_args = combine_lists( input=model_args, default=list("TmbData"=data_list, "RunDir"=working_dir, "Version"=settings$Version,
    "RhoConfig"=settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method) )
  tmb_list = do.call( what=VAST::make_model, args=model_args )  # VAST::

  # Run the model or optionally don't
  if( run_model==TRUE ){
    # Optimize object
    message("\n### Estimating parameters")
    optimize_args = combine_lists( input=optimize_args, default=list(obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
      savedir=working_dir, bias.correct=settings$bias.correct, newtonsteps=newtonsteps,
      bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=settings$vars_to_correct),
      control=list(eval.max=10000,iter.max=10000,trace=1)) )
    if(silent==TRUE) tmb_list$Obj$env$beSilent()
    parameter_estimates = do.call( what=TMBhelper::Optimize, args=optimize_args )

    # Extract standard outputs
    Report = tmb_list$Obj$report()
    ParHat = tmb_list$Obj$env$parList( parameter_estimates$par )

    # Build and output
    Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
      "data_list"=data_list, "tmb_list"=tmb_list, "parameter_estimates"=parameter_estimates, "Report"=Report,
      "ParHat"=ParHat, "year_labels"=year_labels, "years_to_plot"=years_to_plot, "settings"=settings,
      "extrapolation_args"=extrapolation_args, "model_args"=model_args, "optimize_args"=optimize_args)
  }else{
    # Build and output
    Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
      "data_list"=data_list, "tmb_list"=tmb_list, "year_labels"=year_labels, "years_to_plot"=years_to_plot,
      "settings"=settings, "extrapolation_args"=extrapolation_args, "model_args"=model_args)
  }
  return( Return )
}
