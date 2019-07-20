
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
#' @inheritParams TMBhelper::fit_tmb
#' @param extrapolation_args tagged list of optional arguments to pass to \code{FishStatsUtils::make_extrapolation_info}
#' @param spatial_args tagged list of optional arguments to pass to \code{FishStatsUtils::make_spatial_info}
#' @param optimize_args tagged list of optional arguments to pass to \code{TMBhelper::Optimize}
#' @param model_args tagged list of optional arguments to pass to \code{VAST::make_model}
#' @param run_model Boolean indicating whether to run the model or simply return the inputs and built TMB object
#' @param test_fit Boolean indicating whether to apply \code{VAST::check_fit} before calculating standard errors, to test for parameters hitting bounds etc; defaults to TRUE
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
  extrapolation_args=list(), spatial_args=list(), optimize_args=list(), model_args=list(),
  silent=TRUE, run_model=TRUE, test_fit=TRUE, ... ){

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
  data_frame = data.frame( "Lat_i"=Lat_i, "Lon_i"=Lon_i, "a_i"=a_i, "v_i"=v_i, "b_i"=b_i, "t_i"=t_iz, "c_iz"=c_iz )
  # Decide which years to plot
  year_labels = seq( min(t_iz), max(t_iz) )
  years_to_plot = which( year_labels %in% t_iz )

  # Save record
  dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
  #save( settings, file=file.path(working_dir,"Record.RData"))
  capture.output( settings, file=file.path(working_dir,"settings.txt"))

  # Build extrapolation grid
  message("\n### Making extrapolation-grid")
  extrapolation_args = combine_lists( input=extrapolation_args, default=list(Region=settings$Region, strata.limits=settings$strata.limits, zone=settings$zone) )
  extrapolation_list = do.call( what=make_extrapolation_info, args=extrapolation_args )

  # Build information regarding spatial location and correlation
  message("\n### Making spatial information")
  spatial_args = combine_lists( input=spatial_args, default=list(grid_size_km=settings$grid_size_km, n_x=settings$n_x, Method=settings$Method, Lon_i=Lon_i, Lat_i=Lat_i,
    Extrapolation_List=extrapolation_list, DirPath=working_dir, Save_Results=TRUE, fine_scale=settings$fine_scale) )
  #spatial_list = make_spatial_info( grid_size_km=settings$grid_size_km, n_x=settings$n_x, Method=settings$Method, Lon_i=Lon_i, Lat_i=Lat_i,
  #  Extrapolation_List=extrapolation_list, DirPath=working_dir, Save_Results=TRUE, fine_scale=settings$fine_scale )
  spatial_list = do.call( what=make_spatial_info, args=spatial_args )

  # Build data
  message("\n### Making data object") # VAST::
  data_list = VAST::make_data("Version"=settings$Version, "FieldConfig"=settings$FieldConfig, "OverdispersionConfig"=settings$OverdispersionConfig,
    "RhoConfig"=settings$RhoConfig, "VamConfig"=settings$VamConfig, "ObsModel"=settings$ObsModel, "c_iz"=c_iz, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i,
    "s_i"=spatial_list$knot_i-1, "t_iz"=t_iz, "spatial_list"=spatial_list, "Options"=settings$Options, "Aniso"=settings$use_anisotropy,
    Xconfig_zcp=Xconfig_zcp, X_gtp=X_gtp, X_itp=X_itp, Q_ik=Q_ik, ... )

  # Build object
  message("\n### Making TMB object")
  model_args = combine_lists( input=model_args, default=list("TmbData"=data_list, "RunDir"=working_dir, "Version"=settings$Version,
    "RhoConfig"=settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method) )
  tmb_list = do.call( what=VAST::make_model, args=model_args )  # VAST::
  if(silent==TRUE) tmb_list$Obj$env$beSilent()

  # Run the model or optionally don't
  if( run_model==FALSE ){
    # Build and output
    Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
      "data_list"=data_list, "tmb_list"=tmb_list, "year_labels"=year_labels, "years_to_plot"=years_to_plot,
      "settings"=settings, "extrapolation_args"=extrapolation_args, "model_args"=model_args)
    class(Return) = "fit_model"
    return(Return)
  }

  # Optimize object
  message("\n### Estimating parameters")
  optimize_args_phase1 = combine_lists( default=optimize_args, input=list(obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
    savedir=NULL, getsd=FALSE, newtonsteps=0, bias.correct=FALSE, quiet=TRUE,
    control=list(eval.max=10000,iter.max=10000,trace=1), loopnum=2) )
  parameter_estimates = do.call( what=TMBhelper::fit_tmb, args=optimize_args_phase1 )

  # Check fit of model (i.e., evidence of non-convergence based on bounds, approaching zero, etc)
  if(exists("check_fit") & test_fit==TRUE ){
    problem_found = VAST::check_fit( parameter_estimates )
    if( problem_found==TRUE ){
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try\n", call.=FALSE)
    }
  }

  # Restart estimates after checking parameters
  optimize_args_phase2 = combine_lists( input=optimize_args, default=list(obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
    savedir=working_dir, bias.correct=settings$bias.correct, newtonsteps=newtonsteps,
    bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=settings$vars_to_correct),
    control=list(eval.max=10000,iter.max=10000,trace=1), loopnum=1) )
  optimize_args_phase2 = combine_lists( input=list(startpar=parameter_estimates$par), default=optimize_args_phase2 )
  parameter_estimates = do.call( what=TMBhelper::fit_tmb, args=optimize_args_phase2 )

  # Extract standard outputs
  Report = tmb_list$Obj$report()
  ParHat = tmb_list$Obj$env$parList( parameter_estimates$par )

  # Build and output
  Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
    "data_list"=data_list, "tmb_list"=tmb_list, "parameter_estimates"=parameter_estimates, "Report"=Report,
    "ParHat"=ParHat, "year_labels"=year_labels, "years_to_plot"=years_to_plot, "settings"=settings,
    "extrapolation_args"=extrapolation_args, "model_args"=model_args, "optimize_args"=optimize_args)
  class(Return) = "fit_model"
  return( Return )
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{fit_model}}
#' @param ... Not used
#' @return NULL
#' @method print fit_model
#' @export
print.fit_model <- function(x, ...)
{
  cat("fit_model(.) result\n")
  if( "parameter_estimates" %in% names(x) ){
    print( x$parameter_estimates )
  }else{
    cat("`parameter_estimates` not available in `fit_model`\n")
  }
  invisible(x$parameter_estimates)
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param fit Output from \code{\link{fit_model}}
#' @param what String specifying what elements of results to plot;  options include `extrapolation_grid`, `spatial_mesh`, and `results`
#' @param ... Arguments passed to \code{\link{plot_results}}
#' @return NULL
#' @method plot fit_model
#' @export plot
#' @export
plot.fit_model <- function(x, what="results", ...)
{
  ## Plot extrapolation-grid
  if( length(grep(what, "extrapolation_grid")) ){
    cat("\n### Running `plot.make_extrapolation_info`\n")
    plot( x$extrapolation_list )
    return(invisible(NULL))
  }

  ## Plot extrapolation-grid
  if( length(grep(what, c("spatial_info","inla_mesh"))) ){
    cat("\n### Running `plot.make_spatial_info`\n")
    plot( x$spatial_list )
    return(invisible(NULL))
  }

  # diagnostic plots
  if( length(grep(what, "results")) ){
    cat("\n### Running `plot_results`\n")
    ans = plot_results( x, ... )
    return(invisible(ans))
  }

  stop( "input `what` not matching available options" )
}

#' Extract summary of spatial estimates
#'
#' @title Extract spatial estimates
#' @param fit Output from \code{\link{fit_model}}
#' @param what Boolean indicating what to summarize; only option is `density`
#' @param ... Not used
#' @return NULL
#' @method summary fit_model
#' @export
summary.fit_model <- function(x, what="density", ...)
{
  ans = NULL

  if( tolower(what) == "density" ){
    # Load location of extrapolation-grid
    ans[["extrapolation_grid"]] = print( x$extrapolation_list, quiet=TRUE )

    # Load density estimates
    if( "D_gcy" %in% names(x$Report)){
      ans[["Density_array"]] =  x$Report$D_gcy
      if( !( x$settings$fine_scale==TRUE | x$spatial_list$Method=="Stream_network" ) ){
        index_tmp = x$spatial_list$NN_Extrap$nn.idx[ which(x$extrapolation_list[["Area_km2_x"]]>0), 1 ]
        ans[["Density_array"]] = ans[["Density_array"]][ index_tmp,,,drop=FALSE]
      }
      dimnames(ans[["Density_array"]]) = list( rownames(ans[["extrapolation_grid"]]), paste0("Category_",1:dim(ans[["Density_array"]])[[2]]), x$year_labels )
      # Expand as grid
      Density_dataframe = expand.grid("Grid"=1:dim(ans[["Density_array"]])[[1]], "Category"=dimnames(ans[["Density_array"]])[[2]], "Year"=dimnames(ans[["Density_array"]])[[3]])
      Density_dataframe = cbind( Density_dataframe, ans[["extrapolation_grid"]][Density_dataframe[,'Grid'],], "Density"=as.vector(ans[["Density_array"]]) )
      ans[["Density_dataframe"]] = Density_dataframe
      rownames(Density_dataframe) = NULL
      cat("\n### Printing head of and tail `Density_dataframe`, and returning data frame in output object")
      print(head(Density_dataframe))
      print(tail(Density_dataframe))
    }else{
      stop( "`summary.fit_model` not implemented for the version of `VAST` being used" )
    }
  }

  if( tolower(what) %in% c("parhat","estimates") ){
    ans[["estimates"]] = x$ParHat
    cat("\n### Printing slots of `ParHat`, and returning list in output object")
    print(names(x$ParHat))
  }

  if( is.null(ans) ){
    stop( "`summary.fit_model` not implemented for inputted value of argument `what`" )
  }

  # diagnostic plots
  return(invisible(ans))
}



