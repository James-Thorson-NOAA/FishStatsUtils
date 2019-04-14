
#' Make list of settings
#'
#' \code{make_settings} makes a list of settings for a given purpose
#'
#' This function assembles a default set of user-decisions for a specified modelling purpose. The default settings are guessed based on generic guidance, and should be carefully reviewed for real-world purposes. If the user supplies values for individual settings e.g. \code{FieldConfig}, then these values override the defaults that are provided by interpreting \code{purpose}
#'
#' @param purpose character indicating what purpose is intended for the model, and therefore what default settings are perhaps appropriate. Only currently implemented for \code{purpose="index"} or \code{purpose="condition_and_density"}.
#' @inheritParams VAST::make_data
#' @inheritParams make_extrapolation_info
#' @inheritParams make_spatial_info
#' @inheritParams Convert_LL_to_UTM_Fn
#' @param use_anisotropy Boolean indicating whether to estimate two additional parameters representing geometric anisotropy
#' @param vars_to_correct a character-vector listing which parameters to include for bias-correction, as passed to \code{TMBhelper::Optimize}
#'
#' @return Tagged list containing default settings for a given purpose, use \code{names} on output to see list of settings.
#'
#' @family wrapper functions
#' @seealso \code{?VAST} for general documentation, \code{?make_settings} for generic settings, \code{?fit_model} for model fitting, and \code{?plot_results} for generic plots
#'
#' @export
make_settings = function( n_x, Region, purpose="index", fine_scale=TRUE,
  strata.limits=data.frame('STRATA'="All_areas"), zone=NA, FieldConfig, RhoConfig,
  OverdispersionConfig, ObsModel, bias.correct, Options, use_anisotropy,
  vars_to_correct, Version ){

  # Get version
  if(missing(Version)) Version = FishStatsUtils::get_latest_version()

  # Index standardization
  if( purpose=="index" ){
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(Options)) Options =  c("SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl" )
  }

  # Condition and density
  if( purpose=="condition_and_density" ){
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( c(2,2,"IID",0,0,"IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=0, "Epsilon2"=0)
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,4)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(Options)) Options =  c("SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl" )
  }

  # Check for bad input
  if( !( purpose %in% c("index","condition_and_density")) ){
    stop("'purpose' is currently set up only for index-standardization models and correlations between condition and density")
  }

  # Other defaults
  grid_size_km = 25
  Method = "Mesh"
  if(missing(use_anisotropy)) use_anisotropy = TRUE

  # Default naming
  names(RhoConfig) = c("Beta1","Beta2","Epsilon1","Epsilon2")

  # Bundle and export
  settings = list("Version"=Version, "n_x"=n_x, "Region"=Region, "strata.limits"=strata.limits, "zone"=zone, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
    "OverdispersionConfig"=OverdispersionConfig, "ObsModel"=ObsModel, "vars_to_correct"=vars_to_correct, "Options"=Options, "grid_size_km"=grid_size_km,
    "Method"=Method, "use_anisotropy"=use_anisotropy, "fine_scale"=fine_scale, "bias.correct"=bias.correct )
  return(settings)
}
