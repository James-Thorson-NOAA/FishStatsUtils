
#' Make list of settings
#'
#' \code{make_settings} makes a list of settings for a given purpose
#'
#' This function assembles a default set of user-decisions for a specified modelling purpose. The default settings are guessed based on generic guidance, and should be carefully reviewed for real-world purposes. If the user supplies values for individual settings e.g. \code{FieldConfig}, then these values override the defaults that are provided by interpreting \code{purpose}
#'
#' @param purpose character indicating what purpose is intended for the model, and therefore what default settings are perhaps appropriate. Only currently implemented for \code{purpose="index"}.
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

  # Check defaults
  if( purpose=="index" ){
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=3, "Epsilon1"=3, "Omega2"=3, "Epsilon2"=3)
    }
    if(missing(RhoConfig)){
      RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    }else{
      names(RhoConfig) = c("Beta1","Beta2","Epsilon1","Epsilon2")
    }
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(Options)) Options =  c("SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl" )
  }else{
    stop("'purpose' is currently set up only for index-standardization models")
  }

  # Other defaults
  grid_size_km = 25
  Method = "Mesh"
  if(missing(use_anisotropy)) use_anisotropy = TRUE

  # Bundle and export
  settings = list("Version"=Version, "n_x"=n_x, "Region"=Region, "strata.limits"=strata.limits, "zone"=zone, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
    "OverdispersionConfig"=OverdispersionConfig, "ObsModel"=ObsModel, "vars_to_correct"=vars_to_correct, "Options"=Options, "grid_size_km"=grid_size_km,
    "Method"=Method, "use_anisotropy"=use_anisotropy, "fine_scale"=fine_scale, "bias.correct"=bias.correct )
  return(settings)
}
