
#' @title
#' Make list of settings
#'
#' @description
#' \code{make_settings} makes a list of settings for a given purpose
#'
#' @inheritParams VAST::make_data
#' @inheritParams make_extrapolation_info
#' @inheritParams make_spatial_info
#'

#' @export
make_settings = function( n_x, Region, purpose="index", fine_scale=TRUE,
  strata.limits=data.frame('STRATA'="All_areas"), zone=NA ){

  # Get version
  Version = FishStatsUtils::get_latest_version()

  # Check defaults
  if( purpose=="index" ){
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      FieldConfig = c("Omega1"=3, "Epsilon1"=3, "Omega2"=3, "Epsilon2"=3)
    }
    RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    ObsModel = c(1,1)
    vars_to_correct = c( "Index_cyl" )
    bias.correct = TRUE
    Options =  c("SD_site_logdensity"=FALSE, "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE )
  }else{
    stop("'purpose' is currently set up only for index-standardization models")
  }

  # Other defaults
  grid_size_km = 25
  Method = "Mesh"
  use_anisotropy = TRUE

  # Bundle and export
  settings = list("Version"=Version, "n_x"=n_x, "Region"=Region, "strata.limits"=strata.limits, "zone"=zone, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
    "OverdispersionConfig"=OverdispersionConfig, "ObsModel"=ObsModel, "vars_to_correct"=vars_to_correct, "Options"=Options, "grid_size_km"=grid_size_km,
    "Method"=Method, "use_anisotropy"=use_anisotropy, "fine_scale"=fine_scale, "bias.correct"=bias.correct )
  return(settings)
}
