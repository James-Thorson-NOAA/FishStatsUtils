
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  if( getOption("repos")["CRAN"] == "@CRAN@" ){
    options(repos = c("CRAN" = "http://cran.us.r-project.org"))
  }
  if( !"INLA" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: INLA...")
    #utils::install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
    utils::install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  }
  #if( !"TMB" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing TMB...")
  #  devtools::install_github("kaskr/adcomp/TMB")
  #}
  #if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: TMBhelper...")
  #  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  #}
}

#' Copy of FishStatsUtils::plot_loadings
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::plot_loadings} to see list of arguments and usage
#' @export
PlotLoadings = function( ... ){
  .Deprecated( new="FishStatsUtils::plot_loadings" )
  FishStatsUtils::plot_loadings( ... )
}

