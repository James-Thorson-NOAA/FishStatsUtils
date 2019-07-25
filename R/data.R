
#' Data to demonstrate and test stream-network options
#'
#' Data sufficient to demonstrate stream-network options in VAST, using \code{Method="Stream_network"},
#'  and also helpful to do integrated testing of this feature
#'
#' \itemize{
#'   \item observations data-frame of biological sampling data
#'   \item network data-frame defining stream-network connectivity, used to define Ornstein-Uhlenbeck covariance function
#'   \item habitat data-frame containing density covariates
#' }
#'
#' @name stream_network_eel_example
#' @docType data
#' @author Merrill Rudd
#' @usage data(stream_network_eel_example)
#' @keywords data
NULL

#' Data to demonstrate and test covariate effects
#'
#' Data sufficient to demonstrate spatially varying coefficient models
#'
#' \itemize{
#'   \item sampling_data data-frame of biological sampling data and associated covariate measurements
#'   \item Region region for model demo
#'   \item strata.limits user-specified stratification of results
#'   \item X_xtp 3D array of covariates, specified at knots given the use of \code{fine_scale=FALSE}
#' }
#'
#' @name GOA_pcod_covariate_example
#' @docType data
#' @author Dave McGowan
#' @usage data(GOA_pcod_covariate_example)
#' @keywords data
NULL
