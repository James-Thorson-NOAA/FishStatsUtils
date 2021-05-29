

#' Format density covariate matrix
#'
#' \code{amend_output} uses a formula interface to generate covariates
#'
#' @export
amend_output <-
function( TmbData,
          Report,
          Sdreport = NULL,
          year_labels = NULL,
          category_names = NULL ){

  # Local functions
  add_dimnames = function( Report, report_names, dimnames ){
    for( i in seq_along(report_names) ){
      if( report_names[i] %in% names(Report) ){
        dimnames(Report[[report_names[i]]]) = dimnames
      }
    }
    return(Report)
  }

  # Defaults
  if( "treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options) ){
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }else{
    treat_missing_as_zero = FALSE
  }

  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:ncol(Report$D_xt) )
    category_names = "singlespecies"
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$D_xct)[3] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$D_xct)[2] )
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$D_xcy)[3] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$D_xcy)[2] )
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version 8.0.0 through 9.3.0
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$D_gcy)[3] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$D_gcy)[2] )
  }
  if( "D_gct" %in% names(Report)){
    # VAST Version >= 9.4.0
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$D_gct)[3] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$D_gct)[2] )
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$dhat_ktp)[2] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$dhat_ktp)[3] )
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
    if( is.null(year_labels) ) year_labels = paste0( "Time_", 1:dim(Report$dpred_ktp)[2] )
    if( is.null(category_names) ) category_names = paste0( "Category_", 1:dim(Report$dpred_ktp)[3] )
  }

  # Determine year-category pairs with no data
  Num_gct = rep(1,TmbData$n_g) %o% abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3)
  # Drop from maps
  if( treat_missing_as_zero==TRUE ){
    Report$D_gct = ifelse(Num_gct==0, 0, Report$D_gct)
  }else{
    Report$D_gct = ifelse(Num_gct==0, NA, Report$D_gct)
  }

  # Add labels
  Report = add_dimnames( Report = Report,
                         report_names = c("R1_gct","R2_gct","D_gct","Epsilon1_gct","Epsilon2_gct","eta1_gct","eta2_gct"),
                         dimnames = list(NULL, "Category"=category_names, "Time"=year_labels) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Omega1_gc","Omega2_gc"),
                         dimnames = list(NULL, "Category"=category_names) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Xi1_gcp","Xi2_gcp"),
                         dimnames = list(NULL, "Category"=category_names, NULL) )

  # Not used yet
  if( !is.null(Sdreport) ){
    SD_stderr = TMB:::as.list.sdreport( Sdreport, what="Std. Error", report=TRUE )
    SD_estimate = TMB:::as.list.sdreport( Sdreport, what="Estimate", report=TRUE )
  }

  # Check for bad entries
  return( Report )
}

