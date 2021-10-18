

#' Amend output from VAST for user convenience
#'
#' \code{amend_output} add labels, units, and performs logical operations,
#' e.g., adds zeros as needed, to simplify user and downstream interpretation.
#'
#' @export
amend_output <-
function( TmbData,
          Report,
          extrapolation_list,
          Sdreport = NULL,
          year_labels = NULL,
          category_names = NULL,
          strata_names = NULL ){

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

  # Local function
  process_labels = function( labels, prefix, length ){
    if(is.null(labels)){
      labels = paste0( prefix, "_", seq_len(length) )
    }else{
      if(length(labels)!=length) stop("Check labels")
    }
    return(labels)
  }

  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    year_labels = process_labels( year_labels, "Time", ncol(Report$D_xt) )
    category_names = "singlespecies"
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_xct)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_xct)[2] )
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_xcy)[3] )
    category_names = process_labels( category_names, "Category_", dim(Report$D_xcy)[2] )
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version 8.0.0 through 9.3.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_gcy)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_gcy)[2] )
  }
  if( "D_gct" %in% names(Report)){
    # VAST Version >= 9.4.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_gct)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_gct)[2] )
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    year_labels = process_labels( year_labels, "Time", dim(Report$dhat_ktp)[2] )
    category_names = process_labels( category_names, "Category", dim(Report$dhat_ktp)[3] )
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
    year_labels = process_labels( year_labels, "Time", dim(Report$dpred_ktp)[2] )
    category_names = process_labels( category_names, "Category", dim(Report$dpred_ktp)[3] )
  }
  strata_names = process_labels( strata_names, "Stratum", dim(Report$Index_ctl)[3] )

  # Determine year-category pairs with no data
  Num_gct = rep(1,TmbData$n_g) %o% abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3)
  # Drop from maps
  if( treat_missing_as_zero==TRUE ){
    Report$D_gct = ifelse(Num_gct==0, 0, Report$D_gct)
  }else{
    Report$D_gct = ifelse(Num_gct==0, NA, Report$D_gct)
  }

  # Add labels for all variables plotted using `plot_maps`
  Report = add_dimnames( Report = Report,
                         report_names = c("P1_gct","P2_gct","R1_gct","R2_gct","D_gct","Epsilon1_gct","Epsilon2_gct","eta1_gct","eta2_gct"),
                         dimnames = list(NULL, "Category"=category_names, "Time"=year_labels) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Omega1_gc","Omega2_gc"),
                         dimnames = list(NULL, "Category"=category_names) )
  Report = add_dimnames( Report = Report,
                         report_names = "Xi1_gcp",
                         dimnames = list(NULL, "Category"=category_names, "Covariate"=colnames(TmbData$X1_ip)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Xi2_gcp",
                         dimnames = list(NULL, "Category"=category_names, "Covariate"=colnames(TmbData$X2_ip)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Phi1_gk",
                         dimnames = list(NULL, "Covariate"=colnames(TmbData$Q1_ik)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Phi2_gk",
                         dimnames = list(NULL, "Covariate"=colnames(TmbData$Q2_ik)) )

  # Add labels for other useful variables
  Report = add_dimnames( Report = Report,
                         report_names = c("Index_ctl","effective_area_ctl","mean_D_ctl"),
                         dimnames = list("Category"=category_names, "Time"=year_labels, "Stratum"=strata_names) )

  # Modify Sdreport
  if( !is.null(Sdreport) ){
    # See plot_biomass_index for how to efficiently extra SEs
  }

  # Add units
  units(Report$Index_ctl) = units(TmbData$b_i / TmbData$a_i * extrapolation_list$Area_km2[1])
  units(Report$D_gct) = units(TmbData$b_i)
  units(Report$mean_D_ctl) = units(TmbData$b_i)

  # Check for bad entries
  return( Report )
}

