
#' Plot spatial clusters
#'
#' \code{plot_clusters} plots the clusters on a map
#'
#' Code to plot clusters to look for biogeographic shifts
#'
#' @inheritParams fastcluster::hclust.vector
#' @inheritParams stats::cutree
#'
#' @param var_name Name of object from \code{fit$Report} that is used.  Only
#'        implemented options are \code{"D_gct","Omega1_gc","Omega2_gc","Epsilon1_gct","Epsilon2_gct"}.
#' @param transform_var function to apply to \code{fit$Report[[var_name]]} prior to clustering.
#'        I recommend using log-transform for density, and otherwise using \code{transform_var=identity}.
#' @param method Distance metric. Default \code{method="ward"} is very fast,
#'        but can instead use \code{method="bcdist"} for small problems which calculates
#'        Bray-Curtist dissimilarity using \code{\link[ecodist]{bcdist}} and then applies
#'        Ward clustering
#' @param replace_Inf_with_NA Boolean whether to replace \code{Inf} or \code{-Inf} values
#'        with \code{NA} prior to clustering, as useful sometimes when \code{var_name="D_gct"},
#'        \code{transform_var=log} and replacing nonencounters with zero.
#' @param map_list output from \code{\link{make_map_info output}}
#'
#' @export
plot_clusters <-
function( fit,
          var_name = "D_gct",
          transform_var = log,
          k = 4,
          method = "ward",
          year_labels = fit$year_labels,
          map_list = NULL,
          working_dir = paste0(plotdirgetwd(),"/"),
          file_name = paste0("Class-",var_name),
          replace_Inf_with_NA = TRUE,
          ... ){

  # Informative error
  if( !(var_name %in% c("D_gct","Omega1_gc","Omega2_gc","Epsilon1_gct","Epsilon2_gct")) ){
    stop("Check `var_name`")
  }

  # Change labels
  fit$Report = amend_output( fit,
                             year_labels = year_labels )

  # Make map_list if necessary
  if( missing(map_list) ){
    map_list = make_map_info( Region = fit$settings$Region,
                              spatial_list = fit$spatial_list,
                              Extrapolation_List = fit$extrapolation_list )
  }

  # Extract object
  Y_gct = transform_var( strip_units(fit$Report[[var_name]]) )
  if( length(dim(Y_gct))==2 ){
    Y_gct = Y_gct %o% array(1,dimnames=list("Time"=1))
  }

  #
  if( replace_Inf_with_NA==TRUE ){
    Y_gct = ifelse( abs(Y_gct)==Inf, NA, Y_gct )
  }

  # Change shape
  Y_z = reshape2:::melt.array( data=Y_gct, varnames=names(dimnames(Y_gct)) )
  Y_zc = reshape2::acast(Y_z, formula = Time + Site ~ Category )

  # Warnings
  if( nrow(Y_zc) > 100000 ){
    warning("`plot_clusters` will likely not work due to large size")
    return( list("Y_zc"=Y_zc) )
  }
  if( nrow(Y_zc) > 10000 ) warning("`plot_clusters` will go slowly due to large size")

  # Apply clustering
  if( method == "bcdist" ){
    # Option 1 -- breaks with large sample size
    Dist_zz = ecodist::bcdist(Y_zc)  # dist or ecodist::bcdist
    Hclust = hclust(Dist_zz, method="ward.D2" )
    # Option 2 -- memory issues
    #cluster::agnes( Y_zc )
  }else{
    # Option 3 -- faster and lower memory
    Hclust = fastcluster::hclust.vector( Y_zc, method=method )
  }

  # Cut and convert shape
  Class_z = cutree( Hclust, k = k )
  if( prod(dim(Y_gct)[c(1,3)]) != length(Class_z) ) stop("Check `plot_clusters`")
  Class_gt = array(Class_z, dim=dim(Y_gct)[c(1,3)], dimnames=dimnames(Y_gct)[c(1,3)])

  # Make plot
  plot_variable(
    Y_gt = Class_gt,
    map_list = map_list,
    file_name = file_name,
    working_dir = working_dir,
    #format = format,
    panel_labels = colnames(Class_gt),
    ...
  )

  # Return stuff
  Return = list("Y_zc"=Y_zc, "Class_gt"=Class_gt)
  return( invisible(Return) )
}
