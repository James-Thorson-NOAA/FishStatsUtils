
#' Plot spatial clusters
#'
#' \code{plot_clusters} plots the clusters on a map
#'
#' Code to plot clusters to look for biogeographic shifts
#'
#' @param method Distance metric. Default \code{method="ward"} is very fast,
#'        but can instead use \code{method="bcdist"} for small problems
#'        \code{\link[ecodist]{bcdist}}
#'
#' @export
plot_clusters <-
function( fit,
          map_list = NULL,
          var_name = "D_gct",
          transform_var = log,
          k = 4,
          method = "ward",
          working_dir = paste0(getwd(),"/"),
          format = "points",
          file_name = paste0("Class-",var_name),
          ... ){

  # Informative error
  if( !(var_name %in% c("D_gct","Omega1_gc","Omega2_gc","Epsilon1_gct","Epsilon2_gct")) ){
    stop("Check `var_name`")
  }

  # Change labels
  fit$Report = amend_output(fit)

  # Make map_list if necessary
  if( missing(map_list) ){
    map_list = make_map_info( Region = fit$settings$Region,
                              spatial_list = fit$spatial_list,
                              Extrapolation_List = fit$extrapolation_list )
  }

  # Extract object
  Y_gct = transform_var( fit$Report[[var_name]] )
  if( length(dim(Y_gct))==2 ){
    Y_gct = Y_gct %o% array(1,dimnames=list("Time"=1))
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
    working_dir = run_dir,
    #format = format,
    panel_labels = colnames(Class_gt),
    ...
  )

  # Return stuff
  Return = list("Y_zc"=Y_zc, "Class_gt"=Class_gt)
  return( invisible(Return) )
}
