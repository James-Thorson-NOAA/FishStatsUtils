

#' Format habitat covariate matrix
#'
#' \code{make_covariates} uses a formula interface to generate covariates
#'
#' This function generates 3D arrays \code{Cov_gtp} and \code{Cov_itp} required by \code{VAST::make_data} to incorporate density covariates.
#' The user must supply a data frame \code{covariate_data} of covariate values, with columns named Lat, Lon, and Year,
#' as well as values for all covariates as additional named columns.
#' This data frame is then used as a "look-up table", and is matched against variables listed in \code{formula}.
#'
#' Specifically, for every observation \code{i} at location \code{Lat_i[i]} and \code{Lon_i[i]} in year \code{t_i[t]}, the nearest
#' Lat-Lon observation in that year is identified in \code{covariate_data}, and covariate
#' values in that row of \code{covariate_data} are assigned to observation \code{i}.
#' Similarly, for every extrapolation-grid cell \code{g} at location \code{spatial_list$latlon_g[i,]} in each year,
#' the nearest row of \code{covariate_data} in that year
#' is used to assign covariate values. \code{make_covariates} then formats these covariate values appropriately and returns them.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Similar specification to \code{\link{stats::lm}}
#' @param covariate_data data frame of covariate values with columns \code{Lat}, \code{Lon}, and \code{Year}, and other columns matching names in \code{formula}; \code{Year=NA} can be used for covariates that do not change among years (e.g., depth)
#'
#' @return Tagged list of useful output
#' \describe{
#'   \item{Cov_gtp}{3-dimensional array for use in \code{VAST::make_data}}
#'   \item{Cov_itp}{3-dimensional array for use in \code{VAST::make_data}}
#' }

#' @export
make_covariates = function( formula, covariate_data, Year_i, spatial_list, extrapolation_list ){

  # Errors
  if( !is.data.frame(covariate_data) ) stop("Please ensure that `covariate_data` is a data frame")
  if( !all(c("Lat","Lon","Year") %in% names(covariate_data)) ){
    stop( "`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`" )
  }

  # transform data inputs
  sample_data = data.frame( "Year"=Year_i, "Lat"=spatial_list$latlon_i[,'Lat'], "Lon"=spatial_list$latlon_i[,'Lon'] )
  covariate_names = setdiff( names(covariate_data), names(sample_data) )

  # set of years needed
  Year_Set = min(Year_i):max(Year_i)

  # extract latitude and longitude for extrapolation grid
  latlon_g = spatial_list$latlon_g

  # Create data frame of necessary size
  DF_zp = NULL
  #DF_ip = cbind( sample_data, covariate_data[rep(1,nrow(sample_data)),covariate_names] )
  DF_ip = data.frame( sample_data, covariate_data[rep(1,nrow(sample_data)),covariate_names] )
  colnames(DF_ip) = c( names(sample_data) ,covariate_names )
  #DF_ip[,covariate_names] = NA

  # Loop through data and extrapolation-grid
  for( tI in seq_along(Year_Set) ){

    # Subset to same year
    tmp_covariate_data = covariate_data[ which(Year_Set[tI]==covariate_data[,'Year'] | is.na(covariate_data[,'Year'])), , drop=FALSE]
    if( nrow(tmp_covariate_data)==0 ){
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years" )
    }
    #
    Which = which(Year_Set[tI]==sample_data[,'Year'])
    # Do nearest neighbors to define covariates for observations, skipping years without observations
    if( length(Which) > 0 ){
      NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=sample_data[Which,c("Lat","Lon")], k=1 )
      # Add to data-frame
      nearest_covariates = tmp_covariate_data[ NN$nn.idx[,1], covariate_names, drop=FALSE ]
      DF_ip[Which, covariate_names] = nearest_covariates
    }

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=latlon_g[,c("Lat","Lon")], k=1 )
    # Add rows
    nearest_covariates = tmp_covariate_data[ NN$nn.idx[,1], covariate_names, drop=FALSE ]
    newrows = cbind("Year"=Year_Set[tI], latlon_g, nearest_covariates )
    DF_zp = rbind( DF_zp, newrows )
  }

  # Convert to dimensions requested
  DF = rbind( DF_ip, DF_zp )

  # Make model.matrix
    # To ensure identifiability given betas (intercepts), add intercept to formula
    # and then remove that term from model.matrix. This will not fix identifiability
    # issues arising when both conditions are met:
    # factor(Year) has an interaction with another factor, and
    # betas vary among years (are not constant)
  Model_matrix = model.matrix( update.formula(formula, ~.+1), data=DF )
  Columns_to_keep = which( attr(Model_matrix,"assign") != 0 )
  coefficient_names = attr(Model_matrix,"dimnames")[[2]][Columns_to_keep]
  X = Model_matrix[,Columns_to_keep,drop=FALSE]

  # Make X_ip
  X_ip = X[ 1:nrow(DF_ip), , drop=FALSE ]
  X_itp = aperm( X_ip %o% rep(1,length(Year_Set)), perm=c(1,3,2) )
  if( any(is.na(X_itp)) ) stop("Problem with `X_itp` in `make_covariates(.)")

  # Make X_gpt and then permute dimensions
  X_gpt = NULL
  indices = nrow(X_ip)
  for( tI in seq_along(Year_Set) ){
    indices = max(indices) + 1:nrow(latlon_g)
    if( max(indices)>nrow(X) ) stop("Check problem in `make_covariates`")
    X_gpt = abind::abind( X_gpt, X[ indices, , drop=FALSE ], along=3 )
  }
  X_gtp = aperm( X_gpt, perm=c(1,3,2) )
  if( any(is.na(X_gtp)) ) stop("Problem with `X_gtp` in `make_covariates(.)")

  # warnings
  if( any(apply(X_gtp, MARGIN=2:3, FUN=sd)>10 | apply(X_itp, MARGIN=2:3, FUN=sd)>10) ){
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }

  # return stuff
  Return = list( "X_gtp"=X_gtp, "X_itp"=X_itp, "covariate_names"=covariate_names,
    "coefficient_names"=coefficient_names )
  return( Return )
}

