

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
#' If all covariates as "static" (not changing among years),
#' then set Year = NA to cause values to be duplicated internally for all values of Year.
#' If using a mix of static and dynamic covariates,
#' then duplicate rows for static covariates for every value of Year
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
make_covariates <-
function( formula,
          covariate_data,
          contrasts,
          Year_i,
          spatial_list ){

  # Check for bad entries
  if( !is.data.frame(covariate_data) ) stop("Please ensure that `covariate_data` is a data frame")
  if( !all(c("Lat","Lon","Year") %in% names(covariate_data)) ){
    stop( "`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`" )
  }

  # set of years needed
  Year_Set = min(Year_i):max(Year_i)

  # Make `covariate_df` by expanding for rows with Year=NA
  covariate_df = covariate_data[ which(!is.na(covariate_data[,'Year'])), ]
  for( tI in seq_along(Year_Set) ){
    newrows = covariate_data[ which(is.na(covariate_data[,'Year'])), ]
    newrows[,"Year"] = rep( Year_Set[tI], nrow(newrows) )
    covariate_df = rbind( covariate_df, newrows )
  }

  # Make model.matrix
    # To ensure identifiability given betas (intercepts), add intercept to formula
    # and then remove that term from model.matrix. This will not fix identifiability
    # issues arising when both conditions are met:
    # factor(Year) has an interaction with another factor, and
    # betas vary among years (are not constant)
  Model_matrix = model.matrix( update.formula(formula, ~.+1), data=covariate_df,  contrasts.arg=contrasts )
  Columns_to_keep = which( attr(Model_matrix,"assign") != 0 )
  coefficient_names = attr(Model_matrix,"dimnames")[[2]][Columns_to_keep]
  X = Model_matrix[,Columns_to_keep,drop=FALSE]
  dimnames(X) = list(NULL, coefficient_names)

  # transform data inputs
  sample_i = data.frame( "Year"=Year_i, "Lat"=spatial_list$latlon_i[,'Lat'], "Lon"=spatial_list$latlon_i[,'Lon'] )
  #covariate_names = setdiff( names(covariate_data), names(sample_data) )

  # extract latitude and longitude for extrapolation grid
  latlon_g = spatial_list$latlon_g

  # Create data frame of necessary size
  X_gtp = array( NA, dim=c(nrow(latlon_g),length(Year_Set),ncol(X)), dimnames=list(NULL,Year_Set,colnames(X)) )
  X_ip = array( NA, dim=c(nrow(sample_i),ncol(X)), dimnames=list(NULL,colnames(X)) )

  # Loop through data and extrapolation-grid
  for( tI in seq_along(Year_Set) ){

    # Subset to same year
    tmp_covariate_df = covariate_df[ which(Year_Set[tI]==covariate_df[,'Year']), , drop=FALSE]
    tmp_X = X[ which(Year_Set[tI]==covariate_df[,'Year']), , drop=FALSE]
    if( nrow(tmp_covariate_df)==0 ){
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years" )
    }

    # Fill in values in X_ip
    Which = which(Year_Set[tI]==sample_i[,'Year'])
    # Do nearest neighbors to define covariates for observations, skipping years without observations
    if( length(Which) > 0 ){
      NN = RANN::nn2( data=tmp_covariate_df[,c("Lat","Lon")], query=sample_i[Which,c("Lat","Lon")], k=1 )
      # Fill in values
      X_ip[Which, ] = tmp_X[ NN$nn.idx[,1], , drop=FALSE]
    }

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_df[,c("Lat","Lon")], query=latlon_g[,c("Lat","Lon")], k=1 )
    # Add rows
    X_gtp[,tI,] = tmp_X[ NN$nn.idx[,1], , drop=FALSE ]
  }

  # Make X_itp
  X_itp = aperm( X_ip %o% rep(1,length(Year_Set)), perm=c(1,3,2) )

  # Check for obvious problems
  if( any(is.na(X_itp)) ) stop("Problem with `X_itp` in `make_covariates(.)")
  if( any(is.na(X_gtp)) ) stop("Problem with `X_gtp` in `make_covariates(.)")

  # warnings
  if( any(apply(X_gtp, MARGIN=2:3, FUN=sd)>10 | apply(X_itp, MARGIN=2:3, FUN=sd)>10) ){
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }

  # return stuff
  Return = list( "X_gtp"=X_gtp, "X_itp"=X_itp, # "covariate_names"=covariate_names,
    "coefficient_names"=coefficient_names )
  return( Return )
}

