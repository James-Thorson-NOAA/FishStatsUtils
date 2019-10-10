

#' Format habitat covariate matrix
#'
#' \code{make_covariates} uses a formula interface to generate covariates
#'
#' This function generates 3D arrays required by \code{VAST::make_data} to incorporate density covariates. The user must supply a data frame of covariate values \code{covariate_data}, with Lat, Lon, and Year columns, and then covariates at each observation \code{i} are matched against the nearest supplied covariate value in that year. Similarly, covariate values at each extrapolation-grid cell in each year, \code{X_gtp}, are matched against the nearest covariate value in that year.
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
make_covariates = function( formula=~0, covariate_data, Year_i, spatial_list, extrapolation_list ){

  # Errors
  if( !is.data.frame(covariate_data) ) stop("Please ensure that `covariate_data` is a data frame")
  if( !all(c("Lat","Lon","Year") %in% names(covariate_data)) ){
    stop( "`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`" )
  }

  #
  sample_data = data.frame( "Year"=Year_i, "Lat"=spatial_list$latlon_i[,'Lat'], "Lon"=spatial_list$latlon_i[,'Lon'] )

  # set of years needed
  Year_Set = min(Year_i):max(Year_i)

  # extract latitude and longitude for extrapolation grid
  latlon_g = spatial_list$latlon_g

  # Create data frame of necessary size
  DF_zp = NULL
  DF_ip = cbind( sample_data, covariate_data[rep(1,nrow(sample_data)),-match(names(sample_data),colnames(covariate_data))] )
  DF_ip[,-match(names(sample_data),colnames(DF_ip))] = NA

  # Loop through data and extrapolation-grid
  for( tI in seq_along(Year_Set) ){

    # Subset to same year
    tmp_covariate_data = covariate_data[ which(Year_Set[tI]==covariate_data[,'Year'] | is.na(covariate_data[,'Year'])), ]
    if( nrow(tmp_covariate_data)==0 ){
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years" )
    }
    #
    Which = which(Year_Set[tI]==sample_data[,'Year'])
    # Do nearest neighbors to define covariates for observations, skipping years without observations
    if( length(Which) > 0 ){
      NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=sample_data[Which,c("Lat","Lon")], k=1 )
      # Add to data-frame
      newcolumns = tmp_covariate_data[ NN$nn.idx[,1], -match(names(sample_data),colnames(tmp_covariate_data)), drop=FALSE ]
      DF_ip[Which, -match(names(sample_data),colnames(DF_ip))] = newcolumns
    }

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=latlon_g[,c("Lat","Lon")], k=1 )
    # Add rows
    newcolumns = tmp_covariate_data[ NN$nn.idx[,1], -match(names(sample_data),colnames(tmp_covariate_data)), drop=FALSE ]
    newrows = cbind("Year"=Year_Set[tI], latlon_g, newcolumns )
    DF_zp = rbind( DF_zp, newrows )
  }
  if( any(is.na(DF_ip)) ) stop("Problem with `DF_ip` in `make_covariates(.)")

  # Convert to dimensions requested
  DF = rbind( DF_ip, DF_zp )
  X = model.matrix( update.formula(formula, ~.+0), data=DF )[,,drop=FALSE]

  # Make X_ip
  X_ip = X[ 1:nrow(DF_ip), , drop=FALSE ]
  X_itp = aperm( X_ip %o% rep(1,length(Year_Set)), perm=c(1,3,2) )

  # Make X_gpt and then permute dimensions
  X_gpt = NULL
  indices = nrow(X_ip)
  for( tI in seq_along(Year_Set) ){
    indices = max(indices) + 1:nrow(latlon_g)
    if( max(indices)>nrow(X) ) stop("Check problem in `make_covariates`")
    X_gpt = abind::abind( X_gpt, X[ indices, , drop=FALSE ], along=3 )
  }
  X_gtp = aperm( X_gpt, perm=c(1,3,2) )

  # warnings
  if( any(apply(X_gtp, MARGIN=2:3, FUN=sd)>10 | apply(X_itp, MARGIN=2:3, FUN=sd)>10) ){
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }

  # return stuff
  Return = list( "X_gtp"=X_gtp, "X_itp"=X_itp )
  return( Return )
}

