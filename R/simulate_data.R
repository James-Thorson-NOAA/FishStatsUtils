
#' @title
#' Simulate new dynamics and sampling data
#'
#' @description
#' \code{simulate_data} conducts a parametric bootstrap to simulate new data and potentially simulate new population dynamics and associated variables
#'
#' @param fit output form \code{fit_model(.)}
#' @param type integer stating what type of simulation to use.  \code{type=1} simulates new data conditional upon estimated fixed and random effects.  \code{type=2} simulates new random effects conditional upon fixed effects, and new data conditional upon both.  \code{type=3} simulates new fixed and random effects from the joint precision matrix, and new data conditional upon these values.
#' @param random_seed integer passed to \code{\link[base]{set.seed}} whenever \code{type=3}, where the default value \code{random_seed=NULL} resets the random-number seed. Argument no effect when \code{type!=3} because TMB has no interface for setting the random-number seed in C++
#'

#' @return Report object containing new data and population variables including
#' \describe{
#'   \item{b_i}{New simulated data}
#'   \item{D_gcy}{Density for each grid cell g, category c, and year y}
#'   \item{Index_cyl}{Index of abundance for each category c, year y, and stratum l}
#' }

#' @export
simulate_data = function( fit, type=1, random_seed=NULL ){

  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
    set.seed( random_seed )
    z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L = Matrix::Cholesky(prec, super=TRUE)
    z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z = as.matrix(z)
    return(mu + z)
  }

  # Extract stuff
  Obj = fit$tmb_list$Obj
  if( !is.null(random_seed) & type!=3 ){
    stop("Specifying argument `random_seed` in `simulate_data(.) only works when `type=3`")
  }

  # Simulate conditional upon fixed and random effect estimates
  if( type==1 ){
    Obj$env$data$Options_list$Options['simulate_random_effects'] = FALSE
    Return = Obj$simulate( complete=TRUE )
  }

  # Simulate new random effects and data
  if( type==2 ){
    Obj$env$data$Options_list$Options['simulate_random_effects'] = TRUE
    Return = Obj$simulate( complete=TRUE )
  }

  # Simulate from predictive distribution, and then new data
  if( type==3 ){
    # Informative error messages
    if( !("jointPrecision" %in% names(fit$parameter_estimates$SD)) ){
      stop("jointPrecision not present in fit$parameter_estimates$SD; please re-run with `getJointPrecision=TRUE`")
    }

    # Sample from joint distribution
    newpar = rmvnorm_prec( mu=Obj$env$last.par.best, prec=fit$parameter_estimates$SD$jointPrecision, n.sims=1, random_seed=random_seed )[,1]

    # Simulate
    Obj$env$data$Options_list$Options['simulate_random_effects'] = FALSE
    Return = Obj$simulate( par=newpar, complete=TRUE )
  }

  # Return
  return( Return )
}


