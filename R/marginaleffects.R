#' @export
get_coef.fit_model = function(x, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  x$ParHat[[param]]
}

#' @export
get_vcov.fit_model = function(x, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  whichrows = which( names(x$parameter_estimates$par) == param )
  if( is.null(x$parameter_estimates$SD) ){
    out = NULL
  }else{
    out = array(x$parameter_estimates$SD$cov.fixed[whichrows,whichrows],dim=rep(length(whichrows),2))
  }
  return(out)
}

#' @export
set_coef.fit_model = function(x, newpar, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  if( length(newpar) != length(x$ParHat[[param]]) ){
    stop("Check length of 'newpar'")
  }
  x$ParHat[[param]][] <- newpar
  return(x)
}

#' @export
get_predict.fit_model = function(x, newdata, covariate, center=FALSE, ...){
  # update formula (following logic in make_covariates)
  if(covariate=="X1"){
    formula = update.formula(x$X1_formula, ~.+1)
    param = "gamma1_cp"
    data = x$effects$covariate_data_full
  }
  if(covariate=="X2"){
    formula = update.formula(x$X2_formula, ~.+1)
    param = "gamma2_cp"
    data = x$effects$covariate_data_full
  }
  if(covariate=="Q1"){
    formula = update.formula(x$Q1_formula, ~.+1)
    param = "lambda1_k"
    data = x$effects$catchability_data_full
  }
  if(covariate=="Q2"){
    formula = update.formula(x$Q2_formula, ~.+1)
    param = "lambda2_k"
    data = x$effects$catchability_data_full
  }

  # build original model.frame
  frame0 = model.frame( formula=formula, data=data )
  terms0 = terms( frame0 )
  xlevels = .getXlevels( terms0, frame0 )

  # get new design matrix
  terms1 = delete.response( terms0 )
  frame1 = model.frame( terms1, newdata, xlev=xlevels )
  X_ip = model.matrix( terms1, frame1 )

  # Drop intercept (following logic in make_covariates)
  Columns_to_keep = which( attr(X_ip,"assign") != 0 )
  X_ip = X_ip[,Columns_to_keep,drop=FALSE]

  # Multiply and center
  gamma_cp = get_coef(x, covariate=covariate)
  yhat_ic = X_ip %*% t(gamma_cp)
  if(center==TRUE) yhat_ic = yhat_ic - outer(rep(1,nrow(yhat_ic)),colMeans(yhat_ic))

  # Return
  if( covariate %in% c("X1","X2") ){
    out = expand.grid( rowid=seq_along(yhat_ic[,1]), category=x$category_names )
    out$predicted = as.vector(yhat_ic)
  }
  if( covariate %in% c("Q1","Q2") ){
    out = data.frame( rowid=seq_along(yhat_ic[,1]), predicted=yhat_ic[,1] )
  }
  return(out)
}
