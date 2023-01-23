#' @method get_coef fit_model
#' @export
get_coef.fit_model = function(model, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  model$ParHat[[param]]
}

#' @method get_vcov fit_model
#' @export
get_vcov.fit_model = function(model, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  whichrows = which( names(model$parameter_estimates$par) == param )
  if( is.null(model$parameter_estimates$SD) ){
    out = NULL
  }else{
    out = array(model$parameter_estimates$SD$cov.fixed[whichrows,whichrows],dim=rep(length(whichrows),2))
  }
  return(out)
}

#' @method set_coef fit_model
#' @export
set_coef.fit_model = function(model, newpar, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  #if( length(newpar) != length(model$ParHat[[param]]) ){
  #  stop("Check length of 'newpar'")
  #}
  # Make substitution
  if( param%in%names(fit$tmb_list$Map) & any(is.na(fit$tmb_list$Map[[param]])) ){
    newvec = newpar[fit$tmb_list$Map[[param]]]
    model$ParHat[[param]][] <- ifelse( is.na(newvec), model$ParHat[[param]][], newvec )
  }else{
    model$ParHat[[param]][] <- newpar
  }
  return(model)
}

#' @method get_predict fit_model
#' @export
get_predict.fit_model = function(model, newdata, covariate, center=FALSE, ...){
  # update formula (following logic in make_covariates)
  if(covariate=="X1"){
    formula = update.formula(model$X1_formula, ~.+1)
    param = "gamma1_cp"
    data = model$effects$covariate_data_full
  }
  if(covariate=="X2"){
    formula = update.formula(model$X2_formula, ~.+1)
    param = "gamma2_cp"
    data = model$effects$covariate_data_full
  }
  if(covariate=="Q1"){
    formula = update.formula(model$Q1_formula, ~.+1)
    param = "lambda1_k"
    data = model$effects$catchability_data_full
  }
  if(covariate=="Q2"){
    formula = update.formula(model$Q2_formula, ~.+1)
    param = "lambda2_k"
    data = model$effects$catchability_data_full
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
  gamma_cp = get_coef(model, covariate=covariate)
  yhat_ic = X_ip %*% t(gamma_cp)
  if(center==TRUE) yhat_ic = yhat_ic - outer(rep(1,nrow(yhat_ic)),colMeans(yhat_ic))

  # Return
  if( covariate %in% c("X1","X2") ){
    out = expand.grid( rowid=seq_along(yhat_ic[,1]), category=model$category_names )
    out$predicted = as.vector(yhat_ic)
  }
  if( covariate %in% c("Q1","Q2") ){
    out = data.frame( rowid=seq_along(yhat_ic[,1]), predicted=yhat_ic[,1] )
  }
  return(out)
}
