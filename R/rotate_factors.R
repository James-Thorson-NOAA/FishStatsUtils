
#' Rotate results
#'
#' \code{rotate_factors} rotates results from a factor model
#'
#' @param Cov_jj Covariance calculated from loadings matrix
#' @param L_pj Loadings matrix for `p` categories and `j` factors (calculated from \code{Cov_jj} if it is provided)
#' @param Psi_sjt Array of factors (1st dimension: spatial knots;  2nd dimension: factors;  3rd dimension:  time)
#' @param RotationMethod Method used for rotation, Options: "PCA" (recommended) or "Varimax"
#' @param testcutoff tolerance for numerical rounding when confirming that rotation doesn't effect results

#' @return tagged list of outputs
#' \describe{
#'   \item{L_pj_rot}{Loadings matrix after rotation}
#'   \item{Psi_rot}{Factors after rotation}
#'   \item{Hinv}{Object used for rotation}
#'   \item{L_pj}{Loadings matrix}
#' }

#' @export
rotate_factors = function( Cov_jj=NULL, L_pj=NULL, Psi_sjt=NULL, RotationMethod="PCA", testcutoff=1e-10,
  quiet=FALSE ){

  # If missing time, add a third dimension
  if( length(dim(Psi_sjt))==2 ){
    Psi_sjt = array( Psi_sjt, dim=c(dim(Psi_sjt),1) )
  }

  # Local functions
  approx_equal = function(m1,m2,denominator=mean(m1+m2),d=1e-10) (2*abs(m1-m2)/denominator) < d
  trunc_machineprec = function(n) ifelse(n<1e-10,0,n)
  Nknots = dim(Psi_sjt)[1]
  Nfactors = ncol(L_pj)
  Nyears = nrow(L_pj)

  # Optional inputs
  if( !is.null(Cov_jj) ){
    if(quiet==FALSE) message( "Re-calculating L_pj from Cov_jj")
    if( sum(eigen(Cov_jj)$values>testcutoff)<ncol(Cov_jj) ){
      stop("Calculating L_pj from Cov_jj in 'Rotate_Fn' only works well when Cov_jj is full rank")
    }
    L_pj = t(chol(Cov_jj))[,1:Nfactors]
  }else{
    if( !is.null(L_pj) ){
      if(quiet==FALSE) message("Using L_pj for loadings matrix")
    }else{
      stop( "Must provide either L_pj or Cov_jj" )
    }
  }

  # Varimax
  if( RotationMethod=="Varimax" ){
    Hinv = varimax( L_pj, normalize=FALSE )
    H = solve(Hinv$rotmat)
    L_pj_rot = L_pj %*% Hinv$rotmat
    if( !is.null(Psi_sjt) ){
      Psi_rot = array(NA, dim=dim(Psi_sjt))
      # My new factors
      for( n in 1:Nknots ){
        Psi_rot[n,,] = H %*% Psi_sjt[n,,]
      }
    }
  }

  # PCA
  if( RotationMethod=="PCA" ){
    Cov_tmp = L_pj%*%t(L_pj)
    Cov_tmp = 0.5*Cov_tmp + 0.5*t(Cov_tmp) # Avoid numerical issues with complex eigen-decomposition due to numerical underflow
    Eigen = eigen(Cov_tmp)
    Eigen$values_proportion = Eigen$values / sum(Eigen$values)
    Eigen$values_cumulative_proportion = cumsum(Eigen$values) / sum(Eigen$values)
    # Check decomposition
    #all(approx_equal( Eigen$vectors%*%diag(Eigen$values)%*%t(Eigen$vectors), L_pj%*%t(L_pj)))
    # My attempt at new loadings matrix
    L_pj_rot = (Eigen$vectors%*%diag(sqrt(trunc_machineprec(Eigen$values))))[,1:Nfactors,drop=FALSE]
    rownames(L_pj_rot) = rownames(L_pj)
    # My new factors
    H = corpcor::pseudoinverse(L_pj_rot) %*% L_pj
    Hinv = list("rotmat"=solve(H))
    if( !is.null(Psi_sjt) ){
      Psi_rot = array(NA, dim=dim(Psi_sjt))
      for( n in 1:Nknots ){
        Psi_rot[n,,] = H %*% Psi_sjt[n,,]
      }
    }
  }

  # Flip around
  for( j in 1:dim(L_pj_rot)[2] ){
    if( !is.null(Psi_sjt) ){
      Psi_rot[,j,] = Psi_rot[,j,] * sign(sum(L_pj_rot[,j]))
    }
    L_pj_rot[,j] = L_pj_rot[,j] * sign(sum(L_pj_rot[,j]))
  }

  # Check for errors
  # Check covariance matrix
    # Should be identical for rotated and unrotated
  if( !is.na(testcutoff) ){
    if( !all(approx_equal(L_pj%*%t(L_pj),L_pj_rot%*%t(L_pj_rot), d=testcutoff)) ) stop("Covariance matrix is changed by rotation")
    # Check linear predictor
      # Should give identical predictions as unrotated
    if( !is.null(Psi_sjt) ){
      for(i in 1:dim(Psi_sjt)[[1]]){
      for(j in 1:dim(Psi_sjt)[[3]]){
        MaxDiff = max(L_pj%*%Psi_sjt[i,,j] - L_pj_rot%*%Psi_rot[i,,j])
        if( !all(approx_equal(L_pj%*%Psi_sjt[i,,j],L_pj_rot%*%Psi_rot[i,,j], d=testcutoff, denominator=1)) ) stop(paste0("Linear predictor is wrong for site ",i," and time ",j," with difference ",MaxDiff))
      }}
    }
    # Check rotation matrix
      # Should be orthogonal (R %*% transpose = identity matrix) with determinant one
      # Doesn't have det(R) = 1; determinant(Hinv$rotmat)!=1 ||
    Diag = Hinv$rotmat %*% t(Hinv$rotmat)
    diag(Diag) = ifelse( diag(Diag)==0,1,diag(Diag) )
    if( !all(approx_equal(Diag,diag(Nfactors), d=testcutoff)) ) stop("Rotation matrix is not a rotation")
  }

  # Return stuff
  Return = list( "L_pj_rot"=L_pj_rot, "Hinv"=Hinv, "L_pj"=L_pj )
  if(RotationMethod=="PCA") Return[["Eigen"]] = Eigen
  if( !is.null(Psi_sjt) ){
    Return[["Psi_rot"]] = Psi_rot
  }
  return( Return )
}
