
#' Plot factor-decomposition of covariance
#'
#' \code{plot_factors} plots factor loadings, average spatial factors, and spatio-temporal factors
#'
#' @inheritParams plot_overdispersion
#' @inheritParams summarize_covariance
#' @inheritParams plot_maps
#' @param Year_Set plotting-names for time dimension
#' @param mapdetails_list output from \code{FishStatsUtils::MapDetails_Fn}
#' @param Dim_year Plotting dimension (row,column) for plot of years (default: square with sufficient size for number of years)
#' @param Dim_species Plotting dimension (row,column) for plot of categories (default: square with sufficient size for number of categories)
#' @param plotdir directory for saving plots
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param ... additional arguments passed to \code{plot_maps(.)} when plotting spatio-temporal terms Epsilon

#' @export
plot_factors = function( Report, ParHat, Data, SD=NULL, Year_Set=NULL, category_names=NULL, RotationMethod="PCA",
  mapdetails_list=NULL, Dim_year=NULL, Dim_species=NULL, plotdir=paste0(getwd(),"/"), land_color="grey", zlim=NA, ... ){

  # Extract Options and Options_vec (depends upon version)
  if( all(c("Options","Options_vec") %in% names(Data)) ){
    Options_vec = Data$Options_vec
    Options = Data$Options
  }
  if( "Options_list" %in% names(Data) ){
    Options_vec = Data$Options_list$Options_vec
    Options = Data$Options_list$Options
  }

  #### Deals with backwards compatibility for FieldConfig
  # Converts from 4-vector to 3-by-2 matrix
  if( is.vector(Data$FieldConfig) && length(Data$FieldConfig)==4 ){
    Data$FieldConfig = rbind( matrix(Data$FieldConfig,ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c(-2,-2) )
  }
  # Converts from 3-by-2 matrix to 4-by-2 matrix
  if( is.matrix(Data$FieldConfig) & all(dim(Data$FieldConfig)==c(3,2)) ){
    Data$FieldConfig = rbind( Data$FieldConfig, "Epsilon_time"=c(-3,-3) )
  }
  # Checks for errors
  if( !is.matrix(Data$FieldConfig) || !all(dim(Data$FieldConfig)==c(4,2)) ){
    stop("`FieldConfig` has the wrong dimensions in `plot_factors`")
  }
  # Renames
  dimnames(Data$FieldConfig) = list( c("Omega","Epsilon","Beta","Epsilon_time"), c("Component_1","Component_2") )

  # Fill in missing inputs
  if( "D_xct" %in% names(Report) ){
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xct)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xct)[2]
  }
  if( "D_xcy" %in% names(Report) ){
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
  }
  if( "D_gcy" %in% names(Report) ){
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_gcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_gcy)[2]
  }

  # Dimensions for plotting
  Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
  Dim_year = Dim(length(Year_Set))
  Dim_species = Dim(length(category_names))

  # Extract loadings matrices (more numerically stable than extracting covariances, and then re-creating Cholesky)
  Psi2prime_list = Psiprime_list = Lprime_SE_list = Hinv_list = L_SE_list = Lprime_list = L_list = vector("list", length=8)    # Add names at end so that NULL doesn't interfere

  # Loop through
  for(i in 1:8){

    # Variable names
    Par_name = c("Omega1", "Epsilon1", "Beta1", "EpsilonTime1", "Omega2", "Epsilon2", "Beta2", "EpsilonTime2")[i]

    # Backwards compatible loading of variables and names
    if(Par_name == "Omega1"){ Var_name = "Omegainput1_sf"; Var2_name = "Omegainput1_gf"; L_name = "L_omega1_cf" }
    if(Par_name == "Epsilon1"){ Var_name = "Epsiloninput1_sft"; Var2_name = "Epsiloninput1_gft"; L_name = "L_epsilon1_cf" }
    if(Par_name == "Beta1"){ Var_name = "beta1_ft"; Var2_name = "missing"; L_name = "L_beta1_cf" }
    if(Par_name == "EpsilonTime1"){ Var_name = "Epsiloninput1_sff"; Var2_name = "Epsiloninput1_gff"; L_name = "Ltime_epsilon1_tf" }
    if(Par_name == "Omega2"){ Var_name = "Omegainput2_sf"; Var2_name = "Omegainput2_gf"; L_name = "L_omega2_cf" }
    if(Par_name == "Epsilon2"){ Var_name = "Epsiloninput2_sft"; Var2_name = "Epsiloninput2_gft"; L_name = "L_epsilon2_cf" }
    if(Par_name == "Beta2"){ Var_name = "beta2_ft"; Var2_name = "missing"; L_name = "L_beta2_cf" }
    if(Par_name == "EpsilonTime2"){ Var_name = "Epsiloninput2_sff"; Var2_name = "Epsiloninput2_gff"; L_name = "Ltime_epsilon2_tf" }

    # Continue if component is included
    if( as.vector(Data[["FieldConfig"]])[i] > 0 ){

      # Get loadings matrix
      if( "L_beta1_cf" %in% names(Report) ){
        L_list[[i]] = Report[[L_name]]
      }else{
        L_list[[i]] = calc_cov( L_z=ParHat[[paste0("L_",tolower(Par_name),"_z")]], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c, returntype="loadings_matrix" )
      }
      if( !Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        rownames(L_list[[i]]) = category_names
      }

      # Load SE
      if( class(SD)=="sdreport" ){
        rowindex = grep( paste0("L_",tolower(Par_name),"_z"), rownames(SD$cov.fixed) )
        L_rz = mvtnorm::rmvnorm( n=1e3, mean=ParHat[[paste0("L_",tolower(Par_name),"_z")]], sigma=SD$cov.fixed[rowindex,rowindex] )
        L_rcf = array(NA, dim=c(nrow(L_rz),dim(L_list[[i]])) )
        for( rI in 1:nrow(L_rz) ){
          L_rcf[rI,,] = calc_cov( L_z=L_rz[rI,], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c, returntype="loadings_matrix" )
        }
        Lmean_cf = apply(L_rcf, MARGIN=2:3, FUN=mean)
        Lsd_cf = apply(L_rcf, MARGIN=2:3, FUN=sd)
        L_SE_list[[i]] = Lsd_cf
        rownames(L_SE_list[[i]]) = category_names
      }

      # Get covariance
      if(Var_name%in%names(ParHat)) Psi_sjt = ParHat[[Var_name]]
      if(Var_name%in%names(Report)) Psi_sjt = Report[[Var_name]]
      Psi_gjt = Report[[Var2_name]]
      ## the betas and EpsilonTimeare transposed compared to others so fix that here
      if( Par_name %in% c("Beta1","Beta2") ){
        Psi_sjt <- t(Psi_sjt)
      }
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Psi_sjt = aperm( Psi_sjt, c(1,3,2) )
        Psi_gjt = aperm( Psi_gjt, c(1,3,2) )
      }
      if(is.null(Psi_sjt)){
        stop(paste("Covariance is empty for parameter", Var_name))
      }
      logkappa = unlist(ParHat[c('logkappa1','logkappa2')])[c(1,1,1,1,2,2,2,2)[i]]
      if(Options_vec[8]==0){
        tau = 1 / (exp(logkappa) * sqrt(4*pi));
      }else if(Options_vec[8]==1){
        tau = 1 / sqrt(1-exp(logkappa*2));
      }else stop("Check 'Options_vec[8]' for allowable entries")

      # Rotate stuff
      Var_rot = rotate_factors( L_pj=L_list[[i]], Psi=Psi_sjt/tau, RotationMethod=RotationMethod, testcutoff=1e-4, quiet=TRUE )
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Var_rot$Psi_rot = aperm( Var_rot$Psi_rot, c(1,3,2) )
      }
      Report_tmp = list("D_xct"=Var_rot$Psi_rot, "Epsilon1_sct"=Var_rot$Psi_rot, "Epsilon2_sct"=Var_rot$Psi_rot)
      Lprime_list[[i]] = Var_rot$L_pj_rot
      if( !Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        rownames(Lprime_list[[i]]) = category_names
      }
      Psiprime_list[[i]] = Var_rot$Psi_rot
      Hinv_list[[i]] = Var_rot$Hinv

      # Extract SEs if available
      if( class(SD)=="sdreport" ){
        rowindex = grep( paste0("L_",tolower(Par_name),"_z"), rownames(SD$cov.fixed) )
        L_rz = mvtnorm::rmvnorm( n=1e3, mean=ParHat[[paste0("L_",tolower(Par_name),"_z")]], sigma=SD$cov.fixed[rowindex,rowindex] )
        Lprime_rcf = array(NA, dim=c(nrow(L_rz),dim(L_list[[i]])) )
        for( rI in 1:nrow(L_rz) ){
          tmpmat = calc_cov( L_z=L_rz[rI,], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c, returntype="loadings_matrix" )
          Lprime_rcf[rI,,] = rotate_factors( L_pj=tmpmat, RotationMethod="PCA", testcutoff=1e-4, quiet=TRUE )$L_pj_rot
        }
        Lmean_cf = apply(Lprime_rcf, MARGIN=2:3, FUN=mean)
        Lsd_cf = apply(Lprime_rcf, MARGIN=2:3, FUN=sd)
        Lprime_SE_list[[i]] = Lsd_cf
        rownames(Lprime_SE_list[[i]]) = category_names
      }

      # Extract projected factors is available
      if( !is.null(Psi_gjt) ){
        Var2_rot = rotate_factors( L_pj=L_list[[i]], Psi=Psi_gjt/tau, RotationMethod=RotationMethod, testcutoff=1e-4, quiet=TRUE )
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          Var2_rot$Psi_rot = aperm( Var2_rot$Psi_rot, c(1,3,2) )
        }
        Report2_tmp = list("D_xct"=Var2_rot$Psi_rot, "Epsilon1_sct"=Var2_rot$Psi_rot, "Epsilon2_sct"=Var2_rot$Psi_rot)
        Psi2prime_list[[i]] = Var2_rot$Psi_rot
      }else{
        Report2_tmp = NULL
      }

      # Plot loadings
      Dim_factor = Dim( as.vector(Data[["FieldConfig"]])[i] )
      png( file=paste0(plotdir,"Factor_loadings--",Par_name,".png"), width=Dim_factor[2]*4, height=Dim_factor[1]*4, units="in", res=200 )
        par( mfrow=Dim_factor, mar=c(0,2,2,0) )
        for( cI in 1:as.vector(Data[["FieldConfig"]])[i] ) FishStatsUtils::plot_loadings( L_pj=Var_rot$L_pj_rot, whichfactor=cI )
      dev.off()

      # Plot factors
      if( !is.null(mapdetails_list) & !is.null(Report2_tmp) ){

        # Plot Epsilon
        # Use plot_maps to automatically make one figure per factor
        if( Par_name %in% c("Epsilon1","Epsilon2") ){
          plot_maps(plot_set=c(6,6,NA,6,7,7,NA,7)[i], Report=Report2_tmp, PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]],
            working_dir=plotdir, Year_Set=Year_Set, category_names=paste0("Factor_",1:dim(Var_rot$Psi_rot)[2]),
            legend_x=mapdetails_list[["Legend"]]$x/100, legend_y=mapdetails_list[["Legend"]]$y/100, zlim=zlim, ...)
        }  #

        # Plot Omega
        # Use plot_variable to plot all factors on single figure
        if( Par_name %in% c("Omega1", "Omega2")){
          plot_variable( Y_gt=array(Report_tmp$D_xct[,,1],dim=dim(Report2_tmp$D_xct)[1:2]), map_list=mapdetails_list, working_dir=plotdir,
            panel_labels=paste0("Factor_",1:dim(Var_rot$Psi_rot)[2]), file_name=paste0("Factor_maps--",Par_name) )
        }

        ## Doesn't make sense to make maps of beta factors since they aren't spatial

        # Plot EpsilonTime
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          plot_maps(plot_set=c(6,6,NA,6,7,7,NA,7)[i], Report=Report2_tmp, PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]],
            working_dir=plotdir, category_names=paste0("Factor_",1:dim(Var_rot$Psi_rot)[2]), Panel="Year",
            legend_x=mapdetails_list[["Legend"]]$x/100, legend_y=mapdetails_list[["Legend"]]$y/100, zlim=zlim, ...)
        }  #
      }
    }else{
      Lprime_SE_list[[i]] = L_SE_list[[i]] = L_SE_list[[i]] = Psi2prime_list[[i]] = Psiprime_list[[i]] = Lprime_list[[i]] = L_list[[i]] = "Element not estimated, and therefore empty"
    }
  }

  # Return stuff invisibly
  names(Hinv_list) = names(Psi2prime_list) = names(Psiprime_list) = names(Lprime_SE_list) = names(L_SE_list) = names(Lprime_list) = names(L_list) = c("Omega1", "Epsilon1", "Beta1", "Epsilon1Time1", "Omega2", "Epsilon2", "Beta2", "Epsilon1Time2")
  Return = list("Loadings"=L_list, "Rotated_loadings"=Lprime_list, "Rotated_factors"=Psiprime_list, "Rotated_projected_factors"=Psi2prime_list, "Rotation_matrices"=Hinv_list)
  if( class(SD)=="sdreport" ){
    Return[["Loadings_SE"]] = L_SE_list
    Return[["Rotated_loadings_SE"]] = Lprime_SE_list
  }
  return( invisible(Return) )
}
