#' @title
#' Plot standard maps
#'
#' @description
#' \code{plot_maps} plots a standard set of diagnostic maps
#'
#' @inheritParams FishStatsUtils::plot_variable

#' @param plot_set integer-vector defining plots to create
#' \describe{
#'   \item{plot_set=1}{Probability of encounter/non-encounter}
#'   \item{plot_set=2}{Log-expected positive catch rate}
#'   \item{plot_set=3}{Log-predicted density (product of encounter probability and positive catch rates)}
#'   \item{plot_set=4}{Log-positive catch rates (rescaled)}
#'   \item{plot_set=5}{Log-predicted density (rescaled)}
#'   \item{plot_set=6}{Spatio-temporal variation in encounter probability}
#'   \item{plot_set=7}{Spatio-temporal variation in log-positive catch rates}
#'   \item{plot_set=8}{Linear predictor for encounter probability}
#'   \item{plot_set=9}{Linear predictor for positive catch rates}
#'   \item{plot_set=10}{Coefficient of variation for predicted density (available only if \code{Data_Fn(...,Options=c('SD_site_logdensity'=1,...))}}
#'   \item{plot_set=11}{Covariates that are included in the model}
#'   \item{plot_set=12}{Total biomass across all categories (only useful in a multivariate model)}
#'   \item{plot_set=13}{Covariate effects on encounter probability}
#'   \item{plot_set=14}{Covariate effects on positive catch rates}
#' }
#' @param Report tagged list of outputs from TMB model via \code{Obj$report()}
#' @param Sdreport Standard deviation outputs from TMB model via \code{sdreport(Obj)}
#' @param Panel Whether to plot years for a given category (\code{Panel="Category"}) or categories for a given year ((\code{Panel="Year"})  in each panel figure
#' @param MapSizeRatio Default size for each panel
#' @param Year_Set Year names for labeling panels
#' @param Years2Include integer vector, specifying positions of \code{Year_Set} for plotting (used to avoid plotting years with no data, etc.)
#' @param category_names character vector specifying names for different categories (only used for R package \code{VAST})
#' @param ... arguments passed to \code{FishStatsUtils::plot_variable}
#'
#' @return Mat_xt a matrix (rows: modeled knots; column: modeled year) for plotted output of last element of \code{plot_set}
#'

#' @export
plot_maps <-
function(plot_set=3, Report, PlotDF, Sdreport=NULL, TmbData=NULL, projargs='+proj=longlat',
         Panel="Category", Year_Set=NULL, Years2Include=NULL, category_names=NULL, quiet=FALSE,
         working_dir=paste0(getwd(),"/"), MapSizeRatio, n_cells, ...){

  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    if( is.null(Year_Set) ) Year_Set = 1:ncol(Report$D_xt)
    if( is.null(Years2Include) ) Years2Include = 1:ncol(Report$D_xt)
    category_names = "singlespecies"
    Ncategories = length(category_names)
    Nyears = dim(Report$D_xt)[2]
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xct)[3]
    if( is.null(Years2Include) ) Years2Include = 1:dim(Report$D_xct)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xct)[2]
    Ncategories = dim(Report$D_xct)[2]
    Nyears = dim(Report$D_xct)[3]
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xcy)[3]
    if( is.null(Years2Include) ) Years2Include = 1:dim(Report$D_xcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
    Ncategories = dim(Report$D_xcy)[2]
    Nyears = dim(Report$D_xcy)[3]
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version >= 8.0.0
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_gcy)[3]
    if( is.null(Years2Include) ) Years2Include = 1:dim(Report$D_gcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_gcy)[2]
    Ncategories = dim(Report$D_gcy)[2]
    Nyears = dim(Report$D_gcy)[3]
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$dhat_ktp)[2]
    if( is.null(Years2Include) ) Years2Include = 1:dim(Report$dhat_ktp)[2]
    if( is.null(category_names) ) category_names = 1:dim(Report$dhat_ktp)[3]
    Ncategories = dim(Report$dhat_ktp)[3]
    Nyears = dim(Report$dhat_ktp)[2]
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$dpred_ktp)[2]
    if( is.null(Years2Include) ) Years2Include = 1:dim(Report$dpred_ktp)[2]
    if( is.null(category_names) ) category_names = 1:dim(Report$dpred_ktp)[3]
    Ncategories = dim(Report$dpred_ktp)[3]
    Nyears = dim(Report$dpred_ktp)[2]
  }
  if( missing(MapSizeRatio) ){
    MapSizeRatio = c(3, 3)
  }

  # Errors
  if( Nyears != length(Year_Set) ){
    stop("Problem with `Year_Set`")
  }
  if( Ncategories != length(category_names) ){
    stop("Problem with `category_names`")
  }

  # Loop through plots
  Return = NULL
  for(plot_num in plot_set){

    # Extract elements
    Array_xct = NULL
    plot_code <- c("encounter_prob", "pos_catch", "density", "", "", "epsilon_1", "epsilon_2", "linear_predictor_1", "linear_predictor_2", "density_CV", "covariates", "total_density", "covariate_effects_1", "covariate_effects_2", "omega_1", "omega_2")[plot_num]
    #if( missing(textmargin) ){
    #  textmargin <- c("Probability of encounter", "Density, ln(kg. per square km.)", "Density, ln(kg. per square km.)", "", "", "", "", "", "", "CV of density (dimensionless)", "Covariate value", "Density, ln(kg. per square km.)", "", "")[plot_num]
    #}

    # Extract matrix to plot
    if(plot_num==1){
      # Presence/absence ("Pres")
      if( quiet==FALSE ) message(" # Plotting presence/absense maps")
      if("D_xt"%in%names(Report)) Array_xct = Report$R1_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$R1_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$R1_xcy
      #if("D_gcy"%in%names(Report)) Array_xct = Report$R1_gcy
      if("D_gcy"%in%names(Report)) Array_xct = Report$R1_gcy
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report))) stop("Not implemented for SpatialVAM")
      message( "`plot_num=1` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells" )
    }
    if(plot_num==2){
      # Positive values ("Pos")
      if( quiet==FALSE ) message(" # Plotting positive catch rate maps")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$R2_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(Report$R2_xct)
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$R2_xcy)
      #if("D_gcy"%in%names(Report)) Array_xct = log(Report$R2_gcy)
      if("D_gcy"%in%names(Report)) Array_xct = Report$R2_gcy
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
      message( "`plot_num=2` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells" )
    }
    if(plot_num==3){
      # Density ("Dens")
      if( quiet==FALSE ) message(" # Plotting density maps")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(Report$D_xct)
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$D_xcy)
      if("D_gcy"%in%names(Report)) Array_xct = log(Report$D_gcy)
      if("dhat_ktp" %in% names(Report)) Array_xct = aperm(Report$dhat_ktp[,,cI],c(1,3,2))
      if("dpred_ktp" %in% names(Report)) Array_xct = aperm(Report$dpred_ktp[,,cI],c(1,3,2))
    }
    if(plot_num==4){
      # Positive values rescaled ("Pos_Rescaled")
      stop( "`plot_num=4` is deprecated")
    }
    if(plot_num==5){
      # Density rescaled ("Dens_Rescaled")
      stop( "`plot_num=5` is deprecated")
    }
    if(plot_num==6){
      # Epsilon for presence/absence ("Eps_Pres")
      if( quiet==FALSE ) message(" # Plotting spatio-temporal effects (Epsilon) in 1st linear predictor")
      if("D_xt"%in%names(Report)) Array_xct = Report$Epsilon1_st
      if("D_xct"%in%names(Report)) Array_xct = Report$Epsilon1_sct
      if("D_xcy"%in%names(Report)) Array_xct = Report$Epsilon1_sct
      if("D_gcy"%in%names(Report)) Array_xct = Report$Epsilon1_gct
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==7){
      # Epsilon for positive values ("Eps_Pos")
      if( quiet==FALSE ) message(" # Plotting spatio-temporal effects (Epsilon) in 2nd linear predictor")
      if("D_xt"%in%names(Report)) Array_xct = Report$Epsilon2_st
      if("D_xct"%in%names(Report)) Array_xct = Report$Epsilon2_sct
      if("D_xcy"%in%names(Report)) Array_xct = Report$Epsilon2_sct
      if("D_gcy"%in%names(Report)) Array_xct = Report$Epsilon2_gct
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==8){
      # Linear predictor for probability of encounter
      if( quiet==FALSE ) message(" # Plotting 1st predictor after action of link function")
      if("D_xt"%in%names(Report)) Array_xct = Report$P1_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$P1_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$P1_xcy
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==9){
      # Linear predictor for positive catch rates
      if( quiet==FALSE ) message(" # Plotting 2nd predictor after action of link function")
      if("D_xt"%in%names(Report)) Array_xct = Report$P2_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$P2_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$P2_xcy
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==10){
      # Density ("Dens") CV             # Index_xtl
      if( quiet==FALSE ) message(" # Plotting density maps")
      if( is.null(Sdreport) ) stop("Must supply 'Sdreport' if 'plot_num=10'")
      if("D_xt"%in%names(Report)){
        if( !("log(Index_xtl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'SpatialDeltaGLMM'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xtl)"),], dim=c(dim(Report$D_xt),ncol(Report$Index_tl),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,1,'Std. Error']
      }
      if("D_xct"%in%names(Report)){
        if( !("log(Index_xctl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xctl)"),], dim=c(dim(Report$D_xct),dim(Report$Index_ctl)[3],2), dimnames=list(NULL,NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,,1,'Std. Error']
      }
      if("D_xcy"%in%names(Report)){
        if( !("log(Index_xcyl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xcyl)"),], dim=c(dim(Report$D_xcy),dim(Report$Index_cyl)[3],2), dimnames=list(NULL,NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,,1,'Std. Error']
      }
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report))) stop("'plot_num=10' not implemented for 'SpatialVAM'")
      # Convert to CV
      Array_xct = sqrt( exp(Array_xct^2) - 1 )
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
    }
    if(plot_num==11){
      if( quiet==FALSE ) message(" # Plotting covariates")
      if(is.null(TmbData)) stop( "Must provide `TmbData` to plot covariates" )
      #if(!("X_xtp" %in% names(TmbData))) stop( "Can only plot covariates for VAST version >= 2.0.0" )
      if("X_xtp"%in%names(TmbData)) Array_xct = aperm( TmbData$X_xtp, perm=c(1,3,2) )
      if("X_gtp"%in%names(TmbData)) Array_xct = aperm( TmbData$X_gtp, perm=c(1,3,2) )
      category_names = 1:dim(Array_xct)[2]
    }
    if(plot_num==12){
      # Total density ("Dens")
      if( quiet==FALSE ) message(" # Plotting total density")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(apply(Report$D_xct, FUN=sum, MARGIN=c(1,3)))
      if("D_xcy"%in%names(Report)) Array_xct = log(apply(Report$D_xcy, FUN=sum, MARGIN=c(1,3)))
      if("D_gcy"%in%names(Report)) Array_xct = log(apply(Report$D_gcy, FUN=sum, MARGIN=c(1,3)))
      logsum = function(vec){ max(vec) + log(sum(exp(vec-max(vec)))) }
      if("dhat_ktp" %in% names(Report)) Array_xct = apply(aperm(Report$dhat_ktp,c(1,3,2)), FUN=logsum, MARGIN=c(1,3))
      if("dpred_ktp" %in% names(Report)) Array_xct = apply(aperm(Report$dpred_ktp,c(1,3,2)), FUN=logsum, MARGIN=c(1,3))
    }
    if(plot_num==13){
      # Covariate effects for probability of encounter
      if( quiet==FALSE ) message(" # Plotting covariate effects for 1st linear predictor")
      if("D_xt"%in%names(Report)) stop()
      if("D_xct"%in%names(Report)) stop()
      if("D_xcy"%in%names(Report)) Array_xct = Report$eta1_xct
      if("D_gcy"%in%names(Report)) Array_xct = Report$eta1_gct
      if("dhat_ktp" %in% names(Report)) stop()
      if("dpred_ktp" %in% names(Report)) stop()
    }
    if(plot_num==14){
      # Covariate effects for positive catch rates
      if( quiet==FALSE ) message(" # Plotting covariate effects for 2nd linear predictor")
      if("D_xt"%in%names(Report)) stop()
      if("D_xct"%in%names(Report)) stop()
      if("D_xcy"%in%names(Report)) Array_xct = Report$eta2_xct
      if("D_gcy"%in%names(Report)) Array_xct = Report$eta2_gct
      if("dhat_ktp" %in% names(Report)) stop()
      if("dpred_ktp" %in% names(Report)) stop()
    }
    #if(plot_num==15){
    #  # Spatial effects for probability of encounter
    #  if( quiet==FALSE ) message(" # Plotting spatial effects (Omega) for 1st linear predictor")
    #  if("D_xt"%in%names(Report)) stop()
    #  if("D_xct"%in%names(Report)) stop()
    #  if("D_xcy"%in%names(Report)) Array_xct = Report$Omega1_sc %o% 1
    #  if("D_gcy"%in%names(Report)) Array_xct = Report$Omega1_gc %o% 1
    #  if("dhat_ktp" %in% names(Report)) stop()
    #  if("dpred_ktp" %in% names(Report)) stop()
    #}
    #if(plot_num==16){
    #  # Spatial effects for positive catch rates
    #  if( quiet==FALSE ) message(" # Plotting spatial effects (Omega) for 2nd linear predictor")
    #  if("D_xt"%in%names(Report)) stop()
    #  if("D_xct"%in%names(Report)) stop()
    #  if("D_xcy"%in%names(Report)) Array_xct = Report$Omega2_sc %o% 1
    #  if("D_gcy"%in%names(Report)) Array_xct = Report$Omega2_gc %o% 1
    #  if("dhat_ktp" %in% names(Report)) stop()
    #  if("dpred_ktp" %in% names(Report)) stop()
    #}
    if( is.null(Array_xct)) stop("Problem with `plot_num` in `plot_maps(.)")

    # Plot for each category
    if( tolower(Panel)=="category" ){
      if(length(dim(Array_xct))==2) Nplot = 1
      if(length(dim(Array_xct))==3) Nplot = dim(Array_xct)[2]
      for( cI in 1:Nplot){
        if(length(dim(Array_xct))==2) Return = Mat_xt = Array_xct
        if(length(dim(Array_xct))==3) Return = Mat_xt = array(as.vector(Array_xct[,cI,]),dim=dim(Array_xct)[c(1,3)])

        # Do plot
        #if( is.null(mfrow)) mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
        #if(add==FALSE) par( mfrow=mfrow )
        file_name = paste0(plot_code, ifelse(Nplot>1, paste0("--",category_names[cI]), "") )
        plot_args = plot_variable( Y_gt=Mat_xt[,Years2Include,drop=FALSE], map_list=list("PlotDF"=PlotDF, "MapSizeRatio"=MapSizeRatio), projargs=projargs, working_dir=working_dir,
          panel_labels=Year_Set[Years2Include], file_name=file_name, n_cells=n_cells, ... )
      }
    }
    # Plot for each year
    if( tolower(Panel)=="year" ){
      Nplot = length(Years2Include)
      for( tI in 1:Nplot){
        if(length(dim(Array_xct))==2) Mat_xc = Array_xct[,Years2Include[tI],drop=TRUE]
        if(length(dim(Array_xct))==3) Mat_xc = Array_xct[,,Years2Include[tI],drop=TRUE]
        Return = Mat_xc = array( as.vector(Mat_xc), dim=c(dim(Array_xct)[1],Ncategories)) # Reformat to make sure it has same format for everything

        # Do plot
        #if( is.null(mfrow)) mfrow = c(ceiling(sqrt(length(category_names))), ceiling(length(category_names)/ceiling(sqrt(length(category_names)))))
        #if(add==FALSE) par( mfrow=mfrow )
        file_name = paste0(plot_code, ifelse(Nplot>1, paste0("--",Year_Set[Years2Include][tI]), "") )
        plot_args = plot_variable( Y_gt=Mat_xc, map_list=list("PlotDF"=PlotDF, "MapSizeRatio"=MapSizeRatio), projargs=projargs, working_dir=working_dir,
          panel_labels=category_names, file_name=file_name, n_cells=n_cells, ... )
      }
    }
  }
  if( is.null(Return) & quiet==FALSE ) message(" # No plots selected in `plot_set`")

  return( invisible(Return) )
}
