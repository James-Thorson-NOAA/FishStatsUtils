
#' @title
#' Plot index of abundance
#'
#' @description
#' \code{plot_biomass_index} plots an index proportion to population abundance
#'
#' @inheritParams plot_maps
#' @param TmbData Formatted data inputs, from `VAST::Data_Fn(...)`
#' @param DirName Directory for saving plot and table
#' @param PlotName Name for plot
#' @param interval_width width for confidence intervals
#' @param strata_names names for spatial strata
#' @param category_names names for categories (if using package `VAST`)
#' @param use_biascorr Boolean, whether to use bias-corrected estimates if available
#' @param plot_legend Add legend for labelling colors
#' @param total_area_km2 Total area for calculating a design-based estimator using one design-stratum (only recommended for model exploration)
#' @param plot_log Boolean, whether to plot y-axis in log-scale
#' @param width plot width in inches
#' @param height plot height in inches
#' @param ... Other inputs to `par()`
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{Table}{table of index estimates by stratum and year, e.g., for including in an assessment model}
#' }
#'

#' @export
plot_biomass_index <-
function( TmbData, Sdreport, Year_Set=NULL, Years2Include=NULL, DirName=paste0(getwd(),"/"), PlotName="Index", interval_width=1,
  strata_names=NULL, category_names=NULL, use_biascorr=TRUE, plot_legend=TRUE, total_area_km2=NULL, plot_log=FALSE,
  width=NULL, height=NULL, create_covariance_table=FALSE, ... ){

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")
  if(!is.null(category_names) && length(category_names)!=TmbData$n_c ) stop("`category_names` must have same length as `TmbData$n_c`")
  if(!is.null(Year_Set) && length(Year_Set)!=TmbData$n_t ) stop("`Year_Set` must have same length as `TmbData$n_t`")
  if(!is.null(strata_names) && length(strata_names)!=TmbData$n_l ) stop("`strata_names` must have same length as `TmbData$n_l`")

  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    ParName = "Index_tl"
    TmbData[['n_c']] = 1
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    ParName = "Index_ctl"
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    ParName = "Index_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }
  if( "Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialVAM
    ParName = "Index_tp"
    TmbData[["n_l"]] = 1
    TmbData[["n_c"]] = TmbData[["n_p"]]
  }

  # Add t_i if missing (e.g., from VAST V2.8.0 through V9.3.0)
  if( !("t_i" %in% names(TmbData)) ){
    TmbData$t_i = TmbData$t_iz[,1]
  }

  # Add in t_yz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_yz" %in% names(TmbData)) ){
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol=1)
  }

  # Fill in missing
  mfrow = c( ceiling(sqrt(TmbData$n_c)), ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c))) )
  if( is.null(width)) width = mfrow[2] * 3
  if( is.null(height)) height = mfrow[1] * 3
  if( is.null(Year_Set) ) Year_Set = 1:TmbData$n_t
  if( is.null(Years2Include) ) Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c

  # Logical check
  if( "unbiased"%in%names(Sdreport) ){
    if( all(is.na(Sdreport$unbiased$value)) ){
      stop("You appear to be using bias-correction, but all values are NA. Please report problem to package author.")
    }
  }

  # Defaults
  if( "treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options) ){
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }else{
    treat_missing_as_zero = FALSE
  }

  # Objects
  SD = TMB::summary.sdreport(Sdreport)
  if( !"report" %in% names(as.list(args(TMB:::as.list.sdreport))) ){
    warning( "package `TMB` should be updated to easily access standard errors")
  }
  SD_stderr = TMB:::as.list.sdreport( Sdreport, what="Std. Error", report=TRUE )
  SD_estimate = TMB:::as.list.sdreport( Sdreport, what="Estimate", report=TRUE )
  if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
    SD_estimate_biascorrect = TMB:::as.list.sdreport( Sdreport, what="Est. (bias.correct)", report=TRUE )
  }
  if( any(is.na(SD_estimate)) | any(is.na(SD_stderr)) ){
    stop( "Problem: Standard errors contain NAs")
  }

  # Extract index (using bias-correctino if available and requested)
  if( ParName %in% c("Index_tl","Index_ctl","Index_cyl")){
    Index_ctl = log_Index_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    # Index
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Index_ctl[] = SD[which(rownames(SD)==ParName),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(Index_ctl)) ){
      message("Using bias-corrected estimates for abundance index (natural-scale)...")
    }else{
      message("Not using bias-corrected estimates for abundance index (natural-scale)...")
      Index_ctl[] = SD[which(rownames(SD)==ParName),c('Estimate','Std. Error')]
    }
    # Log-Index
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      log_Index_ctl[] = SD[which(rownames(SD)==paste0("ln_",ParName)),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(log_Index_ctl)) ){
      message("Using bias-corrected estimates for abundance index (log-scale)...")
    }else{
      message("Not using bias-corrected estimates for abundance index (log-scale)...")
      log_Index_ctl[] = SD[which(rownames(SD)==paste0("ln_",ParName)),c('Estimate','Std. Error')]
    }
  }
  if( ParName %in% c("Index_tp")){
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==ParName)],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==paste0("ln_",ParName))],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'" )
      }
    }else{
      Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'" )
      }
    }
  }

  # Extract biomass ratio Bratio_cty if available (only available if >= V5.3.0 and using spatial Gompertz model features)
  if( "Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    Bratio_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Bratio_ctl[] = SD[which(rownames(SD)=="Bratio_cyl"),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(Bratio_ctl)) ){
      message("Using bias-corrected estimates for biomass ratio (natural-scale)...")
    }else{
      message("Not using bias-corrected estimates for biomass ratio (natural-scale)...")
      Bratio_ctl[] = SD[which(rownames(SD)=="Bratio_cyl"),c('Estimate','Std. Error')]
    }
  }else{
    Bratio_ctl = NULL
  }
  if( "ln_Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    log_Bratio_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      log_Bratio_ctl[] = SD[which(rownames(SD)=="ln_Bratio_cyl"),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(log_Bratio_ctl)) ){
      message("Using bias-corrected estimates for biomass ratio (log-scale)...")
    }else{
      message("Not using bias-corrected estimates for biomass ratio (log-scale)...")
      log_Bratio_ctl[] = SD[which(rownames(SD)=="ln_Bratio_cyl"),c('Estimate','Std. Error')]
    }
  }else{
    log_Bratio_ctl = NULL
  }

  # Extract Fratio
  if( "Fratio_ct" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    Fratio_ct = array( NA, dim=c(unlist(TmbData[c('n_c','n_t')]),2), dimnames=list(category_names,Year_Set,c('Estimate','Std. Error')) )
    Fratio_ct[] = SD[which(rownames(SD)=="Fratio_ct"),c('Estimate','Std. Error')]
    #Fratio_ct = abind::abind( SD_estimate$Fratio, SD_stderr$Fratio, along=3 )
    #dimnames(Fratio_ct) = list(category_names,Year_Set,c('Estimate','Std. Error'))
  }else{
    Fratio_ct = NULL
  }

  # Calculate design-based
  if( !is.null(total_area_km2) & TmbData$n_c==1 ){
    message( "Calculating naive design-based index -- do not use this, its intended only for comparison purposes" )
    Calc_design = TRUE
    Design_t = tapply( TmbData$b_i/TmbData$a_i, INDEX=TmbData$t_i, FUN=mean ) * total_area_km2 / 1000 # Convert to tonnes
    Design_t = cbind( "Estimate"=Design_t, "Std. Error"=sqrt(tapply(TmbData$b_i/TmbData$a_i,INDEX=TmbData$t_i,FUN=var)/tapply(TmbData$b_i/TmbData$a_i,INDEX=TmbData$t_i,FUN=length))*total_area_km2/1000)
    Design_t = cbind( Design_t, "CV"=Design_t[,'Std. Error'] / Design_t[,'Estimate'] )
  }else{
    Calc_design = FALSE
  }

  # Fix at zeros any years-category combinations with no data
  if( treat_missing_as_zero==TRUE ){
    # Determine year-category pairs with no data
    Num_ct = tapply( TmbData$b_i, INDEX=list(factor(TmbData$c_i,levels=1:TmbData$n_c-1),factor(TmbData$t_i,levels=1:TmbData$n_t-1)), FUN=function(vec){sum(!is.na(vec))} )
    Num_ct = ifelse( is.na(Num_ct), 0, Num_ct )
    # Replace values with 0 (estimate) and NA (standard error)
    Index_ctl[,,,'Estimate'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, 0, Index_ctl[,,,'Estimate'])
    Index_ctl[,,,'Std. Error'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, NA, Index_ctl[,,,'Std. Error'])
    log_Index_ctl[,,,'Estimate'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, -Inf, log_Index_ctl[,,,'Estimate'])
    log_Index_ctl[,,,'Std. Error'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, NA, log_Index_ctl[,,,'Std. Error'])
  }

  # Plot biomass and Bratio
  Plot_suffix = "Biomass"
  if( !is.null(Bratio_ctl) ) Plot_suffix = c( Plot_suffix, "Bratio" )
  for( plotI in 1:length(Plot_suffix) ){
    if( Plot_suffix[plotI]=="Biomass" ){ Array_ctl = Index_ctl; log_Array_ctl = log_Index_ctl }
    if( Plot_suffix[plotI]=="Bratio" ){ Array_ctl = Bratio_ctl; log_Array_ctl = log_Bratio_ctl }
    plot_index( Index_ctl=array(Index_ctl[,,,'Estimate'],dim(Index_ctl)[1:3]),
      sd_Index_ctl=array(log_Index_ctl[,,,'Std. Error'],dim(log_Index_ctl)[1:3]),
      Year_Set=Year_Set, Years2Include=Years2Include, strata_names=strata_names, category_names=category_names,
      DirName=DirName, PlotName=paste0(PlotName,"-",Plot_suffix[plotI],".png"),
      interval_width=interval_width, width=width, height=height, xlab="Year", ylab="Index",
      scale="log", plot_args=list("log"=ifelse(plot_log==TRUE,"y","")), "Yrange"=c(ifelse(plot_log==TRUE,NA,0),NA) )
  }

  # Plot
  if( !is.null(Fratio_ct) ){
    Array_ct = Fratio_ct
      Array_ct = ifelse( Array_ct==0, NA, Array_ct )
    plot_index( Index_ctl=array(Array_ct[,,'Estimate'],dim(Array_ct)[1:2]), sd_Index_ctl=array(Array_ct[,,'Std. Error'],dim(Array_ct)[1:2]),
      Year_Set=Year_Set, Years2Include=Years2Include, strata_names=strata_names, category_names=category_names,
      DirName=DirName, PlotName=paste0(PlotName,"-Fratio.png"), scale="uniform",
      interval_width=interval_width, width=width, height=height, xlab="Year", ylab="Fishing ratio" )
  }

  # Plot stock status
  if( !is.null(Bratio_ctl) & !is.null(Fratio_ct) ){
    Par = list( mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i", oma=c(1,2,0,0), mfrow=mfrow, ... )
    Col = colorRampPalette(colors=c("blue","purple","red"))
    png( file=paste0(DirName,"/",PlotName,"-Status.png"), width=width, height=height, res=200, units="in")
      par( Par )
      Array1_ct = Bratio_ctl[,,1,]
      Array1_ct = ifelse( Array1_ct==0, NA, Array1_ct )
      Array2_ct = Fratio_ct
      Array2_ct = ifelse( Array2_ct==0, NA, Array2_ct )
      for( cI in 1:TmbData$n_c ){
        # Calculate y-axis limits
        Xlim = c(0, max(1, Array1_ct[cI,Years2Include,'Estimate']%o%c(1,1) + Array1_ct[cI,Years2Include,'Std. Error']%o%c(-interval_width,interval_width),na.rm=TRUE) )
        Ylim = c(0, max(2, Array2_ct[cI,Years2Include,'Estimate']%o%c(1,1) + Array2_ct[cI,Years2Include,'Std. Error']%o%c(-interval_width,interval_width),na.rm=TRUE) )
        # Plot stuff
        plot(1, type="n", xlim=Xlim, ylim=Ylim, xlab="", ylab="", main=ifelse(TmbData$n_c>1,category_names[cI],"") )
        points( x=Array1_ct[cI,Years2Include,'Estimate'], y=Array2_ct[cI,Years2Include,'Estimate'], col=Col(length(Year_Set))[Years2Include] )
        for( tI in Years2Include ){
          lines( x=rep(Array1_ct[cI,tI,'Estimate'],2), y=Array2_ct[cI,tI,'Estimate']+Array2_ct[cI,tI,'Std. Error']*c(-interval_width,interval_width), col=Col(length(Year_Set))[tI] )
          lines( x=Array1_ct[cI,tI,'Estimate']+Array1_ct[cI,tI,'Std. Error']*c(-interval_width,interval_width), y=rep(Array2_ct[cI,tI,'Estimate'],2), col=Col(length(Year_Set))[tI] )
        }
        abline( v=0.4, lty="dotted" )
        abline( h=1, lty="dotted" )
      }
      legend( "topright", bty="n", fill=c(Col(length(Year_Set))[Years2Include[1]],Col(length(Year_Set))[rev(Years2Include)[1]]), legend=c(Year_Set[Years2Include[1]],Year_Set[rev(Years2Include)[1]]) )
      mtext( side=1:2, text=c("Biomass relative to unfished","Fishing relative to F_40%"), outer=TRUE, line=c(0,0) )
    dev.off()
  }

  # Write to file
  Table = NULL
  for( cI in 1:TmbData$n_c ){
    Tmp = data.frame( "Year"=Year_Set, "Unit"=1, "Fleet"=rep(strata_names,each=TmbData$n_t), "Estimate_metric_tons"=as.vector(Index_ctl[cI,,,'Estimate']), "SD_log"=as.vector(log_Index_ctl[cI,,,'Std. Error']), "SD_mt"=as.vector(Index_ctl[cI,,,'Std. Error']) )
    if( TmbData$n_c>1 ) Tmp = cbind( "Category"=category_names[cI], Tmp)
    Table = rbind( Table, Tmp )
  }
  if(!is.null(total_area_km2)) Table = cbind(Table, "Naive_design-based_index"=Design_t)
  write.csv( Table, file=paste0(DirName,"/Table_for_SS3.csv"), row.names=FALSE)

  # Return stuff
  Return = list( "Table"=Table, "log_Index_ctl"=log_Index_ctl, "Index_ctl"=Index_ctl )

  # Extract and save covariance
  if( "cov"%in%names(Sdreport) & create_covariance_table==TRUE ){
    DF = expand.grid( "Category"=1:TmbData$n_c, "Year"=1:TmbData$n_t, "Stratum"=1:TmbData$n_l )
    Which = which( names(Sdreport$value)==ParName )
    Cov = Sdreport$cov[Which,Which]
    Corr = cov2cor(Cov) - diag(nrow(Cov))
    rowcolDF = cbind( "RowNum"=row(Corr)[lower.tri(Corr,diag=TRUE)], "ColNum"=col(Corr)[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( DF[rowcolDF[,'ColNum'],], DF[rowcolDF[,'RowNum'],] )
    colnames(Table) = paste0(colnames(Table), rep(c(1,2),each=3))
    Table = cbind( Table, "Correlation"=cov2cor(Cov)[lower.tri(Corr,diag=TRUE)], "Covariance"=Cov[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( Table, "Index1"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'ColNum'],],1))], "Index2"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'RowNum'],],1))] )
    WhichZero = which( (Table[,'Index1']*Table[,'Index2']) == 0 )
    Table[WhichZero,c('Correlation','Covariance')] = 0
    Return = c( Return, "Table_of_estimted_covariance"=Table )
  }
  if( !is.null(Bratio_ctl)) Return = c( Return, list("Bratio_ctl"=Bratio_ctl) )
  if( !is.null(log_Bratio_ctl)) Return = c( Return, list("log_Bratio_ctl"=log_Bratio_ctl) )
  if( !is.null(Fratio_ct)) Return = c( Return, list("Fratio_ct"=Fratio_ct) )

  return( invisible(Return) )
}
