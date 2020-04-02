
#' Build objects related to spatial information
#'
#' \code{make_spatial_info} builds a tagged list with all the spatial information needed for \code{Data_Fn}
#'
#' \code{fine_scale=TRUE} is a new feature starting in V8.0.0 which triggers two major changes: (1) projecting Gaussian Markov random fields from knots to sampling and extrapolation-grid locations using bilinear interpolation (i.e., piecewise linear smoothing), and (2) including density covariates individually for extrapolation-grid and sampling locations. \code{fine_scale=FALSE} is designed to be backwards compatible with earlier versions, although V8.0.0 may also require changes to input naming conventions for covariates to specify the same model and attain the same fit.
#'
#' \code{LON_intensity} and \code{LAT_intensity} allow users to specify locations that are used by the k-means algorithm to determine the location of knots, where e.g. users can either hard-code the desired knot locations via these inputs (using \code{n_x} greater than this number of locations), or use the extrapolation-grid to ensure that knots are located proportional to that grid.
#'
#' @param n_x the number of vertices in the SPDE mesh (determines the spatial resolution when Method="Mesh")
#' @param Lon_i Longitude for each sample
#' @param Lat_i Latitude for each sample
#' @param knot_method whether to determine location of GMRF vertices based on the location of samples \code{knot_method=`samples`} or extrapolation-grid cells within the specified strata \code{knot_method='grid'}; default \code{knot_method=NULL} is coerced to \code{knot_method=`samples`}
#' @param Method a character of either "Grid" or "Mesh" where "Grid" is a 2D AR1 process, and "Mesh" is the SPDE method with geometric anisotropy
#' @param fine_scale a Boolean indicating whether to ignore (\code{fine_scale=FALSE}) or account for (\code{fine_scale=TRUE}) fine-scale spatial heterogeneity;  See details for more informatino
#' @param Extrapolation_List the output from \code{Prepare_Extrapolation_Data_Fn}
#' @param grid_size_km the distance between grid cells for the 2D AR1 grid (determines spatial resolution when Method="Grid") when not using \code{Method="Spherical_mesh"}
#' @param grid_size_LL the distance between grid cells for the 2D AR1 grid (determines spatial resolution when Method="Grid") when using \code{Method="Spherical_mesh"}
#' @param Network_sz_LL data frame with "parent_s", "child_s", "dist_s", "Lat", "Lon", default=NULL only needed with Method == "Stream_network"
#' @param ... additional arguments passed to \code{INLA::inla.mesh.create}
#' @inheritParams Calc_Kmeans

#' @return Tagged list containing objects for running a VAST model
#' \describe{
#'   \item{MeshList}{A tagged list with inputs related to the SPDE mesh}
#'   \item{GridList}{A tagged list with inputs related to the 2D AR1 grid}
#'   \item{a_xl}{A data frame with areas for each knot and each strattum}
#'   \item{loc_UTM}{A data frame with the converted UTM coordinates for each sample}
#'   \item{Kmeans}{Output from \code{Calc_Kmeans} with knots for a triangulated mesh}
#'   \item{knot_i}{The knot associated with each sample}
#'   \item{Method}{The Method input (for archival purposes)}
#'   \item{loc_x}{The UTM location for each knot}
#' }
#'

#' @export
make_spatial_info = function( n_x, Lon_i, Lat_i, Extrapolation_List, knot_method=NULL, Method="Mesh",
  grid_size_km=50, grid_size_LL=1, fine_scale=FALSE, Network_sz_LL=NULL,
  iter.max=1000, randomseed=1, nstart=100, DirPath=paste0(getwd(),"/"), Save_Results=FALSE,
  LON_intensity, LAT_intensity, ... ){

  # Deprecated options
  if( Method=="Spherical_mesh" ){
    stop("Method=`Spherical_mesh` is not being maintained, but please write the package author if you need to explore this option")
  }
  if( Method == "Stream_network" & fine_scale == TRUE ){
    stop("Please use fine_scale=FALSE with stream network spatial model; feature fine_scale=TRUE not yet supported for stream network spatial model.")
  }

  # Backwards compatibility for when settings didn't include knot_method, such that settings$knot_method=NULL
  if(is.null(knot_method)) knot_method = "samples"

  # Backwards compatible option for different extrapolation grid
  if( missing(LON_intensity) & missing(LAT_intensity) ){
    if( knot_method=="samples" ){
      LON_intensity = Lon_i
      LAT_intensity = Lat_i
    }
    if( knot_method=="grid" ){
      which_rows = which( Extrapolation_List$Data_Extrap[,'Include']==TRUE & Extrapolation_List[["Area_km2_x"]]>0 & rowSums(Extrapolation_List[["a_el"]])>0 )
      LON_intensity = Extrapolation_List$Data_Extrap[ which_rows, 'Lon']
      LAT_intensity = Extrapolation_List$Data_Extrap[ which_rows, 'Lat']
    }
    if( !(knot_method %in% c("samples","grid")) ) stop("`knot_method` must be either `samples` or `grid`")
  }

  # Convert to an Eastings-Northings coordinate system
  if( Method=="Spherical_mesh" ){
    loc_i = data.frame( 'Lon'=Lon_i, 'Lat'=Lat_i )
    # Bounds for 2D AR1 grid
    Grid_bounds = (grid_size_km/110) * apply(loc_e/(grid_size_km/110), MARGIN=2, FUN=function(vec){trunc(range(vec))+c(0,1)})

    # Calculate k-means centroids
    Kmeans = Calc_Kmeans(n_x=n_x, loc_orig=loc_i[,c("Lon", "Lat")], randomseed=randomseed, ... )

    # Calculate grid for 2D AR1 process
    loc_grid = expand.grid( 'Lon'=seq(Grid_bounds[1,1],Grid_bounds[2,1],by=grid_size_LL), 'Lat'=seq(Grid_bounds[1,2],Grid_bounds[2,2],by=grid_size_LL) )
    Which = sort(unique(RANN::nn2(data=loc_grid, query=Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('Lon','Lat')], k=1)$nn.idx[,1]))
    loc_grid = loc_grid[Which,]
    grid_num = RANN::nn2( data=loc_grid, query=loc_i, k=1)$nn.idx[,1]
  }
  if( Method %in% c("Mesh","Grid","Stream_network") ){
    loc_i = project_coordinates( X=Lon_i, Y=Lat_i, projargs=Extrapolation_List$projargs )
    loc_intensity = project_coordinates( X=LON_intensity, Y=LAT_intensity, projargs=Extrapolation_List$projargs )
    colnames(loc_i) = colnames(loc_intensity) = c("E_km", "N_km")
    # Bounds for 2D AR1 grid
    Grid_bounds = grid_size_km * apply(Extrapolation_List$Data_Extrap[,c('E_km','N_km')]/grid_size_km, MARGIN=2, FUN=function(vec){trunc(range(vec))+c(0,1)})

    # Calculate k-means centroids
    Kmeans = Calc_Kmeans(n_x=n_x, loc_orig=loc_intensity[,c("E_km", "N_km")], randomseed=randomseed, nstart=nstart, DirPath=DirPath, Save_Results=Save_Results )
    NN_i = RANN::nn2( data=Kmeans[["centers"]], query=loc_i, k=1)$nn.idx[,1]

    # Calculate grid for 2D AR1 process
    loc_grid = expand.grid( 'E_km'=seq(Grid_bounds[1,1],Grid_bounds[2,1],by=grid_size_km), 'N_km'=seq(Grid_bounds[1,2],Grid_bounds[2,2],by=grid_size_km) )
    Which = sort(unique(RANN::nn2(data=loc_grid, query=Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('E_km','N_km')], k=1)$nn.idx[,1]))
    loc_grid = loc_grid[Which,]
    grid_num = RANN::nn2( data=loc_grid, query=loc_i, k=1)$nn.idx[,1]
  }

  # Calc design matrix and areas
  if( Method=="Grid" ){
    knot_i = grid_num
    loc_x = loc_grid
  }
  if( Method %in% c("Mesh","Spherical_mesh") ){
    knot_i = NN_i
    loc_x = Kmeans[["centers"]]
  }
  if( Method == "Stream_network" ){
    knot_i = Extrapolation_List$Data_Extrap[,"child_i"]
    loc_x = project_coordinates( X=Network_sz_LL[,"Lon"], Y=Network_sz_LL[,"Lat"], projargs=Extrapolation_List$projargs )
    colnames(loc_x) = c('E_km', 'N_km')
  }

  # Bookkeeping for extrapolation-grid
  if( fine_scale==FALSE ){
    loc_g = loc_x
  }
  if( fine_scale==TRUE ){
    loc_g = Extrapolation_List$Data_Extrap[ which(Extrapolation_List$Area_km2_x>0), c('E_km','N_km') ]
  }

  # Convert loc_x back to location in lat-long coordinates latlon_x
  origargs = "+proj=longlat +ellps=WGS84"
  latlon_x = project_coordinates( X=loc_x[,"E_km"], Y=loc_x[,"N_km"], projargs=origargs, origargs=Extrapolation_List$projargs )[,c("Y","X")]
  colnames(latlon_x) = c("Lat", "Lon")

  # Convert loc_g back to location in lat-long coordinates latlon_g
  latlon_g = project_coordinates( X=loc_g[,"E_km"], Y=loc_g[,"N_km"], projargs=origargs, origargs=Extrapolation_List$projargs )[,c("Y","X")]
  colnames(latlon_g) = c("Lat", "Lon")

  # Bundle lat-lon
  latlon_i = cbind( 'Lat'=Lat_i, 'Lon'=Lon_i )

  # Make mesh and info for anisotropy  SpatialDeltaGLMM::
  MeshList = Calc_Anisotropic_Mesh( Method=Method, loc_x=Kmeans$centers, loc_g=loc_g, loc_i=loc_i, Extrapolation_List=Extrapolation_List, fine_scale=fine_scale, ... )
  n_s = switch( tolower(Method), "mesh"=MeshList$anisotropic_spde$n.spde, "grid"=nrow(loc_x), "spherical_mesh"=MeshList$isotropic_spde$n.spde, "stream_network"=nrow(loc_x) )

  # Make matrices for 2D AR1 process
  Dist_grid = dist(loc_grid, diag=TRUE, upper=TRUE)
  M0 = as( ifelse(as.matrix(Dist_grid)==0, 1, 0), "dgTMatrix" )
  M1 = as( ifelse(as.matrix(Dist_grid)==grid_size_km, 1, 0), "dgTMatrix" )
  M2 = as( ifelse(as.matrix(Dist_grid)==sqrt(2)*grid_size_km, 1, 0), "dgTMatrix" )
  if( Method=="Spherical_mesh" ) GridList = list("M0"=M0, "M1"=M1, "M2"=M2, "grid_size_km"=grid_size_LL)
  if( Method %in% c("Mesh","Grid","Stream_network") ) GridList = list("M0"=M0, "M1"=M1, "M2"=M2, "grid_size_km"=grid_size_km)

  # Make projection matrices
  if( fine_scale==FALSE ){
    A_is = matrix(0, nrow=nrow(loc_i), ncol=n_s)
    A_is[ cbind(1:nrow(loc_i),knot_i) ] = 1
    A_is = as( A_is, "dgTMatrix" )
    A_gs = as( diag(n_x), "dgTMatrix" )
  }
  if( fine_scale==TRUE ){
    A_is = INLA::inla.spde.make.A( MeshList$anisotropic_mesh, loc=as.matrix(loc_i) )
    if( class(A_is)=="dgCMatrix" ) A_is = as( A_is, "dgTMatrix" )
    A_gs = INLA::inla.spde.make.A( MeshList$anisotropic_mesh, loc=as.matrix(loc_g) )
    if( class(A_gs)=="dgCMatrix" ) A_gs = as( A_gs, "dgTMatrix" )
    Check_i = apply( A_is, MARGIN=1, FUN=function(vec){sum(vec>0)})
    Check_g = apply( A_is, MARGIN=1, FUN=function(vec){sum(vec>0)})
    if( any(c(Check_i,Check_g) <= 0 ) ){
      # stop("Problem with boundary")
      # plot(MeshList$anisotropic_mesh)
      # points( x=loc_i[which(Check_i!=3),'E_km'], y=loc_i[which(Check_i!=3),'N_km'], col="red" )
    }
  }

  # Calculate areas
  if( Method != "Stream_network" ){
    PolygonList = Calc_Polygon_Areas_and_Polygons_Fn( loc_x=loc_x, Data_Extrap=Extrapolation_List[["Data_Extrap"]], a_el=Extrapolation_List[["a_el"]])
    if( fine_scale==FALSE ){
      a_gl = PolygonList[["a_xl"]]
    }
    if( fine_scale==TRUE ){
      a_gl = as.matrix(Extrapolation_List[["a_el"]][ which(Extrapolation_List$Area_km2_x>0), ])
    }
  }else{
    PolygonList = NULL
    dist_inp = Network_sz_LL[,"dist_s"]
    dist_inp[which(is.infinite(dist_inp))] <- 0
    a_gl = matrix(dist_inp, nrow=n_x)
  }

  # Moving
  if( fine_scale==TRUE | Method=="Stream_network" ){
    g_e = rep(NA, length(Extrapolation_List[["Area_km2_x"]]))
    g_e[ which(Extrapolation_List[["Area_km2_x"]]>0) ] = 1:length(which(Extrapolation_List[["Area_km2_x"]]>0))
  }else{
    g_e = PolygonList$NN_Extrap$nn.idx[,1]
    g_e[ which(Extrapolation_List[["Area_km2_x"]]==0) ] = NA
  }

  # Return
  Return = list( "fine_scale"=fine_scale, "A_is"=A_is, "A_gs"=A_gs, "n_x"=n_x, "n_s"=n_s, "n_g"=nrow(a_gl), "n_i"=nrow(loc_i),
    "MeshList"=MeshList, "GridList"=GridList, "a_gl"=a_gl, "a_xl"=a_gl, "Kmeans"=Kmeans, "knot_i"=knot_i,
    "loc_i"=as.matrix(loc_i), "loc_x"=as.matrix(loc_x), "loc_g"=as.matrix(loc_g), "g_e"=g_e,
    "Method"=Method, "PolygonList"=PolygonList, "NN_Extrap"=PolygonList$NN_Extrap, "knot_method"=knot_method,
    "latlon_x"=latlon_x, "latlon_g"=latlon_g, "latlon_i"=latlon_i )
  class(Return) = "make_spatial_info"
  return( Return )
}

#' Plot spatial representation used for model
#'
#' @title Plot spatial information
#' @param x Output from \code{\link{make_spatial_info}}
#' @param ... Not used
#' @return NULL
#' @method plot make_spatial_info
#' @export
plot.make_spatial_info <- function(x, ...)
{
  cat("\n### Running `plot.make_spatial_info`\n")
  if( x$Method == "Mesh" ){
    plot(x$MeshList$anisotropic_mesh)
    ans = x$MeshList$anisotropic_mesh
  }else{
    cat( "`plot` not implemented for `Method` used\n" )
  }

  invisible(ans)
}

#' Print spatial representation used for model
#'
#' @title Print spatial information
#' @param x Output from \code{\link{make_spatial_info}}
#' @param ... Not used
#' @return NULL
#' @method print make_spatial_info
#' @export
print.make_spatial_info <- function(x, ...)
{
  cat("\n### Running `print.make_spatial_info`\n")
  cat( paste0("`Method` = ", x$Method,"\n") )
  cat( paste0("`n_x` = ", x$n_x,"\n") )
  cat( paste0("`n_s` = ", x$n_s,"\n") )
  cat( paste0("Extrapolating results to `n_g` = ", x$n_g, " extrapolation-grid cells\n") )

  return(invisible(NULL))
}
