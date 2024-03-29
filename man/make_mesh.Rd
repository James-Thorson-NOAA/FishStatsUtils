% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_mesh.R
\name{make_mesh}
\alias{make_mesh}
\title{Make mesh for distances among points}
\usage{
make_mesh(
  loc_x,
  loc_g,
  loc_i,
  Method,
  Extrapolation_List,
  anisotropic_mesh = NULL,
  fine_scale = FALSE,
  map_data,
  mesh_package = c("INLA", "fmesher"),
  ...
)
}
\arguments{
\item{loc_x}{location (eastings and northings in kilometers, UTM) for each sample or knot}

\item{Method}{spatial method determines ("Mesh" and "Grid" give}

\item{anisotropic_mesh}{OPTIONAL, anisotropic mesh (if missing, its recalculated from loc_x)}

\item{...}{Arguments passed to \code{INLA::inla.mesh.create}}
}
\value{
Tagged list containing distance metrics
}
\description{
\code{make_mesh} builds a tagged list representing distances for isotropic or geometric anisotropic triangulated mesh
}
