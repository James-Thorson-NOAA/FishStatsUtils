% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{map_hypervariance}
\alias{map_hypervariance}
\title{Plot variance of GMRF knots}
\usage{
map_hypervariance(report, Spatial_List, method)
}
\arguments{
\item{report}{Report from TmbObj (e.g., TmbList[["Obj"]])}

\item{Spatial_List}{Output from FishStatsUtils::Spatial_Information_Fn()}

\item{method}{Choose whether to plot anisotropic or isotropic covariance matrix}
}
\description{
\code{map_hypervariance} Plot variance of GMRF knots
}
\examples{
map_hypervariance(report = Save$Report,
                  Spatial_List = Spatial_List,
                  method = "anisotropic")
}
