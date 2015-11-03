% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/misc.R
\name{coordinatize}
\alias{coordinatize}
\title{Convert georeferenced data.frames into projected SpatialPointsDataFrames}
\usage{
coordinatize(dt, latname = "latdec", lonname = "londec")
}
\arguments{
\item{latname}{character column name of the "y" coordinate}

\item{lonname}{character column name of the "x" coordinate}
}
\description{
Convert georeferenced data.frames into projected SpatialPointsDataFrames
}
\examples{
dt<-streamget(201002,qa=TRUE)
}
