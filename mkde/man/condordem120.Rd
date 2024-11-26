\name{condordem120}
\alias{condordem120}
\docType{data}
\title{
A RasterLayer containing a digital elevation model (DEM) for the condor data
}
\description{
Digital elevation model of the terrain within the condor's range with a cell resolution of 120 meters (ASTER GDEM version 2, Tachikawa et. al. 2011).
}
\format{
  A raster package RasterLayer object with 188 rows, 156 columns, and
  a cell size of 120 by 120 meters. Each cell contains an elevation in meters.
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
URL: http://asterweb.jpl.nasa.gov/gdem.asp
}
\references{
Tachikawa, T., Hato, M., Kaku, M. & Iwasaki, A. (2011). The characteristics of ASTER GDEM version 2, IGARSS.
}
\examples{
library(terra)
fpath <- system.file("extdata", "condordem120.RDS", package="mkde")
condordem120 <- terra::readRDS(fpath)
plot(condordem120)
}
\keyword{datasets}
