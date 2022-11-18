\name{condordem}
\alias{condordem}
\docType{data}
\title{
A RasterLayer containing a digital elevation model (DEM) for the condor data
}
\description{
Digital elevation model of the terrain within the condor's range with a cell resolution of 30 meters (ASTER GDEM version 2, Tachikawa et. al. 2011).
}
\usage{data(condordem)}
\format{
  A raster package RasterLayer object with 749 rows, 623 columns, and
  a cell size of 30 by 30 meters. Each cell contains an elevation in meters.
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
library(raster)
data(condordem)
plot(condordem)
}
\keyword{datasets}
