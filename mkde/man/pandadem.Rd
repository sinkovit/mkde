\name{pandadem}
\alias{pandadem}
\docType{data}
\title{
A RasterLayer containing a digital elevation model (DEM) for the panda data
}
\description{
Digital elevation model of the terrain within the panda's range with a cell resolution of 30 meters (ASTER GDEM version 2, Tachikawa et. al. 2011).
}
\usage{data(pandadem)}
\format{
  A raster package RasterLayer object with 58 rows, 78 columns, and
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
data(pandadem)
plot(pandadem)
}
\keyword{datasets}
