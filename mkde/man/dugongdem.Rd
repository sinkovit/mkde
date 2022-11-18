\name{dugongdem}
\alias{dugongdem}
\docType{data}
\title{
A RasterLayer containing bathymetry for the dugong data
}
\description{
A 10-meter resolution bathymetric surface of the core seagrass habitat in depth below the mean sea level within the dugong's range (Sheppard et. al. 2007).
}
\usage{data(dugongdem)}
\format{
  A raster package RasterLayer object with 220 rows, 485 columns, and
  a cell size of 10 by 10 meters. Each cell contains a depth at low tide
  in meters.
}
%\details{
% anything?
%}
\source{
James Sheppard, San Diego Zoo Institute for Conservation Research
}
\references{
Sheppard, J. K., Lawler, I. R., & Marsh, H. (2007). Seagrass as pasture for seacows: Landscape-level dugong habitat evaluation. Estuarine, Coastal and Shelf Science, 71(1), 117-132.
}
\examples{
library(raster)
data(dugongdem)
plot(dugongdem)
}
\keyword{datasets}
