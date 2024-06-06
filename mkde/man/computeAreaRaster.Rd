\name{computeAreaRaster}
\alias{computeAreaRaster}
\title{
Calculate cell areas from elevation raster.
}
\description{
A lower-level function for calculating a cell area matrix (2D array) from an
elevation matrix.  The area is based on the surface area of the terrain
in the cell.  The area matrix is then used in calculating areas of
2.5D MKDES.
}
\usage{
computeAreaRaster(RelevMatrix, RcellSize)	
}
\arguments{
  \item{RelevMatrix}{A 2D array with elevation values.}
  \item{RcellSize}{Size of the cells.  It is assumed to be the same in
    the x and y dimensions.}
}
\details{
  This is a wrapper function for C++ function that calculates the surface area of each raster cell given the cell elevations.  It is not intended to be used directly; instead, the user should call \code{initializeAreaRaster} on an MKDE object.
}
\value{
  A 2D matrix of cell surface areas.
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov}  \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org}
}
\references{
  Jenness J.S. (2004) Calculating landscape surface area from digital
  elevation models. Wildlife Society Bulletin 32: 829-839. \cr
  Jenness, J.S. (2014) Calculating landscape surface area from
  unprojected digital elevation models. In preparation.
}
\examples{
library(terra)
fpath <- system.file("extdata", "condordem.RDS", package="mkde")
condordem <- terra::readRDS(fpath)
cell.sz <- mean(res(condordem))
area.rast <- computeAreaRaster(as.matrix(condordem, wide=TRUE), cell.sz)
}
