\name{setMinimumZfromRaster}
\alias{setMinimumZfromRaster}
\title{
Set minimum z-axis values from a raster
}
\description{
Set the lower bounds in the z-dimension for each location in the x and y
dimensions from a raster.
}
\usage{
setMinimumZfromRaster(mkde.obj, spat.raster)
}
\arguments{
  \item{mkde.obj}{A 2D or 3D MKDE object created with
    \code{initialize3DMKDE} or \code{initialize3DMKDE}, respectively}
  \item{spat.raster}{A RasterLayer object representing the lower bounds
    of the space the animal may occupy in the z-dimension.}
}
\details{
  This function sets the lower bounds of the space the animal may occupy
  in the z-dimension. For example, the \code{ascii.raster.file} argument
  may represent a raster for elevation, depth of the sea floor, or other
  surface.
}
\value{
An updated MKDE list object is returned.
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov}  \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org}
}
\examples{
library(terra)
fpath <- system.file("extdata", "dugongdem.RDS", package="mkde")
dugongdem <- terra::readRDS(fpath)
cell.sz <- mean(res(dugongdem))
ext <- ext(dugongdem)
nx <- ncol(dugongdem)
ny <- nrow(dugongdem)
mkde.obj <- initializeMKDE3D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz,
ny, min(values(dugongdem), na.rm=TRUE), 50.0, 15)
mkde.obj <- setMinimumZfromRaster(mkde.obj, dugongdem-20.0)
}
