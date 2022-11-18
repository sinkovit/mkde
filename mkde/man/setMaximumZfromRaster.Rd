\name{setMaximumZfromRaster}
\alias{setMaximumZfromRaster}
\title{
Initialize maximum z-axis value from a raster
}
\description{
Set the upper bounds in the z-dimension for each location in the x and y
dimensions from a raster.
}
\usage{
setMaximumZfromRaster(mkde.obj, raster.layer)	
}
\arguments{
  \item{mkde.obj}{2D or 3D MKDE object created with
    \code{initialize3DMKDE} or \code{initialize3DMKDE}, respectively}
  \item{raster.layer}{A RasterLayer object representing the lower bounds
    of the space the animal may occupy in the z-dimension.}
}
\details{
This function sets the upper bounds of the space the animal may occupy
in the z-dimension. For example, the \code{ascii.raster.filez} argument
  may represent a raster for elevation for subterranean animals, or
  other surface.
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
library(raster)
data(dugongdem)
cell.sz <- mean(res(dugongdem))
ext <- extent(dugongdem)
nx <- ncol(dugongdem)
ny <- nrow(dugongdem)
mkde.obj <- initializeMKDE3D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny, min(values(dugongdem), na.rm=TRUE), 50.0, 15)
mkde.obj <- setMaximumZfromRaster(mkde.obj, dugongdem)
}
