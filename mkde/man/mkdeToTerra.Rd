\name{mkdeToTerra}
\alias{mkdeToTerra}
\title{
MKDE to RasterLayer or RasterStack
}
\description{
Converts an MKDE into a RasterLayer or RasterStack so that raster
package functions can be used.
}
\usage{
mkdeToTerra(mkde.obj)	
}
\arguments{
  \item{mkde.obj}{A 2D, 2.5D, or 3D MKDE object created with
    \code{initialize2DMKDE} or \code{initialize3DMKDE} and density
    initialized with \code{initializeDensity}}
}
\details{
This method converts the density array in the MKDE oject to an object
of a class from the raster package.  This allows the functions in the
raster package to be used with the MKDEs.
}
\value{
If the MKDE is 2D or 2.5D, a RasterLayer object is returned.  If the
MKDE is 3D, a RasterStack is returned with one layer in the stack per
level.
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

# set up MKDE object
fpath <- system.file("extdata", "pandadem.RDS", package="mkde")
pandadem <- terra::readRDS(fpath)
cell.sz <- mean(res(pandadem))
ext <- ext(pandadem)
nx <- ncol(pandadem)
ny <- nrow(pandadem)
mkde.obj <- initializeMKDE2D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz, ny)

# set up movement data
data(panda)
mv.dat <- initializeMovementData(panda$time, panda$x, panda$y, 
sig2obs=25.0, t.max=185.0)

# estimate density
dens.res <- initializeDensity(mkde.obj, mv.dat)
mkde.obj <- dens.res$mkde.obj
mv.dat <- dens.res$move.dat
mkde.rst <- mkdeToTerra(mkde.obj)
plot(mkde.rst)
}
