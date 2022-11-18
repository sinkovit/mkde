\name{initializeAreaRaster}
\alias{initializeAreaRaster}
\title{
Initialize an area raster for a 2.5D MKDE.
}
\description{
Initialize the surface area for a 2.5D MKDE. 
}
\usage{
initializeAreaRaster(mkde.obj)	
}
\arguments{
  \item{mkde.obj}{An MKDE object created using \code{initialize2DMKDE}}
}
\details{
  After creating the MKDE object and setting the lower bounds in the
  z-dimension using \code{setMinimumZfromRaster}, this function computes
  the surface area of each raster cell and sets the dimension of the
  MKDE object to 2.5.
}
\value{
  Returns a 2.5 MKDE object with an initialized area raster.
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
data(pandadem)
cell.sz <- mean(res(pandadem))
ext <- extent(pandadem)
nx <- ncol(pandadem)
ny <- nrow(pandadem)
mkde.obj <- initializeMKDE2D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz, ny)

mkde.obj <- setMinimumZfromRaster(mkde.obj, pandadem)
mkde.obj <- initializeAreaRaster(mkde.obj)
}
