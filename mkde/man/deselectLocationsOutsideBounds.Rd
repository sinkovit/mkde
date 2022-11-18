\name{deselectLocationsOutsideBounds}
\alias{deselectLocationsOutsideBounds}
\title{
Flag 3D locations with out-of-bounds z-coordinates
}
\description{
For a 3D MKDE object and 3D location data, the z-coordinates in the
location data are checked to make sure they are within the lower and
upper bounds specified in the MKDE list object.
}
\usage{
deselectLocationsOutsideBounds(move.dat, mkde.obj) 	
}
\arguments{
  \item{move.dat}{A move data object created with \code{initializeMovementData}}
  \item{mkde.obj}{An MKDE object created with \code{initialize2DMKDE} or
    \code{initialize3DMKDE}}
}
\details{
If a 2D or 2.5D MDKE object is passed as an argument, no change is made
to the movement data list object.  If a 3D MKDE list object is passed as
an argument, the z-coordinates in the movement data are checked to
determine if they are in range.  If they are not, the corresponding
value in move.dat$use.obs is set to FALSE.  Note that this function is
not called within \code{initialzeDensity}.  If you want to exclude locations
because they are outside of the allowed range in the z-dimension, this
function must be used before computing the density.
}
\value{
An updated move data list object is returned.
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
data(condor)
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

data(condordem)
cell.sz <- mean(res(condordem))
ext <- extent(condordem)
nx <- ncol(condordem)
ny <- nrow(condordem)
mkde.obj <- initializeMKDE3D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny, min(values(condordem), na.rm=TRUE), 30.0, 100)
mkde.obj <- setMinimumZfromRaster(mkde.obj, condordem)
mv.dat <- deselectLocationsOutsideBounds(mv.dat, mkde.obj)
}
