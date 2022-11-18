\name{writeInterpolatedPathVTK}
\alias{writeInterpolatedPathVTK}
\title{
Write path to VTK.
}
\description{
Write the interpolated move path to a VTK file.
}
\usage{
writeInterpolatedPathVTK(move.dat, mkde.obj, description, filename, control)	
}
\arguments{
  \item{move.dat}{A move data object created with \code{initializeMovementData}}
  \item{mkde.obj}{An MKDE object created with \code{initialize3DMKDE}}
  \item{description}{A text description for the file header}
  \item{filename}{A string for the path and file name}
  \item{control}{A list for finer control}
}
\details{
Writes 3D lines between observed locations for move steps. Move steps
where the initial location i has a value of FALSE for
\code{move.dat$use.obs[i]} are omitted. Currently, the list argument control
only has three elements: (1) ``method'' with a default value of
``linear'', (2) method.par which is a list of method parameters, and (3)
z.scale used to scale the z-coordinates. Only the \code{z.scale} should be set
by the user at this time.
}
\value{
No value is returned
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

writeInterpolatedPathVTK(mv.dat, mkde.obj,
"Example California condor move steps", "condor_path.vtk")
}
