\name{computeSizeMKDE}
\alias{computeSizeMKDE}
\title{
Compute the area or volume of an MKDE.
}
\description{
For 2D MKDEs, this function computes the area.  For a 2.5D MKDE it
computes the area based on the cell areas computed from an elevation
raster.  For a 3D MKDE, it computes the volume.
}
\usage{
computeSizeMKDE(mkde.obj, prob)	
}
\arguments{
  \item{mkde.obj}{MKDE list object with density estimate calculated.}
  \item{prob}{Probabilities for the desired contours.}
}
\details{
For a 2D or 2.5D MKDE list object, areas within the contours specified by quant are calculated.  For a 3D MKDE, the volumes within the contours are calculated.
}
\value{
A vector with the areas or volumes is returned.
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov} \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org}
}
\examples{
library(raster)
data(condor)
condor <- condor[1:20,] # simply to make example run more quickly
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

data(condordem120)
cell.sz <- mean(res(condordem120))
ext <- extent(condordem120)
nx <- ncol(condordem120)
ny <- nrow(condordem120)
mkde.obj <- initializeMKDE3D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny, min(values(condordem120), na.rm=TRUE), cell.sz, 25)

# note: we use a raster coarse integration time step so the
# example runs faster
dens.res <- initializeDensity(mkde.obj, mv.dat, integration.step=10.0)
mkde.obj <- dens.res$mkde.obj
mv.dat <- dens.res$move.dat

my.quantiles <- c(0.95, 0.75, 0.50)
res <- computeContourValues(mkde.obj, my.quantiles)
res$volume <- computeSizeMKDE(mkde.obj, my.quantiles)
print(res)
}
