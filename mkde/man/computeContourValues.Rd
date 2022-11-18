\name{computeContourValues}
\alias{computeContourValues}
\title{
Find thresholds for contour intervals
}
\description{
Find the cell or voxel probabilities that correspond to user-specified probability contours
}
\usage{
computeContourValues(mkde.obj, prob)	
}
\arguments{
  \item{mkde.obj}{An MKDE object with density initialized}
  \item{prob}{Probabilities (i.e. proportions) for desired contours of
    the MKDE}
}
\details{
  This function computes threshold cell or voxel probability values
  corresponding to contours for specified proportions of the utilization
  distribution. Note that the arugment \code{prob} specifies the
  cumulative probability of the cells or voxels within the contour
  corresponding to the cell or voxel threshold probability. The cell or
  voxel threshold probabilities may be orders of magnitude smaller than
  the cumulative probabilities provided in the \code{prob} argument.
}
\value{
  A data frame with the probabilities given in the prob argument and
  corresponding thresholds in the MKDE
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
condor <- condor[1:10,] # simply to make example run more quickly
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
print(res)
}
