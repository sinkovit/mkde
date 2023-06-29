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

# Find min/max coordinates and add buffer
xmax = max(condor$x) + 1000
xmin = min(condor$x) - 1000
ymax = max(condor$y) + 1000
ymin = min(condor$y) - 1000

# Calculate grid dimensions
xrange <- xmax - xmin
yrange <- ymax - ymin
cell.sz = 150
nx <- as.integer(xrange/cell.sz)
ny <- as.integer(yrange/cell.sz)

mv.dat <- initializeMovementData(condor$t, condor$x, condor$y, sig2obs=25.0, t.max=185.0)
mkde.obj <- initializeMKDE2D(xmin, cell.sz, nx, ymin, cell.sz, ny)
dens.res <- initializeDensity(mkde.obj, mv.dat)

mkde.obj <- dens.res$mkde.obj
mv.dat <- dens.res$move.dat
my.quantiles <- c(0.95, 0.75, 0.50)
res <- computeContourValues(mkde.obj, my.quantiles)
print(res)
}
