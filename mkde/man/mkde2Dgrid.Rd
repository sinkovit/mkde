\name{mkde2Dgrid}
\alias{mkde2Dgrid}
\title{
Movement-based kernel density estimate (MKDE) in 2D using Rcpp
}
\description{
Provides a function for 2-dimensional MKDEs.
}
\usage{
mkde2Dgrid(mkde.obj, move.dat, t.step, d.thresh)	
}
\arguments{
  \item{mkde.obj}{A 2D or 2.5D MKDE object}
  \item{move.dat}{A move data object}
  \item{t.step}{An integration time step}
  \item{d.thresh}{A kernel density threshold}
}
\details{
  This is  lower-level function that call the C++ function. for
  estimating the movement-based density in 2D.  In practice, users
  should call \code{initializeDensity}. \cr
  The argument d.thresh is a univariate probability density beyond which
  the kernel contribution to the overall MKDE is assumed to be
  negligible. Usually this is set at a very small value and is used to
  prevent calculations from being performed in cells to which the kernel
  makes a negligible contribution.
}
\value{
  An array whose elements are the estimated utilization probabilities
  for each cell.
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov}  \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org} \cr
  Robert Sinkovits, PhD \cr
  San Diego Supercomputer Center \cr
  \email{sinkovit@sdsc.edu } \cr
  Glenn Lockwood, PhD \cr
  San Diego Supercomputer Center \cr
  \email{glock@sdsc.edu} \cr
  Jun Zhu, PhD \cr
  University of Wisconsin-Madison \cr
  \email{jzhu@stat.wisc.edu}
}
\examples{
library(raster)
data(panda)

# Find min/max coordinates and add buffer
xmax = max(panda$x) + 100
xmin = min(panda$x) - 100
ymax = max(panda$y) + 100
ymin = min(panda$y) - 100

# Calculate grid dimensions
xrange <- xmax - xmin
yrange <- ymax - ymin
cell.sz = 30
nx <- as.integer(xrange/cell.sz)
ny <- as.integer(yrange/cell.sz)

mv.dat <- initializeMovementData(panda$time, panda$x, panda$y, t.max=185.0, sig2obs=25.0)
if (all(is.na(mv.dat$sig2xy))) {
  mv.dat <- estVarMKDE(mv.dat)
}

mkde.obj <- initializeMKDE2D(xmin, cell.sz, nx, ymin, cell.sz, ny)
res <- mkde2Dgrid(mkde.obj, mv.dat, 10.0, 1e-20)
}
