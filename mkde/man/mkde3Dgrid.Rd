\name{mkde3Dgrid}
\alias{mkde3Dgrid}
\title{
Movement-based kernel density estimate (MKDE) in 3D using Rcpp
}
\description{
Provides a function for 3-dimensional MKDEs.
}
\usage{
mkde3Dgrid(mkde.obj, move.dat, t.step, d.thresh)	
}
\arguments{
  \item{mkde.obj}{A 3D MKDE object}
  \item{move.dat}{A move data object}
  \item{t.step}{An integration time step}
  \item{d.thresh}{A kernel density threshold}
}
\details{
  This is  lower-level function that call the C++ function. for
  estimating the movement-based density in 3D.  In practice, users
  should call \code{initializeDensity}. \cr
  The argument \code{d.thresh} is a univariate probability density beyond which
  the kernel contribution to the overall MKDE is assumed to be
  negligible. Usually this is set at a very small value and is used to
  prevent calculations from being performed in cells to which the kernel
  makes a negligible contribution.
}
\value{
  An array whose elements are the estimated utilization probabilities
  for each voxel.
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
data(condor)
condor <- condor[1:10,] # simply to make example run more quickly
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

data(condordem120)
cell.sz <- mean(res(condordem120))
ext <- extent(condordem120)
nx <- ncol(condordem120)
ny <- nrow(condordem120)
nz <- ceiling(3000.0/cell.sz)
mkde.obj <- initializeMKDE3D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny, 0.0, cell.sz, nz)

res <- mkde3Dgrid(mkde.obj, mv.dat, 5.0, 1e-20)
}
