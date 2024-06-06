\name{mkde3Dinteraction}
\alias{mkde3Dinteraction}
\title{
Probability of 3D spatial-temporal interaction.
}
\description{
  Metric of 3D spatial-temporal interaction.
}
\usage{
mkde3Dinteraction(mkde.obj, move.dat0, move.dat1, t.step, d.thresh)	
}
\arguments{
  \item{mkde.obj}{An MKDE object created with \code{initialize3DMKDE}}
  \item{move.dat0}{A move data object for the first individual created
    with \code{initializeMovementData}}
  \item{move.dat1}{A move data object for the second individual created
    with \code{initializeMovementData}}
  \item{t.step}{A time step used for numerical integration
    over the movement trajectory}
  \item{d.thresh}{The value of the kernel below which its contibrution
    to the overal density is considered negligible}
}
\details{
  This function assumes that the two individual animals were observed at
  the same times. The voxel values returned in the mkde.obj can be summed 
  to obtain a global measure of spatio-temporal interaction.
}
\value{
  Returns a list with the following elements:
  \item{mkde.obj}{An updated MKDE object containing the voxel-level
    Bhattacharyya coefficients}
  \item{move.dat0}{A move data object for the first individuals with
    updated variance parameters}
  \item{move.dat1}{A move data object for the second individuals with
    updated variance parameters}
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov}  \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org} \cr
  Jun Zhu, PhD \cr
  University of Wisconsin-Madison \cr
  \email{jzhu@stat.wisc.edu}
}
\examples{
library(terra)
data(condor)
condor <- condor[1:4,] # simply to make example run more quickly
mv.dat0 <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)
n <- nrow(condor)
v <- 20.0 # increase from 0 to increase difference between move trajectories
vz <- 5.0
mv.dat1 <- initializeMovementData(condor$time, condor$x+rnorm(n, 0, v), 
condor$y+rnorm(n, 0, v), z.obs=condor$z+rnorm(n, 0, vz), sig2obs=25.0, 
sig2obs.z=81.0, t.max=65.0)

fpath <- system.file("extdata", "condordem120.RDS", package="mkde")
condordem120 <- terra::readRDS(fpath)
# next two lines reduce extent of 2D space to speed execution of example
tmp <- ext(c(range(condor$x) + c(-100, 100), range(condor$y) + c(-100, 100)))
condordem120 <- crop(condordem120, tmp)
cell.sz <- mean(res(condordem120))
ext <- ext(condordem120)
nx <- ncol(condordem120)
ny <- nrow(condordem120)
nz <- ceiling(3000.0/cell.sz)
mkde.obj <- initializeMKDE3D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz,
ny, 0.0, cell.sz, nz)

res <- mkde3Dinteraction(mkde.obj, mv.dat0, mv.dat1, 10.0, 1e-20)
mkde.obj <- res$mkde.obj
mv.dat0 <- res$move.dat0
mv.dat1 <- res$move.dat1
}
