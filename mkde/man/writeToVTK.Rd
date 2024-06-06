\name{writeToVTK}
\alias{writeToVTK}
\title{
Write MKDE to an ASCII Visualization Tool Kit (VTK) structured grid file.
}
\description{
Write the MKDE to a VTK file.
}
\usage{
writeToVTK(mkde.obj, fname, description="3D MKDE", cumprob=FALSE)
}
\arguments{
  \item{mkde.obj}{3D MKDE object created with \code{initialize3DMKDE} and
    density initialized with \code{initializeDensity}}
  \item{fname}{The patch and file name for output VTK file}
  \item{description}{A character string with a brief description that
    will be placed in the header of the VTK file}
  \item{cumprob}{Indicate whether to write the voxel probabilities of
    cumulative probabilities.}
}
\details{
  This function writes a VTK structured grid file for a 3D MKDE that can
  be used in 3D visualization tools such as MayaVI or ParaView
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
library(terra)
data(condor)
condor <- condor[1:20,] # simply to make example run more quickly
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

fpath <- system.file("extdata", "condordem120.RDS", package="mkde")
condordem120 <- terra::readRDS(fpath)
# next two lines reduce extent of 2D space to speed execution of example
tmp <- ext(c(range(condor$x) + c(-100, 100), range(condor$y) + c(-100, 100)))
condordem120 <- crop(condordem120, tmp)

cell.sz <- mean(res(condordem120))
ext <- ext(condordem120)
nx <- ncol(condordem120)
ny <- nrow(condordem120)
mkde.obj <- initializeMKDE3D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz,
ny, min(values(condordem120), na.rm=TRUE), cell.sz, 25)

# note: we use a raster coarse integration time step so the
# example runs faster
dens.res <- initializeDensity(mkde.obj, mv.dat, integration.step=10.0)
mkde.obj <- dens.res$mkde.obj
mv.dat <- dens.res$move.dat

# Save as VTK file
writeToVTK(mkde.obj, "condor_3dMKDE.vtk",
description="Example California condor 3D MKDE")

# Clean up files
unlink("condor_3dMKDE.vtk")
}
