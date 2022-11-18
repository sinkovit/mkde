\name{writeToXDMF}
\alias{writeToXDMF}
\title{
Write 3D MKDE to XDMF XML wrapper and binary data file.
}
\description{
Write the MKDE to a XDMF files.
}
\usage{
writeToXDMF(mkde.obj, fname, nodat="NA", cumprob=FALSE)	
}
\arguments{
  \item{mkde.obj}{3D MKDE object created with \code{initialize3DMKDE}
    and density initialized with \code{initializeDensity}}
  \item{fname}{The path and base file name for output XDMF files}
  \item{nodat}{A no data character string that will be written in place
    of no data values.}
  \item{cumprob}{Indicate whether to write the voxel probabilities of
    cumulative probabilities.}
}
\details{
  This function writes XDMF XML wrapper and binary data file.
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
  \email{jsheppard@sandiegozoo.org} \cr
  Amit Chourasia, MS \cr
  San Diego Supercomputer Center \cr
  \email{amitsdsc.edu}
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

# Write file (note no file extension!)
writeToXDMF(mkde.obj, "condor_3dMKDE")
}
