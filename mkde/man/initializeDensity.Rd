\name{initializeDensity}
\alias{initializeDensity}
\title{
Calculate raster of density values for MKDE.
}
\description{
Given a movement data object (list) and MKDE object (list), estimate the
movement variance parameters and update them in the movement data object
and then compute the density for a 2D, 2.5D, or 3D MKDE.
}
\usage{
initializeDensity(mkde.obj, move.dat, integration.step, d.thresh=1e-25) 	
}
\arguments{
  \item{mkde.obj}{An MKDE object created with \code{initialize2DMKDE} or
    \code{initialize3DMKDE}}
  \item{move.dat}{A move data object created with \code{initializeMovementData}}
  \item{integration.step}{A time step used for numerical integration
    over the movement trajectory}
  \item{d.thresh}{The value of the kernel below which its contibrution
    to the overal density is considered negligible}
}
\details{
This function computes the density for a 2D, 2.5D, or 3D MKDE.  If the
move variance parameters have not been estimated, they will be prior to
computing the density.  The integration time step should be much less
than the maximum time step allowed between observed animal
locations. For a 2.5D MKDE, if the area has been calculated with
\code{initializeAreaRaster}, then the cell probabilities are weighted by
the proportion of the total surface area represented by the cell.
}
\value{
A list is returned with two elements: \code{move.dat} and \code{mkde.obj}.  The first
is an updated move data object and the second is an updated MKDE object.
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
  \email{glock@sdsc.edu}
}
\examples{
library(terra)
data(condor)
condor <- condor[1:20,] # simply to make example run more quickly
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

fpath <- system.file("extdata", "condordem120.RDS", package="mkde")
condordem120 <- terra::readRDS(fpath)
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
}
