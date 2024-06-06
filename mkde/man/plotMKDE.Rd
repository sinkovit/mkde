\name{plotMKDE}
\alias{plotMKDE}
\title{
Make an image plot of an MKDE.
}
\description{
Makes an image plot of an MKDE.  If the MKDE is 3D, the contours will be
based on the entire MKDE, but only one level indexed by z.index will be
plotted.
}
\usage{
plotMKDE(mkde.obj, z.index=1, probs=c(0.99, 0.95, 0.90, 0.75, 0.5, 0.0),
cmap=rev(rainbow(length(probs)-1)), add=FALSE, ...)	
}
\arguments{
  \item{mkde.obj}{A 2D, 2.5D, or 3D MKDE object created with
    \code{initialize2DMKDE} or \code{initialize3DMKDE} and density
    initialized with \code{initializeDensity}}
  \item{z.index}{Index for the z-dimension of the density if a 3D MKDE}
  \item{probs}{Probabilities for image contours.}
  \item{cmap}{Color map for image plot.}
  \item{add}{FALSE to make a new plot, TRUE to add to existing plot}
  \item{...}{Additional graphical parameters.}
}
\details{
A plot of the density for the specified level in the z-dimension (which
should be 1, the defaul value, for a 2D or 2.5D MKDE) is generated.
}
\value{
No value is returned.
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
fpath <- system.file("extdata", "pandadem.RDS", package="mkde")
pandadem <- terra::readRDS(fpath)
cell.sz <- 0.25*mean(res(pandadem))
ext <- ext(pandadem)
nx <- 4*ncol(pandadem)
ny <- 4*nrow(pandadem)
mkde.obj <- initializeMKDE2D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz, ny)

# set up movement data
data(panda)
mv.dat <- initializeMovementData(panda$time, panda$x, panda$y, 
sig2obs=25.0, t.max=185.0)

# estimate density
dens.res <- initializeDensity(mkde.obj, mv.dat)
mkde.obj <- dens.res$mkde.obj
mv.dat <- dens.res$move.dat
plotMKDE(mkde.obj)
}
