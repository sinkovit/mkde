\name{mkde2Dinteraction}
\alias{mkde2Dinteraction}
\title{
Probability of 2D spatial-temporal interaction.
}
\description{
  Probability of 2D spatial-temporal interaction.
}
\usage{
mkde2Dinteraction(mkde.obj, move.dat0, move.dat1, t.step, d.thresh)	
}
\arguments{
  \item{mkde.obj}{An MKDE object created with \code{initialize2DMKDE}}
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
  the same times. The cell values returned in the mkde.obj can be summed 
  to obtain a global measure of spatio-temporal interaction.
}
\value{
  Returns a list with the following elements:
  \item{mkde.obj}{An updated MKDE object containing the cell-level
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
library(raster)
data(panda)
mv.dat0 <- initializeMovementData(panda$time, panda$x, panda$y, 
sig2obs=25.0, t.max=185.0)

n <- nrow(panda)
v <- 20.0 # increase from 0 to increase difference
mv.dat1 <- initializeMovementData(panda$time, panda$x+rnorm(n, 0, v), 
panda$y+rnorm(n, 0, v), sig2obs=25.0, t.max=185.0)

data(pandadem)
cell.sz <- mean(res(pandadem))
ext <- extent(pandadem)
nx <- ncol(pandadem)
ny <- nrow(pandadem)
mkde.obj <- initializeMKDE2D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny)

res <- mkde2Dinteraction(mkde.obj, mv.dat0, mv.dat1, 10.0, 1e-20)
mkde.obj <- res$mkde.obj
mv.dat0 <- res$move.dat0
mv.dat1 <- res$move.dat1
}
