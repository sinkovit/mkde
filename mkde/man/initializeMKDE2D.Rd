\name{initializeMKDE2D}
\alias{initializeMKDE2D}
\title{
Set up a 2D MKDE object.
}
\description{
Define the spatial extent and resolution of a 2D MKDE and create an 2D MKDE list object for use in other functions in the package. 
}
\usage{
initializeMKDE2D(xLL, xCellSize, nX, yLL, yCellSize, nY)	
}
\arguments{
  \item{xLL}{Lower bounds of the grid in the x-dimension}
  \item{xCellSize}{Cell size in the x-dimension}
  \item{nX}{Number of cells in the x-dimension}
  \item{yLL}{Lower bounds of the grid in the y-dimension}
  \item{yCellSize}{Cell size in the y-dimension}
  \item{nY}{Number of cells in the y-dimension}
}
\details{
It is strongly recommended that the same value is used for \code{xCellSize} and \code{yCellSize}.  The grid should be defined so that it covers the area that the animal used, plus a sufficient buffer so that the density is negligable beyond the grid.
}
\value{
A list representing an MKDE object is returned with the following elements:
  \item{dimension}{The dimension of the MKDE; that is, 2.}
  \item{x}{A grid of points along the x-axis where the cell centers occur.}
  \item{y}{A grid of points along the y-axis where the cell centers occur.}
  \item{z}{A grid of points along the z-axis where the cell centers occur.  For a 2D MKDE z = NA.}
  \item{z.min}{A 2D array representing the lower bounds of space in the z-dimension at each x and y coordinate.  Defaults to -Inf.}
  \item{z.max}{A 2D array representing the upper bounds of space in the z-dimension at each x and y coordinate. Defaults to Inf.}
  \item{nx}{Number of cells in the x-dimension.}
  \item{ny}{Number of cells in the y-dimension.}
  \item{nz}{Number of cells in the z-dimension. For a 2D MKDE nz = 1.}
  \item{d}{A 2D array with dimensions (nx, ny) that stores the density. The elements are initialized to NA.}
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
cell.sz <- mean(res(pandadem))
ext <- ext(pandadem)
nx <- ncol(pandadem)
ny <- nrow(pandadem)
mkde.obj <- initializeMKDE2D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz, ny)
}
