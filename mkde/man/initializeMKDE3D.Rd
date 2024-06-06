\name{initializeMKDE3D}
\alias{initializeMKDE3D}
\title{
Set up a 3D MKDE object.
}
\description{
Define the spatial extent and resolution of a 3D MKDE and create an 3D
MKDE list object for use in other functions in the package.
}
\usage{
initializeMKDE3D(xLL, xCellSize, nX, yLL, yCellSize, nY,
zLL, zCellSize, nZ)	
}
\arguments{
  \item{xLL}{Lower bounds of the grid in the x-dimension}
  \item{xCellSize}{Cell size in the x-dimension}
  \item{nX}{Number of cells in the x-dimension}
  \item{yLL}{Lower bounds of the grid in the y-dimension}
  \item{yCellSize}{Cell size in the y-dimension}
  \item{nY}{Number of cells in the y-dimension}
  \item{zLL}{Lower bounds of the grid in the z-dimension}
  \item{zCellSize}{Cell size in the z-dimension}
  \item{nZ}{Number of cells in the z-dimension}
}
\details{
It is strongly recommended that the same value is used for \code{xCellSize} and
\code{yCellSize}.  The grid should be defined so that it covers the volume that
the animal used, plus a sufficient buffer so that the density is
negligable beyond the grid.
}
\value{
A list representing an MKDE object is returned with the following elements:
  \item{dimension}{The dimension of the MKDE; that is, 3.}
  \item{x}{A grid of points along the x-axis where the voxel centers occur.}
  \item{y}{A grid of points along the y-axis where the voxel centers occur.}
  \item{z}{A grid of points along the z-axis where the voxel centers occur.}
  \item{z.min}{A 2D array representing the lower bounds of space in the
    z-dimension at each x and y coordinate.  Defaults to a constant
    value based on the arguments zLL, zCellSize,and nZ.}
  \item{z.max}{A 2D array representing the upper bounds of space in the
    z-dimension at each x and y coordinate. Defaults to a constant value
    based on the arguments zLL, zCellSize,and nZ.}
  \item{nx}{Number of cells in the x-dimension.}
  \item{ny}{Number of cells in the y-dimension.}
  \item{nz}{Number of cells in the z-dimension.}
  \item{d}{A 3D array with dimensions (nx, ny, nz) that stores the
    density. The elements are initialized toNA.}
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
fpath <- system.file("extdata", "condordem.RDS", package="mkde")
condordem <- terra::readRDS(fpath)
cell.sz <- mean(res(condordem))
ext <- ext(condordem)
nx <- ncol(condordem)
ny <- nrow(condordem)
zmin <- min(values(condordem), na.rm=TRUE)
nz <- 30
mkde.obj <- initializeMKDE3D(ext$xmin, cell.sz, nx, ext$ymin, cell.sz,
ny, zmin, cell.sz, nz)
}
