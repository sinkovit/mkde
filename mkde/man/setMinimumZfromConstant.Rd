\name{setMinimumZfromConstant}
\alias{setMinimumZfromConstant}
\title{
Set minimum z-axis value to a constant.
}
\description{
Set the lower bounds in the z-dimension for each location in the x and y
dimensions to a constant value.
}
\usage{
setMinimumZfromConstant(mkde.obj, val)	
}
\arguments{
  \item{mkde.obj}{ 2D or 3D MKDE object created with
    \code{initialize3DMKDE} or \code{initialize3DMKDE}, respectively}
  \item{val}{The value at which the lower bound should be set for all
    locations in the x and y dimensions.}
}
\details{
Obviously, the lower bound must be less than the upper bound.
}
\value{
An updated MKDE list object is returned.
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
library(raster)
data(dugongdem)
cell.sz <- mean(res(dugongdem))
ext <- extent(dugongdem)
nx <- ncol(dugongdem)
ny <- nrow(dugongdem)
mkde.obj <- initializeMKDE3D(ext@xmin, cell.sz, nx, ext@ymin, cell.sz,
ny, min(values(dugongdem), na.rm=TRUE), 50.0, 15)
mkde.obj <- setMinimumZfromConstant(mkde.obj, -20.0)
}
