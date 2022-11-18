\name{initializeMovementData}
\alias{initializeMovementData}
\title{
Initialize a movement data list
}
\description{
This function sets up the movement data for use in other functions in
the package.
}
\usage{
initializeMovementData(t.obs, x.obs, y.obs, z.obs=NULL, sig2obs=0.0,
sig2obs.z=NA, t.max=max(diff(t.obs), na.rm=TRUE))
}
\arguments{
  \item{t.obs}{A vector of times at which the animal locations were observed}
  \item{x.obs}{A vector of x-coordinates at which the animal locations
    were observed}
  \item{y.obs}{A vector of y-coordinates at which the animal locations
    were observed}
  \item{z.obs}{A vector of z-coordinates at which the animal locations
    were observed}
  \item{sig2obs}{Location error variance in the xy-dimensions}
  \item{sig2obs.z}{Location error variance in the z-dimension}
  \item{t.max}{The maximum time allowed between locations if the
    movement step is to be used in computing the density}
}
\details{
If only 2D or 2.5D MKDEs are to be calculated, then z.obs and sig2obs.z
do not have to be provided.
}
\value{
A move data object, in the form of a list, is returned.
  \item{dimension}{The spatial dimension of the movement data.  If
    z-coordinates are passed, the dimension will be 3; otherwise, the
    dimension will be 2.}
  \item{t.obs}{Times of the animal locations.}
  \item{x.obs}{x-coordinates of the animal locations.}
  \item{y.obs}{y-coordinates of the animal locations.}
  \item{z.obs}{z-coordinates of the animal locations.}
  \item{a.obs}{Altitude of the animal; that is, its z-coordinate minus
    the lower bound on the z-coordinate at the corresponding
    xy-coordinates of the animal location.}
  \item{use.obs}{A logical array that indicates whether each location
    should be used in the MKDE calculations.  By default these values
    are set to TRUE.}
  \item{t.max}{The maximum time allowed between locations if the
    movement step is to be used in computing the density}
  \item{sig2xy}{A vector of movemenet variance paramters for the
    xy-dimensions, with NAs as placeholders.  The functions estVarMKDE
    is provided to estimate these parameters, or they can be set
    manually to allow for different movement variances for each move
    step.}
  \item{sig2z}{A vector of movemenet variance paramters for the
    z-dimension, with NAs as placeholders.  The functions estVarMKDE is
    provided to estimate these parameters, or they can be set manually
    to allow for different movement variances for each move step.}
  \item{sig2obs}{A vector of variances for the location observation
    error in the xy-dimensions.  If only one value is provided for the
    variance parameters, the value is repeated to form a vector with one
    element for each location.  Otherwise, a vector of location error
    variances can be passed to allow for different errors for each
    observed location.}
  \item{sig2obs.z}{A vector of variances for the location observation
    error in the z-dimension.  If only one value is provided for the
    variance parameters, the value is repeated to form a vector with one
    element for each location.  Otherwise, a vector of location error
    variances can be passed to allow for different errors for each
    observed location.}
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
data(condor)
mv.dat <- initializeMovementData(condor$time, condor$x,
condor$y, z.obs=condor$z, 65.0, 25.0, sig2obs.z=81.0)
}
