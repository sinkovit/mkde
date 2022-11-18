\name{deselectNonMovementSteps}
\alias{deselectNonMovementSteps}
\title{
Flag non-movements so they are excluded from MKDE estimation
}
\description{
This function deselects move steps where the probability that the
initial and final location are the same location is greater than or
equal to a user-defined threshold probability
}
\usage{
deselectNonMovementSteps(move.dat, p) 	
}
\arguments{
  \item{move.dat}{A move data object created with \code{initializeMovementData}}
  \item{p}{The threshold probability}
}
\details{
  If the probability that the initial and final location are the same
  location is greater than or equal to a user-defined threshold
  probability, the corresponding value in move.dat$use.obs is set to
  FALSE. Note that this function is not called within \code{initialzeDensity}.
  If you want to exclude locations because they the initial location in
  a non-movement step, this function must be used before computing the density.
}
\value{
An updated move data list object is returned.
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
mv.dat <- initializeMovementData(condor$time, condor$x, condor$y, 
z.obs=condor$z, sig2obs=25.0, sig2obs.z=81.0, t.max=65.0)

mv.dat <- deselectNonMovementSteps(mv.dat, 0.05)
}
