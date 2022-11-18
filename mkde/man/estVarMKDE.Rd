\name{estVarMKDE}
\alias{estVarMKDE}
\title{
Estimate move variances for 3D MKDE.
}
\description{
Estimates the move variance for the (x, y) and z dimensions and altitude for a
Brownian bridge MKDE.
}
\usage{
estVarMKDE(move.dat)	
}
\arguments{
  \item{move.dat}{A move data object created with \code{initializeMovementData}}
}
\details{
  Computes maximum-likelihood estimates for move variances. If the MKDE is 2D or 2.5D, only the variance for the xy-dimension is estimated.  If the MKDE is 3D, a separate move variance for the z-dimension is also estimated.
}
\value{
  A move data object updated with move variances
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

mv.dat <- estVarMKDE(mv.dat)
cat(paste("Move variance in xy-dimensions = ",
mv.dat$sig2xy[1], "\n", sep=""))
cat(paste("Move variance in z-dimension = ",
mv.dat$sig2.z[1], "\n", sep=""))
}
