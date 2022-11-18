\name{mkde-package}
\alias{mkde-package}
\alias{mkde}
\docType{package}
\title{
Movement-based kernel density estimates (MKDEs) in 2 or 3 spatial dimensions.
}
\description{
The mkde package enables animal space use to be estimated in three
dimensions (3D) using data collected from biotelemetry tracking
devices. This package addresses a recognized need in modeling animal
space use (Belant et al. 2012) wherein researchers have been limited by
the lack of 3D home range estimators. Animal space use can be
characterized by the (x, y) spatial dimensions as well as a third
z-dimension representing altitude, elevation, or depth for flying,
terrestrial, or aquatic species, respectively. Although many
biotelemetry devices record 3D location data with x, y, and z
coordinates from tracked animals, the third z coordinate is typically
not integrated into studies of animal spatial use. The mkde package
enables users to visually explore the 3D MKDE volumes of animals to more
intuitively understand how they are spatially related to the
environmental covariates and bounding layers within their ranges, such
as bathymetry or topography.

The mkde package builds on previous work on the Brownian bridge
approach for estimating animal utilization distributions (Horne et
al. 2007). This method, in contrast to location-based KDEs, integrates
kernels over time along a movement path interpolated between observed
locations. Benhamou distinguished location-based kernel density
estimators (LKDE) from movement-based kernel density estimators (MKDE),
which includes Brownian bridge and biased random walk models. MKDEs
account for time between consecutively observations in the estimator, do
not requiring independent samples from the UD, and thus more
realistically represent the space used by an animal.

The user inputs animal location data typically obtained by a Global
Positioning System (GPS) or Very High Frequency (VHF) device in which
each observation includes an x-coordinate, a y-coordinate,a
z-coordinate, and time. The observed locations are assumed to be subject
to observation error and are normal random variables. The observation
error variances are either provided by the manufacturers of the
telemetry equipment or estimated from field trials, e.g., Hansen and
Riggs (2008). Often, an animal's movement is limited in the
z-dimension. For example, avian species are generally bounded below by
the earth's surface, whereas marine animals are bounded below by the sea
floor and above by the water's surface. Package functions allow the mkde
user to bound the density in the z-dimension by a(x,y) and b(x,y) with a
constant or a 2D raster.

The mkde package provides a 2.5D approach for computing home range area
that essentially uses a 2D MKDE draped over a 2D elevation raster. The
bias is corrected by calculating and summing the surface area of each
cell of the elevation raster that falls within a desired probability
contour of the 2D MKDE. An algorithm developed by Jenness (2004, 2014)
is used to compute the surface area of each raster cell. This method
uses the cell center coordinates and elevations of the focal cell and
its eight neighboring cells to construct eight triangular facets within
the focal cell. The area of each facet is calculated using Heron's
formula and then summed to obtain the surface area for the focal cell.

Numerous functions are provided to write output files in various formats
(VTK, XDMF, ASCII) for use in other GIS and 3D Visualization
applications.
}
\details{
\tabular{ll}{
Package: \tab mkde\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-23\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
  Jeff A. Tracey (US Geological Survey, San Diego Field Station, Western
  Ecological Research Center) \cr
  James Sheppard (San Diego Zoo Institute for Conservation Research) \cr
  Jun Zhu (Department of Statistics and Department of Entomology,
  University of Wisconsin -- Madison) \cr
  Robert Sinkovits (San Diego Supercomputer Center) \cr
  Amit Chourasia (San Diego Supercomputer Center) \cr
  Glenn Lockwood (San Diego Supercomputer Center) \cr
  Maintainer: Jeff A. Tracey <jatracey@usgs.gov, jeff.a.tracey@gmail.com>
}
\references{
  Tracey, J. A., Sheppard, J., Zhu, J., Wei, F., Swaisgood, R. R. and
  Fisher, R. N. (2014)
  Movement-Based Estimation and Visualization of Space Use in 3D for
  Wildlife Ecology and Conservation. PLoS ONE 9(7): e101205. doi:
  10.1371/journal.pone.0101205 \cr
  Tracy, J. A., Sheppard, J. Lockwood, G., Chourasia, A., Tatineni, M., 
  Fisher, R. N., and Sinkovits, R. (2014) Efficient 3D Movement-Based
  Kernel Density Estimator and Application to Wildlife Ecology. XSEDE 14
  Conference Proceedings, Article No. 14. doi: 10.1145/2616498.2616522\cr
  Belant, J. L., Millspaugh, J. J., Martin, J. A. & Gitzen,
  R. A. (2012). Multi-dimensional space use: the final
  frontier. Frontiers in Ecology & Environment 10, 11-12. \cr
  Benhamou, S. (2011). Dynamic Approach to Space and Habitat Use Based
  on biased random bridges. PLoS ONE 6. \cr
  Hansen, M. C., & Riggs, R. A. (2008). Accuracy, precision, and
  observation rates of global positioning system telemetry collars. The
  Journal of Wildlife Management, 72(2), 518-526. \cr
  Horne, J. S., Garton, E. O., Krone, S. M., Lewis, J. S. (2007). Analyzing
  animal movements using Brownian bridges. Ecology 88, 2354-2363. \cr
  Jenness J. S. (2004) Calculating landscape surface area from digital
  elevation models. Wildlife Society Bulletin 32: 829-839. \cr
  Jenness, J. S. (2014) Calculating landscape surface area from
  unprojected digital elevation models. In preparation.
}
\keyword{package}

