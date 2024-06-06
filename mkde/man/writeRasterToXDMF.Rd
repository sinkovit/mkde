\name{writeRasterToXDMF}
\alias{writeRasterToXDMF}
\title{
Write a 2D raster to XDMF XML wrapper and binary data file.
}
\description{
Write the raster to a XDMF files.
}
\usage{
writeRasterToXDMF(rast, fname, nodat="NA")	
}
\arguments{
  \item{rast}{A RasterLayer object}
  \item{fname}{The path and base file name for output HDF5 files}
  \item{nodat}{A no data character string that will be written in place
    of no data values.}
}
\details{
  This function writes XDMF XML wrapper and binary data file.
}
\value{
  No value is returned
}
\author{
  Jeff A. Tracey, PhD \cr
  USGS Western Ecological Research Center, San Diego Field Station \cr
  \email{jatracey@usgs.gov}  \cr
  James Sheppard, PhD \cr
  San Diego Zoo Institute for Conservation Research \cr
  \email{jsheppard@sandiegozoo.org} \cr
  Amit Chourasia, MS \cr
  San Diego Supercomputer Center \cr
  \email{amit@sdsc.edu}
}
\examples{
library(terra)
fpath <- system.file("extdata", "condordem.RDS", package="mkde")
condordem <- terra::readRDS(fpath)

# Save as XDMF (notice no file extension in file name)
writeRasterToXDMF(condordem, "condor_dem")

# Clean up files
unlink("condor_dem.dat")
unlink("condor_dem.xdmf")
}
