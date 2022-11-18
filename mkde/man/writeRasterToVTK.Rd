\name{writeRasterToVTK}
\alias{writeRasterToVTK}
\title{
Write a 2D raster to XDMF XML wrapper and binary data file.
}
\description{
Write the raster to a XDMF files.
}
\usage{
writeRasterToVTK(elev, r.rst, g.rst, b.rst, descr, fname)	
}
\arguments{
  \item{elev}{A RasterLayer object}
  \item{r.rst}{A RasterLayer object for red}
  \item{g.rst}{A RasterLayer object for green}
  \item{b.rst}{A RasterLayer object for blue}
  \item{descr}{String description to be added to header of VTK file}
  \item{fname}{The path and base file name for output HDF5 files}
}
\details{
  This function writes a raster to VTK format. The raster is colored
  according to the RGB values in \code{r.rst}, \code{g.rst}, and
  \code{b.rst}, respectively. The RGB balues must be an interger from 0
  to 255.
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
library(raster)
data(condordem120)
elev.val <- values(condordem120)
elev.min <- min(elev.val, na.rm=TRUE)
elev.max <- max(elev.val, na.rm=TRUE)

# make a color lookup table
cmap <- data.frame(value=c(0.0, 0.25, 0.5, 0.75, 1.0), 
  R=c(150, 179, 205, 192, 252), 
  G=c(224, 204, 205, 183, 243), 
  B=c(94, 147, 168, 147, 226))
                    
cmap$value <- cmap$value*(elev.max - elev.min) + elev.min
# red
f.R <- approxfun(cmap$value, cmap$R)
elev.r <- raster(condordem120)
values(elev.r) <- round(f.R(elev.val))
# green
f.G <- approxfun(cmap$value, cmap$G)
elev.g <- raster(condordem120)
values(elev.g) <- round(f.G(elev.val))
# blue
f.B <- approxfun(cmap$value, cmap$B)
elev.b <- raster(condordem120)
values(elev.b) <- round(f.B(elev.val))
writeRasterToVTK(condordem120, elev.r, elev.g, elev.b, "Elevation for
California condor Example", "condor_dem.vtk")
}
