\name{dugong}
\alias{dugong}
\docType{data}
\title{
Dugong locations
}
\description{
A data frame containing dugong location data
}
\usage{data(dugong)}
\format{
  A data frame with 426 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{Elapsed time in minutes}
    \item{\code{x}}{x-coordinate (UTM easting in meters)}
    \item{\code{y}}{y-coordinate (UTM northing in meters)}
    \item{\code{z}}{z-coordinate (Depth in meters)}
  }
}
\details{
GPS location data acquired from a wild dugong (Dugong dugon) tracked by James Cook University in Hervey Bay, Australia during August 2004 (Sheppard et. al. 2010)
}
\source{
James Sheppard, San Diego Zoo Institute for Conservation Research
}
\references{
Sheppard, J. Marsh, H., Jones, R.E., & Lawler, I.R. (2010). Dugong habitat use in relation to seagrass nutrients, tides, and diel cycles. Marine Mammal Science 26, 855-879.
}
\examples{
data(dugong)
head(dugong, 30)
}
\keyword{datasets}
