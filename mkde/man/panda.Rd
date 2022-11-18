\name{panda}
\alias{panda}
\docType{data}
\title{
Giant panda locations
}
\description{
A data frame containing giant panda location data
}
\usage{data(panda)}
\format{
  A data frame with 147 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{Elapsed time in minutes}
    \item{\code{x}}{x-coordinate (UTM easting in meters)}
    \item{\code{y}}{y-coordinate (UTM northing in meters)}
    \item{\code{z}}{z-coordinate (height above sea level in meters)}
  }
}
\details{
GPS location data acquired from a wild giant panda (Ailuropoda
melanoleuca) tracked by the Chinese Academy of Sciences and San Diego
Zoo Global in Foping National Nature Reserve, China during August 2008
(Zhang et. al. 2012).
}
\source{
Giant panda (Ailuropoda melanoleuca) biotelemetry data was collected by
Fuwen Wei and others of the Key Laboratory of Animal Ecology and
Conservation Biology, Institute of Zoology, Chinese Academy of Science,
People's Republic of China in collaboration with the San Diego Zoo
Institute for Conservation Research. Giant panda research was funded by
the National Natural Science Foundation of China (31230011), Wildlife
Experimental Platform of Chinese Academy of Sciences, and San Diego Zoo
Global.
}
\references{
Zhang, Z., Sheppard, J.K., Swaisgood, R.R., Wang, G., Nie, Y., Wei, W., Zhao, N. & Wei, F. (2014). Ecological Scale and Seasonal Heterogeneity in the Spatial Behaviors of Giant Pandas. Integrative Zoology, 9: 46-60.
}
\examples{
data(panda)
head(panda, 30)
}
\keyword{datasets}
