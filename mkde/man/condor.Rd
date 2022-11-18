\name{condor}
\alias{condor}
\docType{data}
\title{
California condor locations
}
\description{
A data frame containing California condor location data
}
\usage{data(condor)}
\format{
  A data frame with 421 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{Elapsed time in minutes}
    \item{\code{x}}{x-coordinate (UTM easting in meters)}
    \item{\code{y}}{y-coordinate (UTM northing in meters)}
    \item{\code{z}}{z-coordinate (height above sea level in meters)}
  }
}
\details{
GPS location data acquired from a wild California condor (Gymnogyps californianus) tracked by San Diego Zoo Global around its reintroduction site in Baja California, Mexico during November 2012 (Sheppard et al. 2013).
}
\source{
James Sheppard, San Diego Zoo Institute for Conservation Research
}
\references{
Sheppard, J.K., Walenski, M., Wallace, M.P., Velazco, J.J.V., Porras, C., & Swaisgood, R.R. (2013). Hierarchical dominance structure in reintroduced California condors: correlates, consequences, and dynamics. Behavioral Ecology and Sociobiology. 67: 1-12.
}
\examples{
data(condor)
head(condor, 30)
}
\keyword{datasets}
