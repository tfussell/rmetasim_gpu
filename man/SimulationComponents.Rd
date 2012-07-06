\name{SimulationComponents}
\alias{landscape.advance}
\alias{landscape.carry}
\alias{landscape.extinct}
\alias{landscape.reproduce}
\alias{landscape.survive}
\title{Code components to simulate a landscape}

\description{
  These functions can be used to construct custom simulations of
  landscapes.  Each conducts only a single generations worth of change}
\usage{
  rland <- landscape.advance(rland=l)
  rland <- landscape.carry(rland=l)
  rland <- landscape.extinct(rland=l)
  rland <- landscape.reproduce(rland=l)
  rland <- landscape.survive(rland=l)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{
  \code{landscape.advance()} merely advances the generation counter and selects
  the new generations demographic conditions if such conditions can
  vary.  The other functions implement carrying capacity, local
  extinction, reproduction, and survival/growth, respectively.  The
  function \code{landscape.simulate()} bundles the functionality of these
  components into a single function (and executes it slightly faster all within linked C++ code).
}

\seealso{
  \code{landscape.simulate}
  }

\keyword{misc}
