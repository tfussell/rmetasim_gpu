\name{landscape.to.rtheta}
\alias{landscape.to.rtheta}
\title{Converts genetic marker data in a landscape into a format
suitable for analysis using ANOVA}

\description{
This function converts Rmetasim landscapes into a format that can be
analyzed with ANOVA to calculate Weir's theta.
}
\usage{
  genin <- landscape.to.rtheta(Rland,numi)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
  \item{numi}{number of individuals to sample at random from each
  population/habitat.  The default value of 0 takes all individuals in
  each population.  CAUTION this could take a long time.  }
}
\value{
A matrix with columns representing the population, class, individual,
locus, allelic position, allele id.
}
\seealso{
popstruct
}

\keyword{misc}
