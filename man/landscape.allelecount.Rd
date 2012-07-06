\name{landscape.allelecount}
\alias{landscape.allelecount}
\title{Calculate allele numbers (frequency in the statistical sense) at each locus in each population}

\description{
  Calculate allele counts}
\usage{
  hetmat <- landscape.exp.het(rland,tbl.out=F)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{tbl.out}{return as a (three-dimensional) table if TRUE.  If
  FALSE, return as a dataframe with categorical variables denoting the
  locus, population and allele.}
}
\value{
  Depends on the value of tbl.out.  See above.
}
\seealso{
landscape.allelefreq, landscape.obs.het, landscape.exp.het, landscape.Fwright, landscape.Fst
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.allelefreq(exampleland,tbl.out=TRUE)
  landscape.allelefreq(exampleland,tbl.out=FALSE)
  rm(exampleland)
}
\keyword{misc}
