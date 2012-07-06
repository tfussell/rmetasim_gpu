\name{landscape.read}
\alias{landscape.read}
\title{Load a landscape from a file}

\description{
  Load a Rmetasim landscape from a file into an R object.}
\usage{
  rland <- landscape.read(fn="mylandscape.lnd")}
\arguments{
  \item{fn}{the path and name of the file containing the landscape}
}

\examples{
  ## Needs write access to the current directory, files created!!
  landscape.write(landscape.new.example(), "exampleland.lnd")
  exampleland <- landscape.read("exampleland.lnd")
  exampleland
  rm(exampleland)
}
\keyword{misc}
