\name{landscape.write}
\alias{landscape.write}
\title{Save a landscape to a file}

\description{
  Save a Rmetasim landscape object to a file}
\usage{
  rland <- landscape.write(rland=l,fn="mylandscape.lnd")
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{fn}{the path and name of the file to save the landscape to}
}

\examples{
  ## Needs write access to the current directory, files created!!
  exampleland <- landscape.new.example()
  landscape.write(exampleland, "exampleland.lnd")

  rm(exampleland)
}
\keyword{misc}
