\name{landscape.write.foreign}
\alias{landscape.write.foreign}
\title{Save a landscape to a file in a foreign format}

\description{
  Save a Rmetasim landscape object to a file in a suite of output formats
}
\usage{
  rland <- landscape.write.foreign(rland,fn="foreign",numi=24,fmt="GDA")
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{fn}{the path and name of the file to save the landscape
    to}
  \item{numi}{number of individuals sampled per population for
      inclusion in subsequent analyses}
  \item{fmt}{the output format for the landscape:  Can take the
    following values:"arlequin","arlequinhap","biosys","genpop","R","migrate","rerat","GDA"}
}

\examples{
  ## Needs write access to the current directory, files created!!
  exampleland <- landscape.new.example()
  landscape.write.foreign(exampleland, fn="exampleland.nex", fmt="GDA")
  rm(exampleland)
}
\keyword{misc}
