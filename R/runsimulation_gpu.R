
landscape.simulate.gpu <- function(Rland, numit, seed=-1, compress=FALSE, adj.lambda=0)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- landscape.coerce(Rland)
        .Call("iterate_landscape_gpu",as.integer(numit),Rland,as.integer(compress),as.integer(adj.lambda),PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }
