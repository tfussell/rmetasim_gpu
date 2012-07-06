#returns a landscape with a sample of populations and or individuals within populations
#np is the number of populations, if NULL then all populations, if a vector, gives the populations to sample
#if np is a single number then the number of populations to sample at random
landscape.sample <- function(rland,np=NULL,ns=NULL,pvec=NULL)
  {
    
    if (!is.null(np))
      {
        pops <- sample(unique(landscape.populations(rland)),
                       ifelse(length(unique(landscape.populations(rland)))>np,np,
                              length(unique(landscape.populations(rland)))),
                       ,replace=F)
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pops,]
      }
    else if (!is.null(pvec))
      {
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pvec,]
      }
    if (!is.null(ns))
      {
        pops <- landscape.populations(rland)
        ptbl <- table(pops)
        if ((is.null(pvec) & is.null(np)))
          {
            names <- as.numeric(names(ptbl))
           } else {
             names <- as.numeric(names(which(ptbl>=np)))
           }
        rland$individuals <- rland$individuals[
                                               as.numeric(unlist(sapply(unique(names),
                                                                        function(x,pops,ns)
                                                                        {
                                                                          ss <- ifelse(length(which(pops==x))>ns,
                                                                                       ns,length(which(pops==x)))
                                                                          sample(which(pops==x),ss,replace=F)
                                                                        },pops=pops,ns=ns)))
                                               ,]
      }
    rland
  }


#returns a landscape with a sample of individuals within populations
#
#ns should be a vector of sample sizes to take from each population
#landscape.sample.variable <- function(rland,ns=NULL)
#  {
#    if (!is.null(ns))
#      {
#        pops <- landscape.populations(rland)
#        ptbl <- table(pops)
#        names <- as.numeric(names(which(ptbl>=np)))
#        rland$individuals <-
#          rland$individuals[
#                            as.numeric(unlist(sapply(unique(names),
#                                                     function(x,pops,ns)
#                                                     {
#                                                       ss <- ifelse(length(which(pops==x))>ns[x],
#                                                                    ns[x],length(which(pops==x)))
#                                                       sample(which(pops==x),ss,replace=F)
#                                                     },pops=pops,ns=ns)))
#                            ,]
#      }
#    rland
#  }
