#
# function to read simcoal output and integrate into rmetasim
#
#

#
# Mark's functions to parse arlequin files
#

`%upto%` <- function (from, to) 
if (from <= to) from:to else numeric(0)

coal2rmet.old <- function( file, norm=TRUE){
  dat <- parse.arlequin( file)
  pops <- unique( dat$pop)
  is.nuc <- !is.factor( dat[,2])
  dat[,-1] <- lapply( dat[,-1,drop=F], as.character)
  if( is.nuc) {
    nr <- nrow( dat)
    dat <- dat[ rep( 1:nr, each=2),]
    nl <- (ncol( dat)-1)/2
    dat[ 2*(1:nr), (1:nl)*2] <- dat[ 2*(1:nr)-1, (1:nl)*2+1]
    dat <- dat[ ,-2*(1:nl)-1]
  }
    
  dat[,-1] <- lapply( dat[,-1,drop=F], as.factor)
  n.locs <- ncol( dat)-1
  loctable <- lapply( dat[,-1,drop=F], function( x) {
    m <- t( table( dat[,1], as.integer( x)))
    dimnames( m)[[1]] <- levels( x)
    m
  })
  
  if( norm)
    loctable <- lapply( loctable, function( x) x / rep( colSums( x), each=nrow( x)))
  loctable
    
}

parse.arlequin.old <- function( file){
  dat <- readLines( file)
  data <- NULL
  i.pop <- 0
  while( !is.na( droppo <- grep( 'SampleData *=', dat)[1])) {
    i.pop <- i.pop + 1
    dat <- dat[ -(1:droppo)]
    endo <- grep( '^ *\\} *$', dat)[1]
    dati <- dat[ 1 %upto% (endo-1)]
    blanks <- grep( '^( ||t)*$', dati)
    if( length( blanks))
      dati <- dati[ -blanks]
    is.indiv <- regexpr( '_', dati)>0
    if( is.nuc <- !all( is.indiv))
      dati <- paste( dati[ is.indiv], dati[ !is.indiv], sep=' ')

    fo <- textConnection( dati)
    on.exit( close( fo))
    data.i <- cbind( pop=i.pop, read.table( fo, header=FALSE, row=NULL)[ ,-(1:2),drop=FALSE])
    close( fo)
    on.exit()
    
    if( is.null( data))
      data <- data.i
    else
      data <- rbind( data, data.i)
  }
  
  if( is.nuc) {
    nl <- (ncol( data)-1)/2
    data <- data[ ,c (1, 1+c( matrix( 1:(2*nl), nrow=2, byrow=T)))]
  }
  data
}

# EG: coal2rmet( 'tossm_0.arp')

#
#
# Allan's function to glue the parsed data into rmetasim format
# This function will clobber existing individuals and loci in rland object
#
#
landscape.coalinput.old <- function(rland, npp=200, arlseq = 'seq.arp', seqsitemut= 1e-7, arlms = 'ms.arp', msmut= 5e-4)
  {
    if (is.null(arlseq)&is.null(arlms))
      {
        print ("you must specify some type of coalescent input")
        rland
      } else #give it a shot
      {
        if (!is.null(arlseq))
          clocseq <- coal2rmet(arlseq)
        else
          clocseq <- NULL

###
### temporary hack to make imported sequences work
###
        if(!is.null(clocseq))
          {
            clocseq <- clocseq[1]
          }
###
###
        
        if (!is.null(arlms))
          clocms <- coal2rmet(arlms)
        else
          clocms <- NULL

        cloc <- c(clocseq,clocms)
        
        for (loc in 1:length(cloc))
          {
            states <- rownames(cloc[[loc]])
            ltype <- c(1,2)[length(grep("T",states[1]))+1]
            freqs <- apply(as.matrix(cloc[[loc]]),1,mean)
            if (ltype==1)
              {
                rland <- landscape.new.locus(rland, type=ltype, ploidy=2, mutationrate = msmut,
                                   numalleles=length(states), frequencies = freqs,
                                   states = as.numeric(states), transmission = 0)
              } else {
                rland <- landscape.new.locus(rland, type=ltype, ploidy=1, mutationrate = seqsitemut,
                                   numalleles=length(states), frequencies = freqs,
                                   states = states, allelesize=nchar(states[1]),
                                   transmission = 1)
              }
          }
        
                                        #
                                        #
        
        S <- rland$demography$localdem[[1]]$LocalS  #needs to be changed to localdemk
        R <- rland$demography$localdem[[1]]$LocalR  #needs to be changed to localdemk
                                        #here are the eigenvectors for the local demographies
        ev <- eigen((R+diag(dim(R)[1]))%*%S)$vectors[,1]

        if (length(npp)==1)
          NperPop <- rep(npp,rland$intparam$habitats) #population size per population, make into a parameter
        else
          NperPop <- npp
        indlist <- vector("list",length(NperPop))
        
        for (p in 1:length(NperPop))
          {
            im <- matrix(NA,nrow=NperPop[p],ncol=3+sum(landscape.ploidy(rland)))
            colcnt <- 4
            for (loc in 1:length(landscape.ploidy(rland)))
              {
                probs <- cloc[[loc]][,p]
                indices <- sapply(rland$loci[[loc]]$alleles,function(x){x$aindex})
                for (al in 1:landscape.ploidy(rland)[loc])
                  {
                    im[,colcnt] <- sample(indices,NperPop[p],replace=T,prob=probs)
                    colcnt <- colcnt+1
                  }
              }
            im[,1] <- sample((0:(rland$intparam$stages-1))+((p-1)*rland$intparam$stages),NperPop[p],replace=T,prob=ev)
            im[,2] <- rep(0,NperPop[p])
            im[,3] <- im[,2]
            indlist[[p]] <- im
          }
        individuals <- do.call("rbind",indlist)
        individuals <- individuals[order(individuals[,1]),]
        rland$individuals <- matrix(as.integer(individuals),nrow=dim(individuals)[1])

        rland
      }
  }
