
                                        #
#These files implement some mismatch distribution calculations.  Incomplete of course.
#base differences  matrix of the sequences in the landscape
#
basediff<-function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (Rland$intparam$locusnum>=lnum)
        if (Rland$loci[[lnum]]$type==253)
          {
            sl<-landscape.locus.states(lnum,Rland);
            rmat<-matrix(0,nrow=length(sl[[1]]),ncol=length(sl[[1]]));
            for (i in 1:length(sl[[1]]))
              for (j in i:length(sl[[1]]))
                {
                  if (i!=j)
                    {
                      vi<-strsplit(sl$state[[i]],NULL)[[1]]
                      vj<-strsplit(sl$state[[j]],NULL)[[1]]
                      rmat[j,i]<-length(vi)-sum(vi==vj);
                      rmat[i,j]<-rmat[j,i];
                    }
                }
            list(sl[[1]],rmat);
          }
  }
#
# produce a table of mismatches for a particular locus
#
landscape.mismatchdist<-function(lnum=1,Rland)
  {
    bd<-basediff(lnum,Rland);
    sl<-bd[[1]];
    dmat<-bd[[2]];
    lt<-landscape.locus(lnum,Rland);
    itbl<-table(lt[,(landscape.democol()+1):ncol(lt)]);
    ttbl<-as.table(table(c(0,seq(max(dmat))))*0);
    for (n in names(itbl))
      {
#        print(paste("Working on: ",n))
        mtbl<-as.table(table(dmat[seq(along=sl)[sl==as.numeric(n)],])*itbl[[n]]);
        for (cn in names(mtbl))
          {
            ttbl[[cn]]<-ttbl[[cn]]+mtbl[[cn]];
          }
      }
    ttbl
  }

#mismatch.pop <- function(Rland,pop=c(1:Rland$intparam$habitats))
#  {
#    maxdist <- 0
#    totrows <- 0
#    totcol <- 4
#    poplst <- NULL
#    listcnt <- 1
#    retdf <- NULL
#    for (i in pop)
#      {
#        popland <- Rland
#        popland$individuals <- Rland$individuals[landscape.populations(Rland)==i,]
#        for (j in 1:length(popland$loci))
#          {
#            tbldf <- as.data.frame(mismatchdist(lnum=j,popland))
#            tbldf$locus <- rep(j,nrow(tbldf))
#            tbldf$population <- rep(i,nrow(tbldf))
#            retdf <- rbind(retdf,tbldf)
#          }
#      }
#   names(retdf)[1] <- c("ntdiff")
#
#    retdf
#  }

#nucdiversity <- function(lnum=1, Rland)
#  {
#    if ((is.landscape(Rland))&&(Rland$loci[[lnum]]$type=253)) #is this a sequence from a valid landscape?
#        {
#          tbl <- mismatchdist(lnum,Rland)
#          quants <- as.numeric(names(tbl))
#          props <- tbl/sum(tbl)
#          sum(quants*props)/nchar(Rland$loci[[lnum]]$alleles[[1]]$state)
#        }
#    else
#      {
#        print("must pass a sequence type from a valid landscape")
#        NULL
#      }
#  }

#segsites <- function(lnum=1,Rland)
#  {
#    if ((is.landscape(Rland))&&(Rland$loci[[lnum]]$type=253)) #is this a sequence from a valid landscape?
#        {
#          m <- matrix(unlist(lapply(l$loci[[1]]$alleles,function(x) {strsplit(x$state,split='')[[1]]})),
#                      ncol=nchar(Rland$loci[[lnum]]$alleles[[1]]$state), byrow=TRUE)#
#
#          length(which(apply(m,2,function(x) length(unique(x))>1)))
#        }
#    else
#      {
#        print("must pass a sequence type from a valid landscape")
#        NULL
#      }
#  }
