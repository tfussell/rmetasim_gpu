#
#
#
#routines that operate on loci, both the actual states and the indices stored in
#individual genotypes

#returns a individual x ploidy matrix of aindices
landscape.locus <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          Rland$individuals[,c(rep(TRUE,landscape.democol()),landscape.locusvec(Rland)==lnum)]
        }
  }

#returns a individual x ploidy matrix of states
landscape.states <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          lmat <- as.data.frame(Rland$individuals[,c(rep(TRUE,landscape.democol()),landscape.locusvec(Rland)==lnum)])
          st <- landscape.locus.states(lnum,Rland)
          lmat[,landscape.democol()+1] <- st$state[sapply(lmat[,landscape.democol()+1],function(x,aindex){which(aindex==x)},aindex=st$aindex)]
          if (landscape.ploidy(Rland)[lnum]==2)
            {
              lmat[,landscape.democol()+2] <- st$state[sapply(lmat[,landscape.democol()+2],function(x,aindex){which(aindex==x)},aindex=st$aindex)]
            }
          lmat
        }
  }

#returns a vector of ploidys for all loci
landscape.ploidy<- function(Rland)
  {
    ploidy<-c();
    for (i in 1:Rland$intparam$locusnum)
      {
        ploidy<-c(ploidy,Rland$loci[[i]]$ploidy);
      }
    ploidy
  }

#returns a vector of locus ids
landscape.locusvec<- function(Rland)
  {
    p<-landscape.ploidy(Rland);
    lv<-c();
    for (i in  1:Rland$intparam$locusnum)
      {
        lv<-c(lv,rep(i,p[i]));
      }
    lv
  }
#
#takes a locus and returns the states and their indices
#
landscape.locus.states<-function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          ain<-c();
          sta<-c();
          locin <- landscape.locus(lnum,Rland)[,c(-1:-landscape.democol())]
#          print(locin)
          ainds <- unique(c(locin))
#          print(ainds)
          for (i in 1:length(Rland$loci[[lnum]]$alleles))
            {
              if (Rland$loci[[lnum]]$alleles[[i]]$aindex %in% ainds)
                {
                  ain<-c(ain,Rland$loci[[lnum]]$alleles[[i]]$aindex);
                  sta<-c(sta,Rland$loci[[lnum]]$alleles[[i]]$state);
                }
            }
          list(aindex=ain,state=sta);
        }
  }

#indxfreq <-function(lnum=1,Rland)
#  {
#    lv<-landscape.locus(lnum,Rland)[,c(rep(FALSE,landscape.democol()),
#                                           rep(TRUE,(ncol(locus(lnum,Rland))-landscape.democol())))];
#    if (landscape.ploidy(Rland)[lnum]==1)
#      {
#        table(landscape.populations(Rland),lv)
#      }
#    else
#      {
#        lv2<-c(lv[,1],lv[,2]);
#        table(rep(landscape.populations(Rland),2),lv2)
#      }
#    
#  }





