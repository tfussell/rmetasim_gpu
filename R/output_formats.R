#write.landscape.GDA <- function(rland, numi=24, fn = "landscape.GDA")
#  {
#    if (is.landscape(rland))
#      .Call("writeGDA",fn,rland,numi,PACKAGE = "rmetasim")
#  }


landscape.write.foreign <- function(rland, numi=24, fn = "foreign", fmt="GDA")
{
  if  (is.landscape(rland))
    {
      if (is.null(numi)) {numi <- (-1)}
      if (fmt %in% c("GDA","Gda","gda"))
        {
          .Call("writeGDA",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("Arlequin","arlequin","ArlequinDip","Arlequindip","arlequindip"))
        {
          .Call("writeArlequinDip",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("ArlequinHap","Arlequinhap","arlequinhap"))
        {
          .Call("writeArlequinHap",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("BIOSYS-1","BIOSYS","Biosys","biosys","biosys-1"))
        {
          .Call("writeBIOSYS",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("GenPop","Genpop","genpop"))
        {
          .Call("writeGenPop",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("R","r"))
        {
          .Call("writeR",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("Migrate","MigrateDiploid","migrate","migratediploid"))
        {
          .Call("writeMigrateDip",fn,rland,numi,PACKAGE = "rmetasim")
        }
      if (fmt %in% c("ReRat","rerat"))
        {
          .Call("writeGenPop",fn,rland,numi,PACKAGE = "rmetasim")
        }
    }
  else
    {
      print ("rland not in landscape format")
    }
}
