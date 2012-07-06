landscape.read <- function(fn = "filename")
  {
    if (file.exists(fn))
      {
        .Call("read_landscape",fn,PACKAGE = "rmetasim")
      }
    else
      {
        print (paste("Filename: ",fn,"does not appear to exist"))
        NULL
      }
  }


landscape.write <- function(rland, fn = "filename")
  {
    if (is.landscape(rland))
      .Call("write_landscape",fn,rland,PACKAGE = "rmetasim")
  }

landscape.clean <- function(rland)
  {
    if (is.landscape(rland))
      .Call("clean_landscape",rland,PACKAGE = "rmetasim")
  }


