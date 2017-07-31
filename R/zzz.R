.onAttach <- function(libname = find.package("DataflowR"), pkgname = "DataflowR"){
  
  localpath <- system.file("localpath", package = "DataflowR")
  fdir <- matrix(utils::read.delim(localpath, header = FALSE, 
                                   stringsAsFactors = FALSE))[[1]][1]
  
  packageStartupMessage(paste("DataflowR Data Directory:",
                              fdir, "\n", "To change, modify:", localpath))
                        
  options("fdir" = fdir)
}
