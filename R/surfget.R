#'@name surfget
#'@title Retrieve an interpolated Dataflow surface
#'@description Retrieve an interpolated Dataflow surface from data directory
#'@export
#'@importFrom raster stack
#'@param rnge list of length one or two specifying date range in yyyymm format
#'@param param character parameter name
#'@param fdir character file path to local data directory
#'@examples \dontrun{ 
#'surfs <- surfget(rnge = c(201402, 201410), param = "chlext")
#'surfs <- surfget(rnge = 201509, param = "sal")
#'}

surfget <- function(rnge, param, fdir = getOption("fdir")){
  if(length(rnge) == 1){
    rnge <- c(rnge, rnge)
  }
  
  basepath <- file.path(fdir, "DF_Surfaces")
  
  flist <- list.files(basepath, include.dirs = T, full.names = T)
  flist <- suppressWarnings(flist[
      which(as.numeric(substring(basename(flist), 1, 6)) >= rnge[1])
    ])
  flist <- suppressWarnings(flist[
      which(as.numeric(substring(basename(flist), 1, 6)) <= rnge[2])
    ])
  flist <- list.files(flist, paste0("^", param, ".tif$"), full.names = TRUE, include.dirs = TRUE)
  
  suppressWarnings(raster::stack(flist))
}