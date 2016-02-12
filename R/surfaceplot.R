#'@name surfplot
#'@title Plotting Interpolated Surfaces
#'@description accounts for corner case parameter spellings, variable specific contour breaks are (need to be) defined
#'@author Joseph Stachelek
#'@param rnge numeric string of no more than two dates in yyyymm format
#'@param params character. string of parameter fields to plot
#'@param fdir character file path to local data directory
#'@param yext numeric length 2 y extent
#'@param xext numeric length 2 x extent
#'@return output plots to plot window
#'@importFrom rgdal readOGR
#'@importFrom rasterVis levelplot
#'@importFrom latticeExtra layer
#'@importFrom sp spplot
#'@importFrom raster raster stack reclassify calc writeRaster
#'@export
#'@examples \dontrun{
#'surfplot(rnge=c(200707), params = c("cdom"))
#'surfplot(201513, c("ph", "c6turbidity", "c6chl", "c6cdom"),
#' yext = c(2786102, 2797996), xext = c(557217, 567415))
#'}

surfplot <- function(rnge = c(201402, 201404), params = c("c6chl", "sal"), fdir = getOption("fdir"), yext = c(2772256, 2798000), xext = c(518000.2, 566000)){
  print(fdir)

  fbcoast <- rgdal::readOGR(dsn = file.path(fdir, "DF_Basefile", "FBcoast_big.shp"), layer = "FBcoast_big", verbose = TRUE)
  
  if(length(rnge) == 1){
    rnge <- c(rnge, rnge)
  }
  namesalias <- read.table(text="
                       chlorophyll.a c6chl
                       c6chla c6chl
                       ")
  #define breaks
  brks <- read.table(text = "
        sal list(seq(0,40,2))
        c6chl list(seq(50,200,10))
        chlext list(seq(0,5,0.5),seq(10,30,5))
        chlext_hi list(seq(0,5,0.5),seq(10,30,5))
        chlext_low list(seq(0,5,0.5),seq(10,30,5))
        temp list(seq(14,36,2))
        c6temp list(seq(14,36,2))
        c6cdom list(seq(80,360,40))
        ph list(seq(7.3,8.1,0.1))
        c6turbidity list(seq(0,45,5))")
  
  dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
  
  minrnge<-min(which(substring(basename(dirlist),1,6)>=rnge[1]))
  maxrnge<-max(which(substring(basename(dirlist),1,6)<=rnge[2]))
  rlist<-list.files(dirlist[minrnge:maxrnge],full.names=T,include.dirs=T,pattern="\\.tif$")
  plist <- tolower(sub("[.][^.]*$","",basename(rlist)))
  
  for(n in 1:length(plist)){
    if(any(plist[n]==namesalias[,1])){
      plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
      #print(names(dt)[n])
    }
  }
  
  rlist<-rlist[which(!is.na(match(plist,params)))]
  plist<-plist[which(!is.na(match(plist,params)))]
  
  for(i in 1:length(rlist)){
    #i<-1
    
    my.at <- unlist(eval(parse(text = as.character(brks[which(plist[i] == brks[,1]), 2]))))
      
    print(rasterVis::levelplot(raster::raster(rlist[i]),ylim =yext, xlim = xext,par.settings=rasterVis::PuOrTheme(),at=my.at,margin=FALSE,auto.key=FALSE,scales=list(draw=FALSE),main=paste(as.character(plist[i]),unlist(strsplit(rlist[i],"/"))[length(unlist(strsplit(rlist[i],"/")))-1]))+latticeExtra::layer({sp::SpatialPolygonsRescale(sp::layout.north.arrow(),offset=c(563000,2775000),scale=4400)})+ latticeExtra::layer(sp::sp.polygons(rgdal::readOGR(dsn = file.path(fdir, "DF_Basefile", "FBcoast_big.shp"), layer = "FBcoast_big", verbose = TRUE), fill="green",alpha=0.6)))
    
    }
}

#'@name avmap
#'@title create a difference map compared to average
#'@param yearmon survey of interest to compare against average
#'@param params variable name
#'@param tofile logical save output to disk?
#'@param percentcov numeric account for the different raster extents by setting the percent of all surveys required before a pixel is included in difference from average computations
#'@param tolerance numeric number of monts on either side of yearmon to include in the set of surfaces averaged . defaults to 1.
#'@param fdir character file path to local data directory 
#'@description takes a survey date as input and searches the DF_Surfaces folder for maps of the same parameter within the range of 1-2 months from yearmon for each year. These surfaces are averaged and compared to the surface from yearmon. 
#'@export
#'@importFrom raster raster stack reclassify calc writeRaster
#'@examples \dontrun{
#'avmap(yearmon = 201505, params = "sal", tofile = FALSE, percentcov = 0.6, tolerance = 1, fdir = fdir)
#'avmap(yearmon = 201505, params = "sal", tofile = FALSE)
#'}

avmap <- function(yearmon = 201505, params = "sal", tofile = TRUE, percentcov = 0.6, tolerance = 1, fdir = getOption("fdir")){

#generate-flist=========================================================#
  flist.full <- list.files(file.path(fdir, "DF_Surfaces"), pattern = "*.grd", recursive = T, include.dirs = T, full.names = T)
  flist <- flist.full[basename(flist.full) == paste(toupper(params), ".grd", sep = "") | basename(flist.full) == paste(tolower(params), ".grd", sep = "")]
  
  sdates <- data.frame(matrix(unlist(strsplit(dirname(flist), "/")), nrow = length(flist), byrow = T))
  sdates <- substring(sdates[,ncol(sdates)], 1, 6)
  
  cursurf <- raster::raster(flist[which(sdates == yearmon)])
  flist <- flist[-which(sdates == yearmon)]
  sdates <- sdates[-which(sdates == yearmon)]
    
  curmon <-as.numeric(substring(yearmon, 5, 6))
  
  flist <- flist[as.numeric(substring(sdates, 5, 6)) <= curmon + tolerance & as.numeric(substring(sdates, 5, 6)) >= curmon - tolerance]
  remlist <- which(is.na(flist))
  if(length(remlist) > 0){
    flist <- flist[-remlist]
  }

#stack-and-calculate====================================================#
  rstack <- raster::stack(flist)
  rstack <- raster::reclassify(rstack, c(-Inf, 0, NA))
  rmean <- raster::calc(rstack, fun = mean, na.rm = T)
  rlen <- sum(!is.na(rstack))
  
  rmean[rlen < (percentcov * length(flist))] <- NA
  res <- cursurf - rmean
  
#plotting===============================================================#
  sdates <- data.frame(matrix(unlist(strsplit(dirname(flist), "/")), nrow = length(flist), byrow = T))
  sdates <- substring(sdates[,ncol(sdates)], 1, 6)
  
  sp::plot(rmean, main = paste("Average", params, sdates[1], "-", sdates[length(sdates)], sep = " "))
  sp::plot(cursurf - rmean, main = "Difference from Average")

#save-to-file===========================================================#
  if(tofile == TRUE){
    #raster::writeRaster(rmean, "meansurf.tif", format = "GTiff", overwrite = T)
    #raster::writeRaster(cursurf, "cursurf.tif", format = "GTiff", overwrite = T)
    raster::writeRaster((cursurf - rmean), file.path(fdir, "DF_Surfaces", yearmon, paste0("diff", params, ".tif")), format = "GTiff", overwrite = T)
  }
  
  res
}

#'@name grassmap
#'@title Publication quality maps with GRASS
#'@description not finished yet
#'@author Joseph Stachelek
#'@param fpath file.path to geotiff file
#'@param rnge numeric string of 1, 2, or more dates dates in yyyymm format. A length 1 rnge will produce a single plot, a length 2 rnge will produce a series of plots bookended by the two dates, a rnge object with more than 2 dates will produce a series of plots exactly corresponding to the dates provided.
#'@param params character vector of parameter fields to plot legends and color ramps are defined for sal, chlext, and diffsal
#'@param fdir character file path to local data directory
#'@param cleanup logical remove intermediate rasters and shapefiles?
#'@param rotated logical rotate canvas to fit Florida Bay more squarely? This requires the i.rotate extension to be installed and addons configured (not working).
#'@param labelling logical lablel output with yearmon?
#'@param numrow numeric number of rows in a multipanel
#'@param numcol numeric number of columns in a multipanel
#'@param mapextent numeric vector of length 4
#'@param basin character basin name
#'@param print_track logical print dataflow track?
#'@return output plots to the QGIS_plotting folder
#'@details Probably need to implement this as a seperate package. Set param to "diffsal" to plot outpot of avmap function. Will output an imagemagick plot to the working directory and a pdf plot to the file.path(getOption("fdir"), "QGIS_Plotting") folder. The optional fpath argument only supports pointing to a single geotiff.
#'@import rgrass7
#'@import maptools
#'@import rgeos
#'@export
#'@examples \dontrun{
#'grassmap(rnge = c(201512), params = c("sal"), basin="Manatee Bay")
#'grassmap(rnge = c(201512), params = c("sal"))
#'grassmap(rnge = c(201512), params = c("diffsal"))
#'grassmap(rnge=c(201407),params=c("chlext"))
#'grassmap(rnge=c(201507),params=c("sal"), mapextent = c(494952.6, 564517.2, 2758908, 2799640))
#'grassmap(rnge=c(201509),params=c("chlext"), mapextent = c(494952.6, 564517.2, 2758908, 2799640))
#'grassmap(rnge=c(201512), params = "diffsal", mapextent = c(494952.6, 564517.2, 2758908, 2799640))
#'grassmap(201513, "chlext", mapextent = c(557217, 567415, 2786102, 2797996), print_track = TRUE)
#'grassmap(fpath = file.path(getOption("fdir"), "DF_Surfaces", "200904", "sal.tif"), params = "sal")
#'
#'#create a new color ramp by editing DF_Basefile/*.file and update figure makefile
#'logramp(n = 9, maxrange = 20) #chlext
#'scales::show_col(viridis::viridis_pal()(9))
#'}

grassmap <- function(fpath = NULL, rnge = NULL , params, mapextent = NA, numrow = NULL, numcol = NULL, fdir = getOption("fdir"), basin = "full", labelling = TRUE, print_track = FALSE, cleanup = TRUE, rotated = TRUE){

  if(as.character(Sys.info()["sysname"]) != "Linux"){
    stop("This function only works with Linux!")
  }
  
  paramkey <- read.table(text = "
sal,salrules.file
chlext,chlextrules.file
chlext_low,chlextrules.file
chlext_hi,chlextrules.file
diffsal,diffsalrules.file",
  sep = ",", stringsAsFactors = FALSE)
  rulesfile <- paramkey[which(params == paramkey[,1]), 2]
  
  create_rlist <- function(rnge, params){
    if(length(rnge) == 1){
      rnge <- c(rnge, rnge)
    }
    
    namesalias <- read.table(text = "
                         chlorophyll.a c6chl
                         c6chla c6chl
                         ")
    
    dirlist <- list.dirs(file.path(fdir, "DF_Surfaces"), recursive = F)
    minrnge <- min(which(substring(basename(dirlist), 1, 6) >= rnge[1]))
    maxrnge <- max(which(substring(basename(dirlist), 1, 6) <= rnge[2]))
    rlist <- list.files(dirlist[minrnge:maxrnge], full.names = T, include.dirs = T, pattern = "\\.tif$")
    plist <- tolower(sub("[.][^.]*$", "", basename(rlist)))
    
    for(n in 1:length(plist)){
      if(any(plist[n] == namesalias[,1])){
        plist[n] <- as.character(namesalias[which(plist[n] == namesalias[,1]), 2])
      }
    }
    
    if(length(rnge) > 2){
      rnge <- rnge[order(rnge)]
      rlist <- list.files(dirlist[which(substring(basename(dirlist), 1, 6) %in% rnge)], full.names = T, include.dirs = T, pattern = "\\.tif$")
      plist <- tolower(sub("[.][^.]*$", "", basename(rlist)))
    }
    
    rlist <- rlist[which(!is.na(match(plist, params)))]
    plist <- plist[which(!is.na(match(plist, params)))]
    list(rlist = rlist, plist = plist)
  }
  
  if(length(fpath) == 0){
    rlist <- create_rlist(rnge = rnge, params = params)$rlist
    plist <- create_rlist(rnge = rnge, params = params)$plist
  }else{
    rlist <- fpath
    plist <- rep(params, each = length(rlist))
  }
  
  fathombasins <- rgdal::readOGR(file.path(fdir, "DF_Basefile/fathom_basins_proj.shp"), layer = "fathom_basins_proj", verbose = FALSE)
  fboutline <- rgdal::readOGR(dsn = file.path(getOption("fdir"), "DF_Basefile/FBcoast_big.shp"), layer = "FBcoast_big", verbose = FALSE)
  
  if(print_track == TRUE){
    surveytrack <- coordinatize(streamget(rnge[1]), latname = "lat_dd", lonname = "lon_dd")
  }
  #interp_pnts <- coordinatize(read.csv(file.path(fdir, "DF_Subsets", paste0(rnge[1], "s.csv"))), latname = "lat_dd", lonname = "lon_dd")
  
  print(rlist)
  
  for(i in 1:length(rlist)){
    if(basin != "full"){
      firstras <- raster::raster(rlist[1])
      firstras <- raster::crop(firstras, fathombasins[fathombasins$NAME == basin,])
    }else{
      firstras <- raster::raster(rlist[i])
    }
    
    if(basin != "full"){
      tempras <- raster::raster(rlist[i])
      tempras <- raster::crop(tempras, fathombasins[fathombasins$NAME == basin,])
    }else{
      tempras <- raster::raster(rlist[i])
    }
    
    rasname <- paste(substring(dirname(rlist[i]), nchar(dirname(rlist[i])) - 5, nchar(dirname(rlist[i]))))
    raspath <- file.path(paste(fdir, "/QGIS_plotting", sep = ""), paste(rasname, ".tif", sep = ""))
    outpath = file.path(paste(fdir, "/QGIS_plotting", sep = ""), paste(rasname, "poly.shp" ,sep = ""))
    
    raster::writeRaster(tempras, raspath, format = "GTiff", overwrite = TRUE)
    shellcmds = paste("gdal_polygonize.py", raspath, "-f","'ESRI Shapefile'", outpath) 
    system(shellcmds)
    outpoly <- rgdal::readOGR(dsn = outpath, layer = paste(rasname, "poly", sep = ""), verbose = TRUE)
    #requireNamespace("maptools")
    require("rgeos")
    require("maptools")#cannot seem to execute below without call to require
    maptools::gpclibPermit()
    outpoly <- maptools::unionSpatialPolygons(outpoly, IDs = rep(1, length(outpoly)))
    outlines <- as(outpoly, 'SpatialLines')
    outlines <- sp::SpatialLinesDataFrame(outlines, data = as.data.frame(1))
    
    #GRASS block===============================================================#
    loc <- rgrass7::initGRASS("/usr/lib/grass70", home = file.path(fdir, "QGIS_plotting"), override = TRUE)
    
    #raster
    firstras <- as(firstras, "SpatialGridDataFrame")
    firstras.g <- rgrass7::writeRAST(firstras, "firstras", flags = c("overwrite", "quiet"))
    
    tempras <- as(tempras, "SpatialGridDataFrame")
    tempras.g <- rgrass7::writeRAST(tempras, "tempras", flags = c("overwrite"))
    rgrass7::execGRASS("g.region", raster = "tempras")
    
     if(params %in% c("chlext", "chlext_low", "chlext_hi")){
       if(length(grep("_log", rulesfile)) != 0){
         rulesfile <- gsub("_log", "", rulesfile)
       }
     }
    
    rgrass7::execGRASS("r.colors", map = "tempras", rules = file.path(fdir, "DF_Basefile", rulesfile))
    
    rgrass7::execGRASS("r.grow", input = "tempras", output = "tempras2", radius = 1.3, flags = c("overwrite", "quiet"))
    
    #Florida Bay outline
    if(is.na(mapextent)){
      fboutline <- raster::crop(fboutline, raster::extent(raster::raster(tempras)))
    }else{
      fboutline <- raster::crop(fboutline, mapextent)
    }
    fbvec.g <- rgrass7::writeVECT(fboutline, "fbvec", v.in.ogr_flags = c("o"))
    #rgrass7::execGRASS("g.region",vector="fbvec")
    rgrass7::execGRASS("v.colors", map = "fbvec", column = "cat", color = "grey")
    
    #survey track
    if(print_track == TRUE){
      trackvec.g <- rgrass7::writeVECT(surveytrack, "trackvec", v.in.ogr_flags = c("o"))
      rgrass7::execGRASS("v.colors", use = "cat", map = "trackvec", color = "grey")
    }
    
    #raster outline
    outvec.g <- rgrass7::writeVECT(outlines, "outvec", v.in.ogr_flags = c("o"))
    rgrass7::execGRASS("g.region", vector = "outvec")
    rgrass7::execGRASS("g.region", raster = "firstras")
    rgrass7::execGRASS("g.region", vector = "fbvec")
    
    if(labelling == TRUE){
#     #compose plotting commands here####
    fileConn <- file(file.path(fdir, "QGIS_plotting", "grassplot.file"))
    writeLines(c("raster tempras2",
                 "vlines outvec",
                 "        color black",
                 "        style dashed",
                 "        end",
                 "vareas fbvec",
                 "        masked y",
                 "        end",
                 paste("text 17% 85% ",rasname,sep=""),
                 "        fontsize 35",
                 "        background white",
                 "        border black",
                 "        end",
                 "end"),fileConn)
    close(fileConn)
    if(print_track == TRUE){
      fileConn <- file(file.path(fdir, "QGIS_plotting", "grassplot.file"))
      writeLines(c("raster tempras2",
                   "vpoints trackvec",
                   "        color black",
                   "        fcolor black",
                   "        symbol basic/cross1",
                   "        size 5",
                   "        end",
                   "vlines outvec",
                   "        color black",
                   "        style dashed",
                   "        end",
                   "vareas fbvec",
                   "        masked y",
                   "        end",
                   paste("text 17% 85% ",rasname,sep=""),
                   "        fontsize 35",
                   "        background white",
                   "        border black",
                   "        end",
                   "end"),fileConn)
      close(fileConn)
    }
    }else{
      fileConn<-file(file.path(fdir,"QGIS_plotting","grassplot.file"))
      writeLines(c("raster tempras2",
                   "vlines outvec",
                   "        color black",
                   "        style dashed",
                   "        end",
                   "vareas fbvec",
                   "        masked y",
                   "        end",
                   "end"),fileConn)
      close(fileConn)
    }
    
    rgrass7::execGRASS("ps.map", input = file.path(paste(fdir, "/QGIS_plotting", sep=""), "grassplot.file"), output = file.path(paste(fdir, "/QGIS_plotting", sep = ""), paste(substring(dirname(rlist[i]), nchar(dirname(rlist[i])) - 5, nchar(dirname(rlist[i]))), ".pdf", sep = "")), flags = "overwrite")

#==================================================================#

    legendalias <- read.table(text = "chlext,Chlorophyll (ug/L)
chlext_low,Chlorophyll (ug/L)
chlext_hi,Chlorophyll (ug/L)
sal,Salinity
diffsal,Salinity minus average", sep = ",", stringsAsFactors = FALSE)
    
    legendname <- legendalias[which(params == legendalias[,1]), 2]
    
    if(params == "sal"){
      paramxcoord <- 2160
      legendunits<-seq(from = 5,to = 54,by = 0.1)
      legendunits_print <- "'5 10 15 20 25 30 35 40'"
      legendunits_spacing <- 220
      legend_xlim <- 270
      legend_crop_extent <- 2404
    }
    if(params %in% c("chlext", "chlext_low", "chlext_hi")){
      paramxcoord <- 1980
      legendunits <- log(seq(from = 0, to = 13.5, by = 0.1) + 1)
      legendunits_print <- "'0.0 0.7 2.0 4.0 7.0 13.0'"
      legendunits_spacing <- 275
      
      if(length(grep("_log", rulesfile)) == 0){
        rulesfile <- paste0(rulesfile,"_log")
      }
      
      legend_xlim <- 300
      legend_crop_extent <- 2404
      #legend_xlim <- 140
      #legend_crop_extent <- 2188
    }
    
    #browser()
    if(params == "diffsal"){
      paramxcoord <- 1960
      legendunits <- seq(from = -30, to = 35, by = 1)
      legendunits_print <- "'-30 -25 -20 -15 -10 -5 0 5 10 15 20'"
      legendunits_spacing <- 120
      legend_xlim <- 270
      legend_crop_extent <- 2404 
    }
    #print legend========================================================#
    
    legras <- raster::raster(tempras)
    legras[1:length(legras)] <- legendunits
    #legras[1:(length(legendunits) + 1)] <- legendunits #old implementation
    tempras <- as(legras, "SpatialGridDataFrame")
    tempras.g <- rgrass7::writeRAST(tempras, "tempras", flags = c("overwrite"))
    
    #browser()
    #hist(as.numeric(unlist(tempras@data)))
    #print(legendunits)
    #print(rulesfile)
    
    rgrass7::execGRASS("r.support", map = "tempras", units = legendname)
    rgrass7::execGRASS("g.region", raster = "tempras")
    rgrass7::execGRASS("r.colors", map = "tempras", rules = file.path(fdir, "DF_Basefile", rulesfile))
    
    rgrass7::execGRASS("r.to.vect", input = "tempras", output = "outvec", type = "area", flags = "overwrite")
    rgrass7::execGRASS("v.hull", input = "outvec", output = "outvec2", flags = "overwrite")
    
    rgrass7::execGRASS("ps.map", input = file.path(paste(fdir, "/QGIS_plotting", sep = ""), "legendplot.file"), output = file.path(paste(fdir, "/QGIS_plotting", sep = ""), "legend", paste("legend", ".pdf", sep = "")), flags = "overwrite")
    
    #==================================================================#
    #browser()
    system(paste(
      "echo", "'", legendname, substring(legendunits_print, 2, nchar(legendunits_print) - 1), legendunits_spacing, legend_xlim, legend_crop_extent, "'", ">> 'single.txt'"))
    
    if(length(rlist) == 1){
      makefile <- file.path(fdir, "DF_Basefile", "Makefile_single")
      system(paste0("make -f ", makefile, 
" testpanel.png BASEDIR=", fdir,
" YEARMON=", paste0(substring(dirname(rlist[i]), nchar(dirname(rlist[i])) - 5, nchar(dirname(rlist[i])))),
" PARAM=", shQuote(legendname),
" LEGENDUNITS=", legendunits_print,
" LEGENDUNITSSPACING=", legendunits_spacing,
" LEGEND_XLIM=", legend_xlim,
" PARAMXCOORD=", paramxcoord,
" LEGEND_CROP_EXTENT=", legend_crop_extent
))
      
      if(cleanup == TRUE){
        system(paste0("make -f ", makefile, " clean"))
      }
    }
    #==================================================================#
    
    if(cleanup == TRUE){
      rmlist <- list.files(file.path(paste(fdir, "/QGIS_plotting", sep = "")), pattern = paste(rasname, "*", sep = ""), include.dirs = TRUE, full.names = TRUE)
      rmlist <- rmlist[-grep("*.pdf", rmlist)]
      file.remove(rmlist)
    }    
  }
  
  #==================================================================#
  #browser()
  system(paste(
    "echo", "'", legendname, substring(legendunits_print, 2, nchar(legendunits_print) - 1), legendunits_spacing, legend_xlim, legend_crop_extent, "'", ">> 'multi.txt'"))
  
  #assumes that all pdfs in QGIS_plotting are to be part of panel
  if(length(rlist) > 1){ # & !is.na(panel.dim)
    makefile <- file.path(fdir, "DF_Basefile", "Makefile_multi")
    system(paste0("make -f ", makefile,
                  " multipanel.png BASEDIR=", fdir,
                  " PARAM=", shQuote(legendname),
                  " LEGENDUNITS=", legendunits_print,
                  " LEGENDUNITSSPACING=", legendunits_spacing,
                  " LEGEND_XLIM=", legend_xlim,
                  " PARAMXCOORD=", paramxcoord,
                  " NROW=", numrow,
                  " NCOL=", numcol,
                  " LEGEND_CROP_EXTENT=", legend_crop_extent
                  ))
    if(cleanup == TRUE){
      system(paste0("make -f ", makefile, " clean"))
    }
  }
  
  if(cleanup == TRUE){
    rmlist <- list.files(file.path(paste(fdir,"/QGIS_plotting", sep = "")), pattern = paste("out", "*", sep = ""), include.dirs = TRUE, full.names = TRUE)
    file.remove(rmlist)
  }
}

# #create florida inset
# library(mapdata)
# data(worldHiresMapEnv)
# map("worldHires","usa",xlim=c(-86,-80),ylim=c(24,31),fill=TRUE,col="gray95")
