#'@name surfplot
#'@title Plotting Interpolated Surfaces
#'@description accounts for corner case parameter spellings, variable specific contour breaks are (need to be) defined
#'@author Joseph Stachelek
#'@param rnge numeric string of no more than two dates in yyyymm format
#'@param params character. string of parameter fields to plot
#'@return output plots to plot window
#'@import rasterVis
#'@import grid
#'@import rgdal
#'@export
#'@examples surfplot(rnge=c(201502),params=c("sal"))

surfplot<-function(rnge=c(201402,201404),params=c("c6chl","sal")){
  
  library(rasterVis)
  library(grid)
  suppressMessages(library(rgdal))
  #library(quickmapr)
        
  fbcoast<-readOGR(dsn=file.path(fdir,"DF_Basefile/FBcoast_dissolve.shp"),layer="FBcoast_dissolve",verbose=FALSE)
  if(length(rnge)==1){
    rnge<-c(rnge,rnge)
  }
  namesalias<-read.table(text="
                       chlorophyll.a c6chl
                       c6chla c6chl
                       ")
  #define breaks
  brks<-read.table(text="
        sal list(seq(0,40,2))
        c6chl list(seq(0,300,50))
        ")
  
  dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
  
  minrnge<-min(which(substring(basename(dirlist),1,6)>=rnge[1]))
  maxrnge<-max(which(substring(basename(dirlist),1,6)<=rnge[2]))
  rlist<-list.files(dirlist[minrnge:maxrnge],full.names=T,include.dirs=T,pattern="\\.grd$")
  plist<-tolower(sub("[.][^.]*$","",basename(rlist)))
  
  for(n in 1:length(plist)){
    if(any(plist[n]==namesalias[,1])){
      plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
      #print(names(dt)[n])
    }
  }
  
  rlist<-rlist[which(!is.na(match(plist,params)))]
  plist<-plist[which(!is.na(match(plist,params)))]
  
  for(i in 1:length(rlist)){
    my.at<-unlist(eval(parse(text=as.character(brks[which(plist[i]==brks[,1]),2]))))
    
    print(levelplot(raster(rlist[i]),ylim=c(2772256,2798000),xlim=c(518000.2,566000),par.settings=PuOrTheme(),at=my.at,margin=FALSE,auto.key=FALSE,scales=list(draw=FALSE),main=paste(as.character(plist[i]),unlist(strsplit(rlist[i],"/"))[length(unlist(strsplit(rlist[i],"/")))-1]))+layer({SpatialPolygonsRescale(layout.north.arrow(),offset=c(563000,2775000),scale=4400)})+ layer(sp.polygons(readOGR(dsn=file.path(fdir,"DF_Basefile/FBcoast_dissolve.shp"),layer="FBcoast_dissolve",verbose=FALSE),fill="green",alpha=0.6)))
    }
}

#'@name avmap
#'@title create a difference map compared to average
#'@param yearmon survey of interest to compare against average
#'@param params variable name
#'@param tofile logical save output to disk?
#'@param percentcov numeric account for the different raster extents by setting the percent of all surveys required before a pixel is included in difference from average computations
#'@param tolerance numeric number of monts on either side of yearmon to include in the set of surfaces averaged . defaults to 1. 
#'@description takes a survey date as input and searches the DF_Surfaces folder for maps of the same parameter within the range of 1-2 months from yearmon for each year. These surfaces are averaged and compared to the surface from yearmon. 
#'@export
#'@examples avmap(yearmon=201505,params="sal",tofile=FALSE,percentcov=0.6,tolerance=1)
avmap<-function(yearmon=201505,params="sal",tofile=TRUE,percentcov=0.6,tolerance=1){
  
  flist.full<-list.files(file.path(fdir,"DF_Surfaces"),pattern="*.grd",recursive=T,include.dirs=T,full.names=T)
  flist<-flist.full[basename(flist.full)==paste(toupper(params),".grd",sep="")|basename(flist.full)==paste(tolower(params),".grd",sep="")]
  
  sdates<-data.frame(matrix(unlist(strsplit(dirname(flist),"/")),nrow=length(flist),byrow=T))
  sdates<-substring(sdates[,ncol(sdates)],1,6)
  
  cursurf<-raster(flist[which(sdates==yearmon)])
  flist<-flist[-which(sdates==yearmon)]
  sdates<-sdates[-which(sdates==yearmon)]
    
  curmon<-as.numeric(substring(yearmon,5,6))
  saveflist<-flist
  flist<-flist[as.numeric(substring(sdates,5,6))<=curmon+tolerance&as.numeric(substring(sdates,5,6))>=curmon-tolerance]
  remlist<-which(is.na(flist))
  if(length(remlist)>0){
  flist<-flist[-remlist]
  }
  
  sdates<-data.frame(matrix(unlist(strsplit(dirname(flist),"/")),nrow=length(flist),byrow=T))
  sdates<-substring(sdates[,ncol(sdates)],1,6)
  
  rstack<-stack(flist)
  rstack<-reclassify(rstack,c(-Inf,0,NA))
  rmean<-calc(rstack,fun=mean,na.rm=T)
  rlen<-sum(!is.na(rstack))
  #browser()
  rmean[rlen<(percentcov*length(flist))]<-NA
  
  plot(rmean,main=paste("Average",params,sdates[1],"-",sdates[length(sdates)],sep=" "))
  plot(cursurf-rmean,main="Difference from Average")
  
  if(tofile==TRUE){
  writeRaster(rmean,"meansurf.tif",format="GTiff",overwrite=T)
  writeRaster(cursurf,"cursurf.tif",format="GTiff",overwrite=T)
  writeRaster(cursurf-rmean,"diffsurf.tif",format="GTiff",overwrite=T)
  }
}

#'@name qmap
#'@title Publication quality maps with QGIS
#'@description not finished yet
#'@author Joseph Stachelek
#'@param rnge numeric string of no more than two dates in yyyymm format
#'@param params character. string of parameter fields to plot
#'@return output plots to the QGIS_plotting folder
#'@export
#'@examples qmap(rnge=c(201502),params=c("sal"))

qmap<-function(rnge=c(201402,201404),params=c("c6chl","sal")){
  
  #detect operating system####
  as.character(Sys.info()["sysname"])
  
  #dynamically update python script with below info?####
  
  if(length(rnge)==1){
    rnge<-c(rnge,rnge)
  }
  namesalias<-read.table(text="
                         chlorophyll.a c6chl
                         c6chla c6chl
                         ")
    dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
    
    minrnge<-min(which(substring(basename(dirlist),1,6)>=rnge[1]))
    maxrnge<-max(which(substring(basename(dirlist),1,6)<=rnge[2]))
    rlist<-list.files(dirlist[minrnge:maxrnge],full.names=T,include.dirs=T,pattern="\\.grd$")
    plist<-tolower(sub("[.][^.]*$","",basename(rlist)))
    
    for(n in 1:length(plist)){
      if(any(plist[n]==namesalias[,1])){
        plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
        #print(names(dt)[n])
      }
    }
    
    rlist<-rlist[which(!is.na(match(plist,params)))]
    plist<-plist[which(!is.na(match(plist,params)))]
    
    for(i in 1:length(rlist)){
    
  }
  
  #run shell script####
  shcmd<-matrix(read.delim(file.path(fdir,"QGIS_plotting","DFlow_QGIS_plotting_shellscript"),header=FALSE,stringsAsFactors=FALSE))[[1]]
  
  system2(shcmd)
        
}


# #create florida inset
# library(mapdata)
# data(worldHiresMapEnv)
# map("worldHires","usa",xlim=c(-86,-80),ylim=c(24,31),fill=TRUE,col="gray95")
