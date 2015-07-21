#'@name chlmap
#'@title mapping chl data using existing surfaces and chl coeff
#'@param yearmon numeric date in yyyymm format
#'@author Joseph Stachelek
#'@export
#'@examples 
#'res<-chlmap(yearmon=200808)

chlmap<-function(yearmon,tofile=FALSE,fdir=getOption("fdir")){
  #library(DataflowR)
  #fdir<-getOption("fdir")
  #yearmon<-201002
  
  params<-c("chlaiv","chla")
  #find coefficients that match yearmon####
  coeflist<-read.csv(file.path(fdir,"DF_GrabSamples","extractChlcoef2.csv"),header=T,na.strings="NA")[,-1]
  names(coeflist)<-tolower(names(coeflist))
  model<-coeflist[coeflist$yearmon==yearmon,"model"]
  
  if(length(model)==0){
    stop("No chl model fit for this survey")
  }
  
  coeflist<-coeflist[coeflist$yearmon==yearmon,4:16]
  namelist<-names(coeflist)[which(!is.na(coeflist))]
  coeflist<-coeflist[!is.na(coeflist)]
  
  
  
  #match raster surfaces to non-NA coefficients####
  dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
  rlist<-list.files(dirlist[substring(basename(dirlist),1,6)==as.character(yearmon)],full.names=T,include.dirs=T,pattern="\\.tif$")
  plist<-tolower(sub("[.][^.]*$","",basename(rlist)))
  
  namesalias<-read.table(text="
                       chlorophyll.a c6chl
                       c6chla c6chl
                       chla chlaiv 
                         ")
  
  for(n in 1:length(plist)){
    if(any(plist[n]==namesalias[,1])){
      plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
    }
  }
  
  pr_order<-match(c(namelist),plist)[!is.na(match(c(namelist),plist))]
  plist<-plist[pr_order]
  rlist<-rlist[pr_order]
  
  #create rlist2 if polynomial coefficients
  if(length(grep("2",namelist))>0){
    poly=TRUE
    rlist2<-rlist[match(sapply(namelist[grep("2",namelist)],function(x) substring(x,1,(nchar(x)-1))),plist)]
  }
  
  if(length(rlist)==0){
    stop("no surfaces matching parameter names")
  }
  
  res<-raster::raster(rlist[1])*coeflist[1]
  if(length(rlist)>1){
  for(i in 2:length(rlist)){
    res<-res+(raster::raster(rlist[i])*coeflist[i])
  }
  }
  
  if(poly==TRUE){
    
  }
  
  res<-res + coeflist[length(coeflist)]
  res<-raster::reclassify(res,c(-Inf,0,0))
  #sp::plot(res)
  
  if(tofile==TRUE){
    raster::writeRaster(res,filename=file.path(fdir,"DF_Surfaces",yearmon,"chlext.tif"),overwrite=TRUE,format = "GTiff")
  }
  res
}
