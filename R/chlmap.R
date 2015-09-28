#'@name chlmap
#'@title Create a chlorophyll concentration surface using streaming data and regression against extracted chlorophyll
#'@details
#'A new interpolation is run after calculating an extracted chlorophyll for all streaming observations. Calculated values that exceed the maximum observed grab sample concentration are discarded.
#'@param yearmon numeric date in yyyymm format
#'@param remove.flags logical
#'@return An extracted chlorophyll surface and an updated FullDataSet file
#'@author Joseph Stachelek
#'@export
#'@examples 
#'\dontrun{
#'res<-chlmap(yearmon=201308)
#'}

chlmap<-function(yearmon, remove.flags = TRUE, fdir=getOption("fdir")){
  #library(DataflowR)
  #fdir<-getOption("fdir")
  #yearmon<-201311
  
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
  namelist_sq<-namelist[grep("2",namelist)]
  
  namesalias<-read.table(text="
                       c6chla c6chl
                       chla chlaiv 
                         ")
  
  for(n in 1:length(namelist)){
    if(any(namelist[n]==namesalias[,2])){
      namelist[n]<-as.character(namesalias[which(namelist[n]==namesalias[,2]),1])
    }
  }
  
  if(length(namelist_sq)>0){
  namelist_sq<-sapply(namelist_sq,function(x) substring(x,1,(nchar(x)-1)))
  for(n in 1:length(namelist_sq)){
    if(any(namelist_sq[n]==namesalias[,2])){
      namelist_sq[n]<-as.character(namesalias[which(namelist_sq[n]==namesalias[,2]),1])
    }
  }
  }
  
  #append a chlext column to the cleaned streaming data and interpolate
  dt<-streamget(yearmon,qa=TRUE)
  chlext<-dt[,namelist[1]]*coeflist[1]
  
  if(length(namelist_sq)>0){
    coeflist_sq<-coeflist[grep("2",namelist)]
    coeflist<-coeflist[-1*grep("2",namelist)]
    namelist<-namelist[-1*grep("2",namelist)]
    for(i in 1:(length(namelist_sq))){
      chlext<-chlext+(dt[,namelist_sq[i]]*coeflist[i])
    }
  }
  
  if(length(namelist)>2){
    for(i in 2:(length(namelist)-1)){
      chlext<-chlext+(dt[,namelist[i]]*coeflist[i])
    }
  }
  chlext<-chlext + coeflist[length(coeflist)]
  
  chlext[chlext<0]<-0
  chlext[chlext>range(grabget(yearmon, remove.flags = remove.flags)$chla,na.rm = TRUE)[2]]<-NA
  dt$chlext<-chlext
  
  
  if(file.exists(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep="")))){
  qafile<-read.csv(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep="")))
  qafile$chlext<-NA
  write.csv(qafile,file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep="")))
  }
  
  dtname<-file.path(fdir,.Platform$file.sep,"DF_FullDataSets",.Platform$file.sep,yearmon,"j.csv",fsep="")
  write.csv(dt,dtname,row.names = FALSE)
  
  if(file.exists(file.path(fdir,paste0("/DF_Subsets/chlext",yearmon,".csv"),fsep=""))){
    file.remove(file.path(fdir,paste0("/DF_Subsets/chlext",yearmon,".csv"),fsep=""))
    file.remove(file.path(fdir,paste0("/DF_Validation/chlext",yearmon,".csv"),fsep=""))
  }
  
  streaminterp(dt,paramlist = "chlext",yearmon = yearmon,tname = file.path(fdir,paste0("/DF_Subsets/chlext",yearmon,".csv"),fsep="") ,vname = file.path(fdir,paste0("/DF_Validation/chlext",yearmon,".csv"),fsep=""), missprop = (1/3))
  
}  
  #####old code from when raster algebra was used#####
  #match raster surfaces to non-NA coefficients####
#   dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
#   rlist<-list.files(dirlist[substring(basename(dirlist),1,6)==as.character(yearmon)],full.names=T,include.dirs=T,pattern="\\.tif$")
#   plist<-tolower(sub("[.][^.]*$","",basename(rlist)))
#   
#   namesalias<-read.table(text="
#                        chlorophyll.a c6chl
#                        c6chla c6chl
#                        chla chlaiv 
#                          ")
#   
#   for(n in 1:length(plist)){
#     if(any(plist[n]==namesalias[,1])){
#       plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
#     }
#   }
#   
#   pr_order<-match(c(namelist),plist)[!is.na(match(c(namelist),plist))]
#   plist<-plist[pr_order]
#   rlist<-rlist[pr_order]
  
  #create rlist2 if polynomial coefficients
#   if(length(grep("2",namelist))>0){
#     poly=TRUE
#     rlist2<-rlist[match(sapply(namelist[grep("2",namelist)],function(x) substring(x,1,(nchar(x)-1))),plist)]
#   }
#   
#   if(length(rlist)==0){
#     stop("no surfaces matching parameter names")
#   }
#   
#   res<-raster::raster(rlist[1])*coeflist[1]
#   if(length(rlist)>1){
#   for(i in 2:length(rlist)){
#     res<-res+(raster::raster(rlist[i])*coeflist[i])
#   }
#   }
#   
#   if(poly==TRUE){
#     
#   }
#   
#   res<-res + coeflist[length(coeflist)]
#   res<-raster::reclassify(res,c(-Inf,0,0))
#   #sp::plot(res)
#   
#   if(tofile==TRUE){
#     raster::writeRaster(res,filename=file.path(fdir,"DF_Surfaces",yearmon,"chlext.tif"),overwrite=TRUE,format = "GTiff")
#   }
#   res

