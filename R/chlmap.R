#'@name chlmap
#'@title mapping chl data using existing surfaces and chl coeff
#'@param yearmon numeric date in yyyymm format
#'@author Joseph Stachelek
#'@examples 
#'chlmap(yearmon=201502)

chlmap<-function(yearmon,fdir=getOption("fdir")){
  #yearmon<-201502
  
  params<-c("chlaiv","chla")
  #find coefficients that match yearmon####
  coeflist<-read.csv(file.path(fdir,"DF_GrabSamples","extractChlcoef2.csv"),header=T,na.strings="NA")[,-1]
  names(coeflist)<-tolower(names(coeflist))
  model<-coeflist[coeflist$yearmon==yearmon,"model"]
  coeflist<-coeflist[coeflist$yearmon==yearmon,4:15]
  namelist<-names(coeflist)[which(!is.na(coeflist))]
  coeflist<-coeflist[!is.na(coeflist)]
  
  #match raster surfaces to non-NA coefficients####
  namesalias<-read.table(text="
                       chlorophyll.a c6chl
                         c6chla c6chl
                         ")
  
  dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
  rlist<-list.files(dirlist[substring(basename(dirlist),1,6)==as.character(yearmon)],full.names=T,include.dirs=T,pattern="\\.grd$")
  plist<-tolower(sub("[.][^.]*$","",basename(rlist)))
  
  for(n in 1:length(plist)){
    if(any(plist[n]==namesalias[,1])){
      plist[n]<-as.character(namesalias[which(plist[n]==namesalias[,1]),2])
      #print(names(dt)[n])
    }
  }
  
  rlist<-rlist[which(!is.na(match(plist,params)))]
  plist<-plist[which(!is.na(match(plist,params)))]
  
  
  
  
  
}