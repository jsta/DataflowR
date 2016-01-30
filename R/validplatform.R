#'@name validdbhydro
#'@title Validate Dataflow output against DBHYDRO data
#'@description Validate Dataflow output against DBHYDRO data
#'@param yearmon numeric 6 character date in YYYYMM format
#'@param params character string of parameter names
#'@param fdir file.path to data folder
#'@param tolerance numeric number of days pre/post streaming to search for in DBHYDRO
#'@param dbfile file.path pointing at a cleaned database of dbhydro measurements
#'@examples \dontrun{
#'dbfile<-file.path("/home","jose","Documents","Science","Data","ENP_MMN","WQSuM_2015.csv")
#'validdbhydro("201502",params="chlext",tolerance=10, dbfile=dbfile)
#'validdbhydro(201410,params="chlext",tolerance=13,dbfile=dbfile)
#'validdbhydro(201407,params="chlext",tolerance=2,dbfile=dbfile)
#'validdbhydro(201404,params="chlext",tolerance=17,dbfile=dbfile)
#'validdbhydro(201311,params="chlext",tolerance=10,dbfile=dbfile)
#'validdbhydro(201308,params="chlext",tolerance=5,dbfile=dbfile)
#'}

validdbhydro<-function(yearmon,params="chlext",tolerance=1,dbfile=file.path(),fdir=getOption("fdir")){
  
#   library(DataflowR)
#   dbfile<-file.path("/home","jose","Documents","Science","Data","ENP_MMN","WQSuM_2015.csv")
#   yearmon<-201212#200804#201106#
#   tolerance<-30#
#   fdir<-getOption("fdir")#
#   params<-"chlext"#
  
  flabsites <- c(paste0("FLAB",c(12,14,15,13,11,23,24,10,"08","09","05","06","01")),"BB51","BB56", "BB50")
  
  print(yearmon)
  
  streamdt<-streamget(yearmon,qa=TRUE)
  
  if(all(is.na(streamdt$datetime))){
    streamdt$datetime<-date456posix(streamdt$date,century="20")
  }
  
  streamdates<-as.POSIXct(unique(strftime(as.POSIXct(streamdt$datetime),format="%Y-%m-%d")))
  streamdates<-streamdates[!is.na(streamdates)]
  
  min_date<-strftime(min(streamdates)-((60*60*24)*tolerance),format="%Y-%m-%d")
  max_date<-strftime(max(streamdates)+((60*60*24)*tolerance),format="%Y-%m-%d")
    
  if(length(dbfile)==0){
    
    date_min <- as.character(date456posix(streamget(yearmon)[1,"date"], century = 20) - (60 * 60 * 24 * tolerance))
    date_max <- as.character(date456posix(streamget(yearmon)[1,"date"], century = 20) + (60 * 60 * 24 * tolerance))
    
  dbdt <- dbhydroR::getwq(flabsites, date_min = date_min, date_max = date_max, test_name = "CHLOROPHYLL-A(LC)")
    ##need to add additional cleaning commands here
  }else{
    dbdt<-read.csv(dbfile,stringsAsFactors = FALSE)
  }
  
  
#   comparedates<-unlist(lapply(as.POSIXct(unique(dbdt$collection.date)), function(x) abs(x-mean(as.POSIXct(streamdates)))))
#   
#   as.POSIXct(unique(dbdt$collection.date))[order(as.POSIXct(unique(dbdt$collection.date)))]
#   
#   ndates<-comparedates[order(comparedates)]
#   
#   hist(unlist(lapply(as.POSIXct(unique(dbdt$collection.date)), finddate)))
#   
#   hist(unlist(lapply(as.POSIXct(unique(dbdt$collection.date)), finddate)))
#   
#   dbdt[
#     hist(abs(as.numeric(as.POSIXct(unique(dbdt$collection.date))-as.POSIXct(streamdates))))
#     
#     <=as.POSIXct(max_date)&as.POSIXct(dbdt$collection.date)>=as.POSIXct(min_date),]
  
  
  dbdt<-dbdt[as.POSIXct(dbdt$collection.date)<=as.POSIXct(max_date)&as.POSIXct(dbdt$collection.date)>=as.POSIXct(min_date),]
  
  if(nrow(dbdt)==0){
    stop("No matching DBHYDRO data found within the specified tolerance interval")
  }
  
  dbdt<-coordinatize(dbdt,latname="latdec",lonname="londec")
  
  dirlist<-list.dirs(file.path(fdir,"DF_Surfaces"),recursive=F)
  dirlist<-dirlist[which(substring(basename(dirlist),1,6)==yearmon)]
  rlist<-list.files(dirlist,full.names=T,include.dirs=T,pattern="\\.tif$")
  rlist<-rlist[basename(rlist)==paste0(params,".tif")]
  
  if(length(rlist)==0){
    stop("No surface for this yearmon/param combination")
  }
  
  r<-raster::raster(rlist)
  
  stats<-ipdw::errorGen(r,dbdt,dbdt@data$chl.a.ug.l,plot=TRUE)[[1]]
  stats$days<-as.numeric(difftime(mean(streamdates),mean(as.POSIXct(dbdt@data$collection.date))))
    
  return(stats)
}
