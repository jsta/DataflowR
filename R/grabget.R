#'@name grabget
#'@title aggregate grab data
#'@export
#'@param rnge list of length two specifying date range in yyyymm format
#'@param remove.flags logical trim dataset based on flags?
#'@param fdir character file path to local data directory
#'@examples \dontrun{ 
#'grabs<-grabget(rnge=c(201402,201410))
#'}
#'rnge<-201308

grabget<-function(rnge,remove.flags=FALSE,fdir=getOption("fdir")){
  #rnge<-200910
  ##i=2
  ##rnge=c(201402,201410)
  if(length(rnge)==1){
    rnge<-c(rnge,rnge)
  }
    
  agglist<-list.files(file.path(fdir,"DF_GrabSamples"),"*.csv",include.dirs=T,full.names=T)
  agglist<-suppressWarnings(agglist[which(as.numeric(substring(basename(agglist),1,6)) >= rnge[1])])
  agglist<-suppressWarnings(agglist[which(as.numeric(substring(basename(agglist),1,6)) <= rnge[2])])
  
  dtlist<-list()
  for(i in 1:length(agglist)){
    #i<-1
    dt <- read.csv(agglist[i],header=T, stringsAsFactors = FALSE)[,-1]
    names(dt)<-tolower(names(dt))
    
    if(any(duplicated(names(dt)))){
      dupname<-names(dt)[duplicated(names(dt))]
      print(dupname)
      print("check duplicated names")
      #for(p in 1:length(dupname)){#remove duplicated columns
      #  maxname<-which.max(lapply(which(names(dt)==dupname[p]),function(x) sum(!is.na(dt[,x]))))
      #  dt<-dt[,-which(names(dt)==dupname[p])[maxname]]
      #}
    }
    if(remove.flags==TRUE){
      if(any(names(dt)=="flags")){
        if(length(which(dt$flags=="s"))>0){
      dt<-dt[-which(dt$flags=="s"),]
        }
      }
    }
    
    dtlist[[i]]<-dt
    #print(paste0(agglist[i],ncol(dt)))#
}

dtall<-do.call(rbind,dtlist)

dtall
}
