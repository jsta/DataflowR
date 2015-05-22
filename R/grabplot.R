#'@name grabplot
#'@title plotting grab data
#'@export
#'@examples grabplot(rnge=201410,params=c("sal","chla"),plottype="permitbarplot")

#grabplot.R - Plotting Grab Data

grabplot<-function(rnge=c(201407,201410),params=c("sal","chla"),plottype="permitbarplot"){#a file with all msrmnts does not yet exist
  
  dt<-grabget(rnge)
  
  dt$location<-factor(dt$location,levels=c("Manatee","Barnes","Long Sound","Blackwater","Lake Surprise","Joe Bay","Pond 5","Little Madeira","Madeira","Seven Palm","Monroe","Terrapin","Whipray","Rankin"))
    
 if(plottype=="permitbarplot"){
   permitbarplot(dt,params)
 } 
}

permitbarplot<-function(dt,params){
  #params <- c("salt","chl.a","tss","pp","tp","tdp","po4","toc","doc","tkn","tdkn","chla","temp","cond","sal","trans","cdom","brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","c6temp")
  #ylab <- "Chl"
  dt$datecode<-cutree(hclust(dist(dt$date)),h=1000)#may need to adjust h
  
  for(param in params){
  ylab<-param
  if(!(param %in% colnames(dt))){
    stop("parameters not found in data frame, check spelling?")
  }
  means<-tapply(dt[,param],list(dt$datecode,dt$location),function(x) mean(x,na.rm=T))
  if(any(is.na(colSums(means)))){
  means<-means[,-which(is.na(colSums(means)))]#remove columns of all NA
  }
  if(paste(param,".1",sep="") %in% colnames(dt)){
  stdev<-tapply(dt[,paste(param,".1",sep="")],list(dt$datecode,dt$location),function (x) mean(x,na.rm=T))
  standardErrors<-(stdev/sqrt(2))
  standardErrors<-standardErrors[,apply(standardErrors,2,function(x) !any(is.na(x)))]#temporary until I fix "$location" errors
  plotTop <- max(means+standardErrors*2,na.rm=TRUE)
  nosd<-FALSE
  }else{
    plotTop<-max(means)
    nosd<-TRUE
  }
  
  barCenters <- barplot(means, beside=T, las=2, ylim=c(0,plotTop),cex.names=0.8,col=terrain.colors(length(unique(dt$datecode))),ylab=ylab)
  if(nosd==FALSE){
  segments(barCenters, means, barCenters, means+standardErrors*2, lwd=2) 
  suppressWarnings(arrows(barCenters, means, barCenters, means+standardErrors*2, lwd=1, angle=90, code=2,length=0.02))
  }
  
  labs<-list()
  for(i in unique(dt$datecode)){
    labs[[i]]<-dt$date[dt$datecode==unique(dt$datecode)[i]][1]
  }
  datelabs<-as.character(unlist(labs))
  #legend(locator(1),datelabs,fill=terrain.colors(5),cex=0.7)
  legend("topleft",datelabs,fill=terrain.colors(length(unique(dt$datecode))),cex=0.7)
  }
  
}
  