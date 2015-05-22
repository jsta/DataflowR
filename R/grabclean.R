#'@name grabclean
#'@title cleaning grab data
#'@description if column names do not match up make sure there is only one "date" column
#'@param yearmon numeric survey date in yyyymm format
#'@param tofile logical save output to file
#'@export
#'@examples res<-grabclean(yearmon=201402,tofile=TRUE)
#'
grabclean<-function(yearmon,tofile=FALSE){
  
  fdir_fd<-file.path(fdir,"DF_FullDataSets")
  flist<-list.files(fdir_fd,include.dirs=T,full.names=T)
  dt<-read.csv(flist[substring(basename(flist),1,6)==yearmon])
  fdir_fd<-file.path(fdir,"DF_GrabSamples","Raw")
  flist<-list.files(fdir_fd,include.dirs=T,full.names=T,pattern=".csv")
  sumpath<-suppressWarnings(flist[which(as.numeric(substring(basename(flist),1,6))==yearmon)])
  
  #get col names####
  nms1<-read.csv(sumpath,sep=",",skip=2,header=F,stringsAsFactors=T,na.strings="",strip.white=T)[1,1:23]
  nms1<-apply(nms1,2,as.character)#LABID thru PP
  nms2<-read.csv(sumpath,sep=",",skip=3,header=F,stringsAsFactors=T,na.strings="",strip.white=T)
  #nms2<-nms2[1,25:28]
  nms2<-nms2[1,24:30]
  nms2<-apply(nms2,2,as.character)#Nitrogen
  nms3<-read.csv(sumpath,sep=",",skip=2,header=F,stringsAsFactors=T,na.strings="",strip.white=T)
  #nms3<-nms3[1,29:38]
  nms3<-nms3[1,32:38]
  nms3<-apply(nms3,2,as.character)#P and N
  nms4<-read.csv(sumpath,sep=",",skip=1,header=F,stringsAsFactors=T,na.strings="",strip.white=T)
  if(ncol(nms4)>38){#is there C6 data?
    nms4<-nms4[1,39:ncol(nms4)]
    nms4<-apply(nms4,2,as.character)
    nms.full<-c(nms1,NA,nms2,nms3,nms4)
  }else{
    nms.full<-c(nms1,NA,nms2,nms3)
    nms.full<-c(nms.full,rep(NA,(77-length(nms.full))))
  }
  
  nms.full<-tolower(nms.full)
  nms.full<-gsub(" ","",nms.full)
  nms.full<-gsub("\\(","",nms.full)
  nms.full<-gsub(")","",nms.full)
  
  #clean grabs####
  grabs<-read.csv(sumpath,sep=",",skip=5,header=F,stringsAsFactors=F,na.strings="",strip.white=T)
  grabs<-grabs[!is.na(grabs[,5]),]
  names(grabs)<-nms.full
  grabs<-grabs[,colSums(is.na(grabs))<nrow(grabs)]#remove columns of all NA
  grabs[,4]<-gsub("/","",grabs[,4])
  
  grabs<-grabs[,!is.na(names(grabs))]
  
  narows<-as.numeric(which(apply(grabs,1,function(x)sum(is.na(x)))>50))
  if(length(narows)>0){
    grabs<-grabs[-unique(narows),]#remove na rows
  }
  grabs<-grabs[!is.na(grabs[,8]),]  #eliminate EBs
  
  if(nchar(as.character(grabs[1,5]))<5){
    names(grabs)[4:5]<-c("date","time")#make sure to match parameter names for stations and dt
    grabs[,4:5]<-apply(grabs[,4:5],2,as.numeric)
    stations<-grabs[,4:5]
  }else{
    names(grabs)[5:6]<-c("date","time")#make sure to match parameter names for stations and dt
    grabs[,5:6]<-apply(grabs[,5:6],2,as.numeric)
    stations<-grabs[,5:6]
  }
  stations<-data.frame(apply(stations,2,function(x) as.numeric(as.character(x))))
  
  #clear any streaming data already entered
  if(ncol(grabs)>24){
  grabs<-grabs[,1:25]
  }
  
  #get averages that correspond to grab samples####
  #stream<-merge(stations,dt,all.x=T)#cuts dt down to match "stations"
  names(dt)<-tolower(names(dt))
  if((which(names(dt)=="date")-1)>0){
  dt<-dt[,-(1:(which(names(dt)=="date")-1))]#remove beginning junk columns
  }
  names(dt)[names(dt)=="chla"]<-"chlaiv"
  
  stream<-merge(stations,dt)#cuts dt down to match "stations"
  stream$date<-as.numeric(stream$date)
  
  #check to make sure streaming data exists for each grab
  nostream<-data.frame(matrix(NA,nrow=1,ncol=ncol(grabs)))
  names(nostream)<-names(grabs)
  cnt<-0
    
  stationsave<-stations
  #stations<-stationsave
  
  #change to create a list of lines to remove rather than updating within loop
  rmlist<-list()
  for(j in 1:nrow(stations)){
    #print(j)
    #j=14
    #print(cbind(stations[j,1],stations[j,2]))
    if(all(is.na(match(paste(stream$date,stream$time),paste(stations[j,1],stations[j,2]))))){
      if(cnt==0){
        nostream[1,]<-grabs[j,]
        cnt=1
      }else{
        nostream<-rbind(nostream,grabs[j,])
        #print(j)
      }
      rmlist[[j]]<-j
      }
  }
  if(length(unlist(rmlist))>0){
  stations<-stations[-unlist(rmlist),]
  grabs<-grabs[-unlist(rmlist),]
  }
  
  
  #identical(stationsave,stations)#should be false
  #nrow(grabs)==13#should be True
  #nrow(stations)==13#should be True
  
  #make sure nrow stations and nrow stream2 match
  stream$date<-as.numeric(stream$date)
  stream<-merge(stations,dt)#cuts dt down to match "stations"
  
  
  stream2<-data.frame(matrix(NA,nrow=nrow(stations),ncol=ncol(stream)))
  for(m in 1:ncol(stream)){
    #m=1
    #length(unique(paste(stream$date,stream$time)))
    #length(round(aggregate(stream[,m],by=list(stream$date,stream$time),mean)[,3],2))
    if(class(stream[,m])=="numeric"){
      stream2[,m]<-round(aggregate(stream[,m],by=list(stream$date,stream$time),mean)[,3],2)
    }else{
      stream2[,m]<-aggregate(stream[,m],by=list(stream$date,stream$time),Mode)[,3]
    }
  }
  names(stream2)<-names(stream)
  
    #stream2<-data.frame(aggregate(stream,by=list(stream$time),mean))#minute averages
  stream3<-stream2[order(stream2$date,stream2$time),]
  sname<-which(names(stream3)=="chlaiv")
  ename<-which(names(stream3)=="lat_dd")
  stream4<-cbind(stations,stream3[,c(sname:ename)])
  grabsfull<-cbind(grabs,stream4)
  
  #add back in grabs with missing streaming data
  if(any(!is.na(nostream[1,]))){
  nostream<-nostream[!is.na(nostream[,4]),]
  if((ncol(grabsfull)-ncol(nostream))!=0){
  padna<-data.frame(matrix(NA,nrow=nrow(nostream),ncol=(ncol(grabsfull)-ncol(nostream))))
  #if(any(match(names(nostream),names(grabsfull))[1:ncol(nostream)]!=1:ncol(nostream))){
  #  stop("problem with column names")
  #}
  nostream<-cbind(nostream,padna)
  }
  names(nostream)<-names(grabsfull)
  test<-rbind(grabsfull[1,],nostream)
  grabsfull<-rbind(grabsfull,nostream)
  }
  
  namestemp<-c("date","time","location","salt","chla","chla.1","tss","tss.1","pp","pp.1","tp","tdp","po4","toc","doc","tkn","tdkn","chlaiv","temp","cond","sal","trans","cdom","brighteners","phycoe","phycoc","c6chl","c6cdom","c6turbidity","c6temp","lon_dd","lat_dd")
  nseq<-seq(1,length(namestemp),1)
  
  namesalias<-read.table(text="chlorophyll.a,c6chl
c6chla,c6chl
spcondms,spcond
turbidity,c6turbidity",sep=",")
  namesalias<-apply(namesalias,2,function(x) as.character(x))
  
  #match dt names to a template that includes all possible columns####
  for(n in 1:ncol(grabsfull)){
    if(any(names(grabsfull)[n]==namesalias[,1])){
      names(grabsfull)[n]<-namesalias[which(names(grabsfull)[n]==namesalias[,1]),2]
    }
  }
  
  #trim extra columns and match order to template
  nmiss<-nseq[!(nseq %in% match(names(grabsfull),namestemp))]
  if(length(nmiss)>0){
    for(j in 1:length(nmiss)){
      grabsfull[,ncol(grabsfull)+1]<-NA
      names(grabsfull)[ncol(grabsfull)]<-namestemp[nmiss[j]]
    }
  }
  grabsfull<-grabsfull[,match(namestemp,names(grabsfull))]#sort to match order of namestemp
  
  grabsfull[,5:ncol(grabsfull)]<-suppressWarnings(apply(grabsfull[,5:ncol(grabsfull)],2,function(x) as.numeric(x)))
  grabsfull$flags<-NA
  
  if(tofile==TRUE){
    write.csv(grabsfull,file.path(fdir,"DF_GrabSamples",paste(yearmon,"j.csv",sep="")))
  }
    
  #if(ncol(read.csv(sumpath,sep=",",skip=4,header=F,stringsAsFactors=T,na.strings="",strip.white=T))>38){
  # final<-read.csv(sumpath,sep=",",skip=4,header=F,stringsAsFactors=T,na.strings="",strip.white=T)
  #names(final)<-nms.full
  #  return(final)
  #}else{
    return(grabsfull)
  #}
}

#go to ExtractedChl.R?

#combine all grabsamples####unfinished because of inconsistent column naming
# grablist<-list.files(path = file.path(paste(getwd(),"/DF_GrabSamples",sep="")), pattern = "*.csv", full.names = T)
# #grablist<-grablist[c(9,27,18,6)]#subset
# nms<-list()
# grabs<-cleangrab(grablist[[1]])
# nms[[1]]<-names(grabs)
# for(i in 2:length(grablist)){
# grabs<-rbind(grabs,cleangrab(grablist[[i]]))
# nms[[i]]<-names(cleangrab(grablist[[i]]))
# }
