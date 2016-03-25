#'@name streamclean
#'@title Cleaning raw streaming Dataflow output
#'@description Cleaning raw streaming Dataflow, C6, Eureka-Manta, and YSI-Exo output
#'@param yearmon numeric designation of survey date formatted as yyyymm
#'@param mmin minimum measurement frequency (# measurements/min)
#'@param c6mmin integer optional minimum c6 measurement frequency if different than mmin
#'@param c6pres logical was C6 data collected?
#'@param eummin integer optional minimum eureka measurement frequency if different than mmin
#'@param eupres logical was eureka data collected?
#'@param exommin integer optional minimum exo measurement frequency if different than mmin
#'@param exopres logical was exo data collected?
#'@param tofile logical save cleaned output to DF_FullDataSets?
#'@param sep character optional predesignation of item seperation character in raw data files
#'@param fdir character file path to local data directory
#'@export
#'@importFrom rgdal readOGR
#'@importFrom zoo zoo na.approx
#'@importFrom sp coordinates CRS spTransform
#'@details Dataflow cleaning drops all minutes that have less measurements than "mmin". C6 data is interpolated to match Dataflow.  Automatically compares salinity against conducitivty/temperature recalculated salinity and replaces if slope of fit is not close to 1. Bad DO columns must sometimes be removed manually. TODO - Add check the make sure that the year of the data (not just the filename) matches the year of yearmon
#'@examples \dontrun{
#'dt <- streamclean(yearmon = 201505, mmin = 7, c6mmin = 10, tofile = FALSE, c6pres = TRUE)
#'dt <- streamclean(yearmon = 201513, mmin = 7, c6pres = TRUE, c6mmin = 12,
#' tofile = FALSE, exommin = 60, exopres = TRUE, eupres = TRUE, eummin = 12)
#' 
#' dt <- streamclean(yearmon = 201601, mmin = 12, c6pres = FALSE, eupres = TRUE, eummin = 12)
#' 
#'}

streamclean <- function(yearmon, mmin, c6mmin = NA, c6pres = TRUE, eummin = NA, eupres = FALSE, exommin = NA, exopres = FALSE, tofile = FALSE, sep = ",", fdir = getOption("fdir")){
  
  options(warn = -1)  
  fdir_fd <- file.path(fdir, "DF_FullDataSets", "Raw", "InstrumentOutput")
  flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
  flist <- flist[substring(basename(flist), 1, 6) == yearmon]
  
  dflist <- list.files(flist, pattern = c("*.txt"), include.dirs = T, full.names = T)
  if(length(dflist) == 0){
    dflist <- list.files(flist, pattern = c("*.TXT"), include.dirs = T, full.names = T)
  }
  if(length(dflist) == 0){
    dflist <- list.files(flist, pattern = c("*DF.csv"), include.dirs = T, full.names = T)
  }
  
  c6list <- list.files(flist, pattern = c("*C6.csv"), include.dirs = T, full.names = T)
  if(length(c6list) == 0){
    c6list <- list.files(flist, pattern = c("*C6.CSV"), include.dirs = T, full.names = T)
  }
  
  eulist <- list.files(flist, pattern = c("*eu.csv"), include.dirs = T, full.names = T)
  if(length(eulist) == 0 ){
    eulist <- list.files(flist, pattern = c("*eu.CSV"), include.dirs = T, full.names = T)
  }
  
  exolist <- list.files(flist, pattern = c("*exo.csv"), include.dirs = T, full.names = T)
  if(length(exolist) == 0 ){
    exolist <- list.files(flist, pattern = c("*exo.CSV"), include.dirs = T, full.names = T)
  }
  
  if(length(c6list) != length(dflist)){
    warning("Differing numbers of Dataflow and C6 input files")
  }
  
  set_empty_to_NULL <- function(x) if(length(x) == 0){NULL}else{x}
  dflist  <- set_empty_to_NULL(dflist)
  eulist  <- set_empty_to_NULL(eulist)
  exolist <- set_empty_to_NULL(exolist)
  c6list  <- set_empty_to_NULL(c6list)  
  
  survey_days <- unique(substring(basename(c(dflist, eulist, exolist, c6list)[sapply(c(dflist, eulist, exolist, c6list), function(x) length(x))]), 1, 8))
  
  reslist <- list()
  for(i in 1:length(survey_days)){
    sep <- ","
    dt <- read.csv(dflist[i], skip = 0, header = F, sep = sep)#start with comma sep
        
    if(suppressWarnings(nchar(gsub("\t","",dt[1,]))<nchar(as.character(dt[1,])))){#switch to tab sep
      sep<-"\t"
      dt<-read.csv(dflist[i],skip=0,header=F,sep=sep,stringsAsFactors=FALSE)
    }
    
    #detect beginning of measurements
    fskip=1
    while(!(!(class(dt[,1])!="integer") | !(class(dt[,1])!="numeric"))){
      dt<-read.csv(dflist[i],skip=fskip,header=F,sep=sep,stringsAsFactors=FALSE)
      #print(dt[1,1])
      if(!any(!is.na(dt[,1])) | mean(nchar(as.character(dt[,1])))<1 | sum(is.na(dt[,1]))>(nrow(dt)/2) | sum(nchar(gsub("_","",as.character(dt[,1])))-nchar(as.character(dt[,1])))!=0){
        dt<-dt[,-1]
      }
      fskip<-fskip+1
      if(fskip>20){
        stop(paste("Cannot find beginning of measurements!",dflist[i]))
      }
    }
    
    sep<-","
#==================================================================#
    
    #remove bad columns
    if(class(dt[,3]) == "integer"){
      dt<-dt[,-3]
      print("removing existing seconds column")
    }#remove existing sec column
    
    #remove bad columns of all 0 or NA
    if(ncol(dt) > 14){
      dtno.na <- dt[complete.cases(dt[,1:12]),]
      dt <- dt[,apply(dtno.na, 2, function(x) abs(sum(as.numeric(x), na.rm = T)) != 0)]
    }
    
    #temp should never be less than 10, these are likely 'bad' DO columns?
    if(mean(as.numeric(dt[,4]),na.rm=T)<10 & mean(as.numeric(dt[,5]),na.rm=T)<10){
      dt<-dt[,-4:-5]
    }
    
    dt <- dt[,apply(dt, 2, function(x) abs(sum(as.numeric(x), na.rm = T)) > 22)]#take out all 0 (22 or 38 is an arbitrary "tolerance" value)
    ones<-apply(dt,2,function(x) sd(as.numeric(x)[as.numeric(x)!=0 & !is.na(as.numeric(x))]))!=0
    ones[is.na(ones)]<-TRUE
    ones[1:2]<-TRUE
    dt<-dt[,ones]
    dt<-dt[,apply(dt,2,function(x) mean(nchar(x),na.rm=T))>=3.0] #(3 is an arbitrary "tolerance" value)
    dt<-dt[,apply(dt[,3:ncol(dt)],2,function(x) length(unique(x))!=3)]
    names(dt)<-c("date","time","chla","temp","cond","sal","trans","cdom","lat_dd","lon_dd")
    
    #convert factors to numeric
    dt<-data.frame(as.matrix(dt))
    factorToNumeric<-function(f) as.numeric(levels(f))[f]
    #check to make sure that there are any factor class columns
    if(any(sapply(dt,class)=="factor")){
      dt<-data.frame(sapply(dt,factorToNumeric))  
    }
    
    #fix lon lat formatting
    if(mean(nchar(as.character(round(dt[,"lat_dd"]))))!=2){
      lat<-dt[,"lat_dd"]
      latdeg<-as.numeric(substr(lat,0,2))
      latmin<-as.numeric(substr(lat,3,8))
      dt[,"lat_dd"]<-latdeg+latmin/60
      lon<-dt[,"lon_dd"]
      londeg<-as.numeric(substr(lon,0,2))
      lonmin<-as.numeric(substr(lon,3,8))
      dt[,"lon_dd"]<-(londeg+lonmin/60)*-1
    }
    
    dt$time<-as.numeric(dt$time)
    dt$date<-as.numeric(dt$date)
    
    #remove rows of all NA values
    dt<-dt[as.numeric(rowSums(is.na(dt)))<ncol(dt)-1,]
    dt<-dt[as.numeric(rowSums(is.na(dt[,c("lat_dd","lon_dd")])))<2,]
    
    #remove unrealistic coordinates
    dt <- dt[abs(dt$lat_dd) > 24.5 & abs(dt$lat_dd) < 25.5, ]
    dt <- dt[abs(dt$lon_dd) > 80.1 & abs(dt$lon_dd) < 82, ]
    
    #check for incomplete minutes
    datelist<-unique(dt$date)
    reslist2<-list()
    for(j in 1:length(datelist)){
      #j<-1
      curdat<-dt[dt$date==datelist[j],]
      gdata<-data.frame(table(curdat$time))
      fdata<-as.numeric(as.character(gdata[gdata$Freq<mmin,1]))#too few measurements
      odata<-as.numeric(as.character(gdata[gdata$Freq>mmin,1]))#too many measurements
      if(length(odata)>0){
        for(k in 1:length(odata)){
          #k<-1
          leng<-nrow(curdat[curdat$time==odata[k],])
          remo<-sample(as.numeric(row.names(curdat[curdat$time==odata[k],])),leng-mmin)
          curdat<-curdat[-match(remo,as.numeric(row.names(curdat))),]
        }}
      curdat<-curdat[!curdat$time %in% fdata,]
      curdat<-curdat[!is.na(curdat$time),]
      curdat$sec<-rep(round(seq(from=0,to=60-60/mmin,by=60/mmin),3),times=nrow(data.frame(table(curdat$time))))
      
      reslist2[[j]]<-curdat
    }
    dt<-do.call("rbind",reslist2)
    
    #detect when dt is measured at fractional seconds
    if(!identical(round(dt$sec),dt$sec)){
      dt$sec<-round(dt$sec)
    }
    
    #create POSIXct datetime column
    yr<-substring(dt$date,nchar(dt$date)-1,nchar(dt$date))
    day<-substring(dt$date,nchar(dt$date)-3,nchar(dt$date)-2)
    mon<-substring(dt$date,1,nchar(dt$date)-4)
    hr<-substring(dt$time,1,nchar(dt$time)-2)
    min<-substring(dt$time,nchar(dt$time)-1,nchar(dt$time))
    
    if(mean(nchar(mon))==1){mon<-paste("0",mon,sep="")}
    dt$datetime<-paste(yr,"-",mon,"-",day,"-",hr,"-",min,"-",dt$sec,sep="")
    rm(min)
    dt$datetime<-as.POSIXct(strptime(dt$datetime,format="%y-%m-%d-%H-%M-%S"))
        
    #clean data frame
    #trim beginning and end based on when data is all zeros
    trimdt<-function(dt){
      j=1  
      for(i in 1:nrow(dt)){
        if(dt[i,1:9][order(dt[i,1:9])][2]>0){
          break
        }
        j=i+1
        }
      k<-nrow(dt)
      for(i in nrow(dt):1){
        if(!is.na(min(dt[i,1:9]))>0){
          break
        }
        k=i-1
        }
      dt[j:k,]
    }
    dt<-trimdt(dt)
    
    #check for correct cond to salinity calculations
    #browser()
    corsal<-DataflowR::cond2sal(dt$cond*1000,dt$temp)
    if((lm(corsal~dt$sal)$coefficients[2]-1)>0.02){
      dt$sal<-corsal
    }
  
    print(paste(basename(dflist[i]),"processed",sep=" "))
    
    #clean and merge C6 here===================================================#
    if(c6pres==TRUE){
            c6dfmatch<-which(substring(basename(c6list),1,8)==substring(basename(dflist[i]),1,8))
            #browser()
            if(length(c6dfmatch) != 0){
              c6<-read.csv(c6list[c6dfmatch],skip=12,header=F)[,1:9]
              names(c6)<-c("datec","brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","depth","c6temp")
              if(!any(!is.na(c6[,"c6temp"]))){
                c6<-c6[,-9]
                names(c6)[8]<-"c6temp"
              }else{
                c6<-c6[,-8]
              }
                                        
              #check for mismatched c6df measurement frequency
              dtfreq<-max(dt$sec[2]-dt$sec[1],dt$sec[3]-dt$sec[2])
              
              #check for missing seconds information
              if(all(is.na(as.POSIXct(strptime(c6$datec,"%m/%d/%y %H:%M:%S"))))){
                c6sec <- unlist(lapply(rle(sapply(c6$datec, function(x) strftime(strptime(x, format = "%m/%d/%Y %H:%M"), format = "%M")))$lengths, function(x) seq(0, 60 - (60/x), length.out = x)))
                c6$datec <- as.POSIXct(strptime(paste0(c6$datec, ":", c6sec), "%m/%d/%Y %H:%M:%S"))
              }else{
              c6$datec<-as.POSIXct(strptime(c6$datec,"%m/%d/%y %H:%M:%S"))
              }
              
              if(!any(!is.na(c6$datec))){
                c6<-read.csv(c6list[c6dfmatch],skip=12,header=F)[,1:9]
                names(c6)<-c("datec","brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","depth","c6temp")
                if(!any(!is.na(c6[,"c6temp"]))){
                c6<-c6[,-9]
                names(c6)[8]<-"c6temp"
                }else{
                  c6<-c6[,-8]
                }
                
                if(nchar(strsplit(as.character(c6[,"datec"]),"/")[[1]][1])==1){
                  if(rle(as.character(c6[,"datec"]))$length[1]!=c6mmin){#account for less than full minute to start
                    c6<-c6[c6$datec!=rle(as.character(c6[,"datec"]))$values[1],]
                  }
                  if(!exists("c6mmin")){
                    c6mmin<-mmin
                  }
                  
                padm.addsec<-function(x,c6mmin){
                    #x<-c6[,"datec"]
                    #pad month
                    x<-as.character(x)
                    x<-paste("0",substring(x,0,1),substring(x,2,nchar(x)),sep="")
                    
                    #add sec
                    sseq<-seq(0,60-(60/c6mmin),60/c6mmin)
                    sseq<-sapply(sseq,function(x) ifelse(nchar(x)==1,paste("0",x,sep=""),x))
                    paste(x,":",sseq,sep="")
                  }
                  c6[,"datec"]<-padm.addsec(c6[,"datec"],c6mmin=c6mmin)
                }
                c6$datec<-as.POSIXct(strptime(c6$datec,"%m/%d/%Y %H:%M:%S"))
                
              }
              c6$sec<-as.numeric(format(c6$datec,'%S'))
              c6freq<-c6$sec[2]-c6$sec[1]
              c6$datec<-as.POSIXct(c6$datec)
              c6zoo<-zoo::zoo(c6,c6$datec)
              
              #if true generate second-wise c6 zoo object
              if(dtfreq!=c6freq){
                seqexpand <- data.frame(seq(as.POSIXct(unique(as.Date(c6$datec))), as.POSIXct(unique(as.Date(c6$datec)) + 1), 1))
                names(seqexpand)<-c("datec")
                seqexpand <- zoo::zoo(seqexpand, seqexpand$datec)
                seqexpand<-merge(seqexpand,c6zoo)
                seqexpand<-seqexpand[min(which(!is.na(seqexpand[,"brighteners"]))):max(which(!is.na(seqexpand[,"brighteners"]))),3:ncol(seqexpand)]
                idx<-colSums(!is.na(seqexpand))>1
                alignset<-zoo::na.approx(seqexpand[,idx])
              
              
              #align based on temp(not implemented yet)
              alignset<-alignset[!is.na(alignset$c6cdom),]
              
              #dtsave<-dt
              #c6save<-c6
              
              #merge c6 and df based on time stamp
              c6<-data.frame(alignset,zoo::index(alignset),row.names=NULL)
              names(c6)[ncol(c6)]<-"datetime"
              #browser()
              dt<-merge(dt,c6,by="datetime",all.x=T)
              datetime<-dt[,"datetime"]
              dt<-data.frame(apply(dt[,2:ncol(dt)],2,function(x) as.numeric(as.character(x))))
              dt$datetime<-datetime                   
              }else{
                
                #remove c6 rows with duplicate datetime stamps
                #test<-c6[!duplicated(c6$datetime),]
                                
                datetime<-dt$datetime
                names(c6)[names(c6)=="datec"]<-"datetime"
                #browser()
                dt<-merge(dt,c6,by="datetime",all.x=T)
                
                dt<-data.frame(apply(dt[,2:ncol(dt)],2,function(x) as.numeric(as.character(x))))
                dt$datetime<-datetime  
              }
              
              print(paste(basename(c6list[c6dfmatch]),"processed",sep=" "))
              }else{
                print("No C6 file detected, constructing empty data matrix")
                emptyc6<-data.frame(matrix(NA,nrow=nrow(dt),ncol=20-ncol(dt)))
                names(emptyc6)<-c("brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","c6temp","sec.y")
                dt <- cbind(dt,emptyc6)
                names(dt)[names(dt) == "sec"] <- "sec.x"
                dt<-dt[,match(names(streamget(201505)),names(dt))[!is.na(match(names(streamget(201505)),names(dt)))]]
              }
    }
    
    #clean and merge eu here===============================================#
    if(eupres == TRUE){
      eudfmatch <- which(substring(basename(eulist), 1, 8) == substring(basename(dflist[i]), 1, 8))
      if(length(eudfmatch) != 0){
        eu <- read.csv(eulist[eudfmatch], header = TRUE, stringsAsFactors = FALSE)#[,c(1:13, 17:18)]
        names(eu) <- tolower(make.names(names(eu)))
        eu$datetime <- paste(sapply(eu$date, function(x) mdy2mmyyyy(x)), eu$time)
        
        #check for mismatched eu-df measurement frequency
        dtfreq <- max(dt$sec.x[2] - dt$sec.x[1], dt$sec.x[3] - dt$sec.x[2])
        
        #check for missing seconds information
        if(all(is.na(as.POSIXct(strptime(eu$datetime, "%m/%d/%Y %H:%M:%S"))))){
          eusec <- unlist(lapply(rle(sapply(eu$date, function(x) strftime(strptime(x, format = "%m/%d/%Y %H:%M"), format = "%M")))$lengths, function(x) seq(0, 60 - (60 / x), length.out = x)))
          eu$date <- as.POSIXct(strptime(paste0(eu$date, ":", eusec), "%m/%d/%Y %H:%M:%S"))
        }else{
          eu$datetime <- as.POSIXct(strptime(eu$datetime,"%m/%d/%Y %H:%M:%S"))
          eu$datetime <- as.POSIXct(strptime(eu$datetime, "%Y-%m-%d %H:%M:%S"))
          eu <- eu[!is.na(eu$datetime) & nchar(eu$datetime) == 10,]
        }
        
        eufreq <- Mode(diff(as.numeric(format(eu$datetime, '%S'))))
        eu$datetime <- as.POSIXct(eu$datetime)
        #browser()
        euzoo <- zoo::zoo(eu, eu$datetime)
        
        #if true generate second-wise c6 zoo object
        if(dtfreq != eufreq){
          #eu <- eu[strftime(eu$datetime, format = "%Y-%m-%d") == unique(strftime(eu$datetime, format = "%Y-%m-%d"))[1],]
          
          seqexpand <- data.frame(
            seq(
            as.POSIXct(as.Date(unique(strftime(eu$datetime, format = "%Y-%m-%d"))[1])),
            as.POSIXct(as.Date(unique(strftime(eu$datetime, format = "%Y-%m-%d"))[1]) + 1), 1))
          
          names(seqexpand) <- c("datetime")
          seqexpand <- zoo::zoo(seqexpand, seqexpand$datetime)
          seqexpand <- merge(seqexpand, euzoo)
          seqexpand <- seqexpand[min(which(!is.na(seqexpand[,2]))) : max(which(!is.na(seqexpand[,2]))),]
          idx <- colSums(!is.na(seqexpand)) < nrow(seqexpand)
          idx[c(1:3, ncol(seqexpand))] <- FALSE
          if(any(colSums(!is.na(seqexpand)) == 0)){
            idx[which(colSums(!is.na(seqexpand)) == 0)] <- FALSE
          }
          
          alignset <- zoo::na.approx(seqexpand[,idx])
          
          #merge eu and df based on time stamp
          eu <- data.frame(alignset, zoo::index(alignset), row.names = NULL)
          names(eu)[ncol(eu)] <- "datetime"
          dt <- merge(dt, eu, by = "datetime", all.x = T)
          datetime <- dt[,"datetime"]
          dt <- data.frame(apply(dt[,2:ncol(dt)], 2, function(x) as.numeric(as.character(x))))
          dt$datetime <- datetime                   
        }else{
          
#           datetime <- dt$datetime
#           names(c6)[names(c6)=="datec"]<-"datetime"
#           #browser()
#           dt<-merge(dt,c6,by="datetime",all.x=T)
#           
#           dt<-data.frame(apply(dt[,2:ncol(dt)],2,function(x) as.numeric(as.character(x))))
#           dt$datetime<-datetime  
        }
#         
         # print(paste(basename(c6list[c6dfmatch]),"processed",sep=" "))
      }else{
#         print("No C6 file detected, constructing empty data matrix")
#         emptyc6<-data.frame(matrix(NA,nrow=nrow(dt),ncol=20-ncol(dt)))
#         names(emptyc6)<-c("brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","c6temp","sec.y")
#         dt<-cbind(dt,emptyc6)
#         names(dt)[names(dt)=="sec"]<-"sec.x"
#         dt<-dt[,match(names(streamget(201505)),names(dt))[!is.na(match(names(streamget(201505)),names(dt)))]]
      }
      print(paste(basename(eulist[eudfmatch]),"processed",sep=" "))
    }
    
    #clean and merge exo here==============================================#
    if(exopres == TRUE){
      exodfmatch <- which(substring(basename(exolist), 1, 8) == substring(basename(dflist[i]), 1, 8))
      if(length(exodfmatch) != 0){
        exo <- read.csv(exolist[exodfmatch], header = T, skip = 12, stringsAsFactors = FALSE)
      
      clean_exo <- function(exo){
        names(exo) <- tolower(gsub("\\.", "", names(exo)))
        exo <- exo[,!(names(exo) %in% c("timefractsec", "sitename","x"))]
        exo$datetime <- as.POSIXct(strptime(paste(
          exo[,grep("date", names(exo))],
          exo[,grep("time", names(exo))]), format = "%m/%d/%Y %H:%M:%S")
        )
        exo$sec <- strftime(exo$datetime, "%S")
        exo <- exo[which(!duplicated(exo$datetime)),]
        exo[,"longitudedegrees"] <- exo[,"longitudedegrees"] * -1
        exo <- exo[exo$latitudedegrees > 24,]
        exo
      }
      
      exo <- clean_exo(exo)
      
      #check for mismatched exodf measurement frequency
      dtfreq <- max(dt$sec.x[2] - dt$sec.x[1], dt$sec.x[3] - dt$sec.x[2])
      
      exofreq <- Mode(diff(as.numeric(format(exo$datetime, '%S'))))
      exozoo <- zoo::zoo(exo, exo$datetime)
      
      #if true generate second-wise exo zoo object
      if(dtfreq != exofreq){
        
        seqexpand <- data.frame(
          seq(
            as.POSIXct(as.Date(unique(strftime(exo$datetime, format = "%Y-%m-%d"))[1])),
            as.POSIXct(as.Date(unique(strftime(exo$datetime, format = "%Y-%m-%d"))[1]) + 1), 1)
          )
        
        names(seqexpand) <- c("datetime")
        seqexpand <- zoo::zoo(seqexpand, seqexpand$datetime)
        seqexpand <- merge(seqexpand, exozoo)
        seqexpand <- seqexpand[min(which(!is.na(seqexpand[,2]))):max(which(!is.na(seqexpand[,2]))),]
        idx <- colSums(!is.na(seqexpand)) < nrow(seqexpand)
        idx[c(1:3, (ncol(seqexpand) - 1), ncol(seqexpand))] <- FALSE
        
        alignset <- zoo::na.approx(seqexpand[,idx])
        exo <- data.frame(alignset, zoo::index(alignset),
                            row.names = NULL)
        names(exo)[ncol(exo)]<-"datetime"
        
        dt <- merge(dt, exo, by = "datetime", all.x = T)
        datetime <- dt[,"datetime"]
        dt <- data.frame(apply(dt[,2:ncol(dt)],2,function(x) as.numeric(as.character(x))))
        dt$datetime <- datetime
      }else{
        stop("not implemented yet")
      }
      
      print(paste(basename(exolist[exodfmatch]),"processed",sep=" "))
      }else{
        stop("No exo file detected, set exopres = FALSE?")
      }
    }
    
              
    #create basin designations here
        
    #define projections
    projstr <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    latlonproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    fathombasins <- rgdal::readOGR(file.path(fdir,"DF_Basefile/fathom_basins_proj.shp"),layer="fathom_basins_proj",verbose=FALSE)
    cerpbasins<-rgdal::readOGR(file.path(fdir,"DF_Basefile/fbfs_zones.shp"),layer="fbfs_zones",verbose=FALSE)
    selectiongrid<-rgdal::readOGR(file.path(fdir,"DF_Basefile/testgrid3.shp"),layer="testgrid3",verbose=FALSE)
    
    #spatial join
    dt<-dt[!is.na(dt$lat_dd)&!is.na(dt$lon_dd),]
    xy<-cbind(dt$lon_dd,dt$lat_dd)
    xy<-data.frame(xy)
    sp::coordinates(dt)<-~lon_dd+lat_dd#will throw error if latlon has NA
    sp::proj4string(dt)<-sp::CRS(latlonproj)
    dt<-sp::spTransform(dt,sp::CRS(projstr))
    fulldataset<-dt
    
    fulldataset.over<-sp::over(fulldataset,selectiongrid)
    fulldataset.over2<-sp::over(fulldataset,fathombasins[,1:2])
    fulldataset.over3<-sp::over(fulldataset,cerpbasins[,2])
    fulldataset.over<-cbind(data.frame(fulldataset),data.frame(fulldataset.over),data.frame(fulldataset.over2),data.frame(fulldataset.over3))
    fulldataset.over$lon_dd<-xy[,1]
    fulldataset.over$lat_dd<-xy[,2]
    
    fulldataset.over<-fulldataset.over[,names(fulldataset.over)!="NA."]
    reslist[[i]]<-fulldataset.over
  }
  
  dt<-do.call("rbind",reslist)
  names(dt)<-tolower(names(dt))

  if(tofile==TRUE){
  #add check to verify yearmon before overwriting
    dtname<-file.path(fdir,.Platform$file.sep,"DF_FullDataSets",.Platform$file.sep,substring(basename(dflist[1]),1,6),"j.csv",fsep="")
    if(!file.exists(dtname)){
      write.csv(dt,dtname,row.names = FALSE)
    }else{
      stop("overwrite file?")
    }
  }
options(warn=0)

dt
}

#'@name streamget
#'@title Retrieve previously cleaned full streaming datasets
#'@description Retrieve previously cleaned full streaming datasets
#'@param yearmon numeric date in yyyymm format
#'@param fdir character file path to local data directory
#'@param qa logical strip flagged data?
#'@export
#'@examples \dontrun{
#'yearmon<-201311
#'dt<-streamget(yearmon)
#'}

streamget<-function(yearmon,qa=TRUE,fdir=getOption("fdir")){
  fdir_fd <- file.path(fdir, "DF_FullDataSets")
  flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
  flist <- flist[substring(basename(flist),1,6) == yearmon]
  dt <- read.csv(flist, stringsAsFactors = FALSE)
  
  if(qa == TRUE && file.exists(file.path(fdir, "DF_FullDataSets", "QA", paste(yearmon, "qa.csv", sep = ""))) && identical(dim(dt), dim(read.csv(file.path(fdir, "DF_FullDataSets", "QA", paste(yearmon, "qa.csv", sep = "")))))){
    
    qafile<-read.csv(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep="")))
    
    if(!any(names(qafile)=="chlext")&any(names(dt)=="chlext")){
      qafile$chlext<-NA
    }
    
    if(!(identical(dim(qafile),dim(dt)))){
      warning("QA file dimensions do not match data dimensions")
    }
    dt[!is.na(qafile)]<-NA
  }
  dt
}

#'@name streamqa
#'@title Supervised quality control of streaming datasets
#'@description Supervised quality control of streaming datasets
#'@param yearmon numeric date in yyyymm format
#'@param setthresh logical set parameter thresholds
#'@param trimends logical look to trim ends of data stream? NOT IMPLEMENTED YET
#'@param paired logical examine relationships between paried parameters?
#'@param fdir file.path to data directory
#'@details loop through parameters giving the opportunity to trim measurement ends, set entire variables to NA, remove variables above/below a threshold
#'@return a matrix of the same size/shape of the fulldataset, with entries specifying where to set to NA, saved to DF_FullDataSets/Raw/IntrumentOutput
#'@export
#'@examples \dontrun{
#'dt<-streamqa(yearmon=201410)
#'}

streamqa<-function(yearmon,setthresh=TRUE,trimends=FALSE,paired=TRUE,fdir=getOption("fdir")){
  #yearmon=200904
  dt<-streamget(yearmon)
  
  dt<-dt[with(dt,order(date,time)),]
  
  if(setthresh==TRUE){
  if(file.exists(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep="")))){
    dtqa<-read.csv(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa.csv",sep=""))) 
  }else{  
    dtqa<-data.frame(matrix(NA,nrow=nrow(dt),ncol=ncol(dt)))
    names(dtqa)<-names(dt)
  }
  
  #explore and set parameter threshold limits
  par(mfrow=c(1,1))
  parset<-c("chla","temp","cond","sal","trans","cdom","brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","c6temp")
  
  parset<-parset[parset %in% names(dt)]
  
  
  for(i in parset){
    #i<-"sal"
  
      if(any(!is.na(dt[,i]))){
      
    #i<-"c6chla"
    plot(dt[,i],ylab=i)
    
    threshlog<-"c"
    thresh<-NA
    while(threshlog!="q"){
    threshlog<-readline(message("Set threshold? Enter an upper and lower range as c(lower,upper) or press 'q' to move to next QA step: ",appendLF=FALSE))
    if(!is.na(threshlog)&threshlog!="q"){
      thresh<-threshlog
    plot(dt[,i],ylab=i,ylim=eval(parse(text=threshlog)))
    }    
  }
  if(!is.na(thresh)){
    thresh<-gsub("c\\(","",thresh)
    thresh<-gsub(")","",thresh)
    thresh<-unlist(lapply(strsplit(thresh,","),as.numeric))
    dtqa[,i][dt[,i]<thresh[1]]<-"r"
    dt[,i][dt[,i]<thresh[1]]<-NA
    dtqa[,i][dt[,i]>thresh[2]]<-"r"
    dt[,i][dt[,i]>thresh[2]]<-NA
  }
  }
  }
  }
  
  #explore paired parameter relationships
  if(paired==TRUE){
  par(mfrow=c(3,1),mar=c(0,4,0,0))
  
    if(any(!is.na(dt[,"temp"]))&any(!is.na(dt[,"c6temp"]))){
  plot(dt[,"temp"],xaxt="n",xlab="")
  plot(dt[,"c6temp"],xaxt="n",xlab="")
  plot(dt[,"temp"],dt[,"c6temp"])
  abline(a=0,b=1,col="red")
  qalogical<-readline(message("Press 'Enter' to continue, '1' to set top panel to NA, '2' to set middle panel to NA: ",appendLF=FALSE))
  if(qalogical==1|qalogical=='1'){
    dt[,"temp"]<-NA
    dtqa[,"temp"]<-"r"
  }
  if(qalogical==2|qalogical=='2'){
    dt[,"c6temp"]<-NA
    dtqa[,"c6temp"]<-"r"
  }
    }
  
    if(any(!is.na(dt[,"chla"]))&any(!is.na(dt[,"c6chla"]))){
  plot(dt[,"chla"],xaxt="n",xlab="")
  plot(dt[,"c6chla"],xaxt="n",xlab="")
  plot(dt[,"chla"],dt[,"c6chla"])
  abline(lm(dt[,"c6chla"]~dt[,"chla"]),col="red")
  qalogical<-readline(message("Press 'Enter' to continue, '1' to set top panel to NA, '2' to set middle panel to NA: ",appendLF=FALSE))
  if(qalogical==1|qalogical=='1'){
    dt[,"chla"]<-NA
    dtqa[,"chla"]<-"r"
  }
  if(qalogical==2|qalogical=='2'){
    dt[,"c6chla"]<-NA
    dtqa[,"c6chla"]<-"r"
  }
    }
    
    if(any(!is.na(dt[,"cdom"]))&any(!is.na(dt[,"c6cdom"]))){
  plot(dt[,"cdom"],xaxt="n",xlab="")
  plot(dt[,"c6cdom"],xaxt="n",xlab="")
  plot(dt[,"cdom"],dt[,"c6cdom"])
  abline(lm(dt[,"c6cdom"]~dt[,"cdom"]),col="red")
  qalogical<-readline(message("Press 'Enter' to continue, '1' to set top panel to NA, '2' to set middle panel to NA: ",appendLF=FALSE))
  if(qalogical==1|qalogical=='1'){
    dt[,"cdom"]<-NA
    dtqa[,"cdom"]<-"r"
  }
  if(qalogical==2|qalogical=='2'){
    dt[,"c6cdom"]<-NA
    dtqa[,"c6cdom"]<-"r"
  }
    }
  
  }
  
  if(trimends==TRUE){#NOT IMPLEMENTED YET
  trim<-function(dt){}
  }
  message("QA finished. Printing to file...")
  message(file.path(fdir,"DF_FullDataSets","QA",paste(yearmon,"qa",".csv",sep="")))
  fdir_fd<-file.path(fdir,"DF_FullDataSets","QA")
  write.csv(dtqa,file.path(fdir_fd,paste(yearmon,"qa",".csv",sep="")),row.names = FALSE)
  
}

#'@name streamparse
#'@title Parse old cleaned streaming files
#'@description Includes checks to ensure that data columns are of type numeric. TODO: check that the fathom basins column is populated
#'@param yearmon numeric yyyymm date
#'@param tofile logical save to file?
#'@param fdir character file path to local data directory
#'@export
#'@examples \dontrun{dt<-streamparse(yearmon=201002)}

streamparse<-function(yearmon,tofile=FALSE,fdir=getOption("fdir")){
  #yearmon<-201109
  fdir_fd<-file.path(fdir,"DF_FullDataSets","Raw")
  flist<-list.files(fdir_fd,include.dirs=T,full.names=T)
  flist<-flist[substring(basename(flist),1,6)==yearmon]
  
  dt<-read.csv(flist)
  names(dt)<-tolower(names(dt))
  namestemp<-tolower(names(streamget(201505)))#[-1])
  
  #remove bad coord columns
  coordnames<-c("lat_dd","long_dd","lon_dd")
  for(i in 1:length(coordnames)){
    #i<-1
    cname<-which(!is.na(match(names(dt),coordnames[i])))
     if(length(cname)!=0){
      if(abs(mean(dt[,coordnames[i]]))>100){
        dt<-dt[,-cname]
      }
    }
  }
  
  #remove unrealistic coordinates
  dt <- dt[abs(dt$lat_dd) > 24.5 & abs(dt$lat_dd) < 25.5, ]
  dt <- dt[abs(dt$lon_dd) > 80.1 & abs(dt$lon_dd) < 82, ]
  
    
  #create translation key
  namesalias<-read.table(text="sec,sec.x
cnd,cond
light,trans
fluor,chla",sep=",")
  
  for(n in 1:length(names(dt))){
  #n<-1
    if(any(names(dt)[n]==as.character(namesalias[,1]))){
      names(dt)[n]<-as.character(namesalias[which(names(dt)[n]==namesalias[,1]),2])
    }
  }
  
  #remove non-matching columns
  dt <- dt[,-which(!is.na(match(names(dt), names(dt)[is.na(match(names(dt), namestemp))])))]
  
  #create extra columns if necessary
  dt[,namestemp[is.na(match(namestemp,names(dt)))]]<-NA
  
  #calculate datetime stamp
  #create POSIXct datetime column
  if(mean(nchar(as.character(dt$date)))>6){
    hr<-substring(dt$time,1,nchar(dt$time)-2)
    min<-substring(dt$time,nchar(dt$time)-1,nchar(dt$time))
    
    dt$datetime<-as.POSIXct(strptime(paste(as.character(dt$date)," ",hr,":",min,":",dt$sec.x,sep=""),format="%m/%d/%Y %H:%M:%S"))
    
  }else{
  
  yr<-substring(dt$date,nchar(dt$date)-1,nchar(dt$date))
  day<-substring(dt$date,nchar(dt$date)-3,nchar(dt$date)-2)
  mon<-substring(dt$date,1,nchar(dt$date)-4)
  hr<-substring(dt$time,1,nchar(dt$time)-2)
  min<-substring(dt$time,nchar(dt$time)-1,nchar(dt$time))
  
  if(mean(nchar(mon))==1){mon<-paste("0",mon,sep="")}
  if(!any(!is.na(dt$sec.x))){
    mmin<-Mode(rle(dt$time)$lengths)
    mminseq<-seq(from=0,to=60-60/mmin,by=60/mmin)
    mmin1<-rle(dt$time)$lengths[1]
    mmin1seq<-mminseq[(length(mminseq)-mmin1+1):length(mminseq)]
    dt$sec.x<-c(mmin1seq,rep_len(mminseq,length.out=nrow(dt)-mmin1))
  }
  
  dt$datetime<-paste(yr,"-",mon,"-",day,"-",hr,"-",min,"-",dt$sec.x,sep="")
  #rm(min)
  dt$datetime<-as.POSIXct(strptime(dt$datetime,format="%y-%m-%d-%H-%M-%S"))
  }
  
  #sort columns to match namestemp
  dt<-dt[,match(namestemp,names(dt))]
  
  #ensure that data columns are numeric
  parset<-c("chla","temp","cond","sal","trans","cdom","brighteners","phycoe","phycoc","c6chla","c6cdom","c6turbidity","c6temp")
  dt[,parset]<-suppressWarnings(apply(dt[,parset],2,function(x) as.numeric(x)))
  
  if(tofile==TRUE){
    #add check to verify yearmon before overwriting
    dtname<-file.path(fdir,.Platform$file.sep,"DF_FullDataSets",.Platform$file.sep,substring(basename(flist[1]),1,6),"j.csv",fsep="")
    if(!file.exists(dtname)){
      write.csv(dt,dtname,row.names = FALSE)
    }else{
      stop("overwrite file?")
    }
  }else{
    dt
  }
}

