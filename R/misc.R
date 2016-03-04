# @name fdir_asgn
# @title load path to local Dataflow directory
# @export
# fdir_asgn<-function(){
# ldir<<-matrix(utils::read.delim(system.file("localpath",package="DataflowR"),header=FALSE,stringsAsFactors=FALSE))[[1]][1]
# dir.create(ldir)
# }

#'@name Mode
#'@title Mode
#'@description Returns the mode of a numeric array
#'@export
#'@param x numeric array
Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#'@name cond2sal
#'@title Conductivity to salinity conversion
#'@description Uses the PSS-78 practical salinity equation
#'@export
#'@param c numeric conducitivity in uS (measurements in mS need to multiplied by 1000)
#'@param t numeric celcius temperature
#'@param P numeric optional pressure (defaults to 0)
#'@examples \dontrun{csal<-cond2sal(c=dt$cond*1000,t=dt$temp)
#'plot(csal, dt$sal)
#'abline(a=0,b=1,col="red",lwd=2)}
#'@details algorithm based off an excel implementation by N. Iricanin
#'@seealso \code{\link[wq]{ec2pss}}

cond2sal<-function (c, t = 25, P = 0) 
{
  a0 = 0.008
  a1 = -0.1692
  a2 = 25.3851
  a3 = 14.0941
  a4 = -7.0261
  a5 = 2.7081
  b0 = 5e-04
  b1 = -0.0056
  b2 = -0.0066
  b3 = -0.0375
  b4 = 0.0636
  b5 = -0.0144
  c0 = 0.6766097
  c1 = 0.0200564
  c2 = 0.0001104
  c3 = -6.9698e-07
  c4 = 1.0031e-09
  D1 = 0.03426
  D2 = 0.0004464
  D3 = 0.4215
  D4 = -0.003107
  e1 = 0.000207
  e2 = -6.37e-08
  e3 = 3.989e-12
  Csw = 42.914
  K = 0.0162
  Ct = round(c * (1 + 0.0191 * (t - 25)), 0)
  R = (Ct/1000)/Csw
  rt = c0 + (t * c1) + (t^2 * c2) + (t^3 * c3) + (t^4 * c4)
  Rp = 1 + (P * e1 + e2 * P^2 + e3 * P^3)/(1 + D1 * t + D2 * 
                                             t^2 + (D3 + D4 * t) * R)
  Rt1 = R/(Rp * rt)
  dS = (b0 + b1 * Rt1^(1/2) + b2 * Rt1^(2/2) + b3 * Rt1^(3/2) + 
          b4 * Rt1^(4/2) + b5 * Rt1^(5/2)) * (t - 15)/(1 + K * 
                                                         (t - 15))
  S = a0 + a1 * Rt1^(1/2) + a2 * Rt1^(2/2) + a3 * Rt1^(3/2) + 
    a4 * Rt1^(4/2) + a5 * Rt1^(5/2) + dS
  
    S[is.na(S<0)]<-NA
  
    S[S<2 & !is.na(S)]<- S[S<2 & !is.na(S)] - a0/(1 + 1.5 * (400 * Rt1) + (400 * Rt1)^2) - 
        (b0 * (t - 15)/(1 + K * (t - 15)))/(1 + (100 * Rt1)^(1/2) + 
                                              (100 * Rt1)^(3/2))
    PSS = round(S, 3)
  
  return(PSS)
}

#'@name date456posix
#'@title Convert numeric dates in mddyy to POSIXct
#'@description Convert numeric dates in mddyy to POSIXct 
#'@param x numeric where the first 1-2 digits specify the month attribute because leading zeros have been stripped. Also detects whether the day attributes has had stripped leading zeros.
#'@param century numeric century recommended choice of "19" or "20"
#'@export
#'@examples
#'x<-51514
#'date456posix(x,century="20")
#'dates<-c("51514","101214","8714","1214")
#'date456posix(dates,century="20")

date456posix <- function(x, century){
  year <- paste0(century, substring(x, (nchar(x) - 1), nchar(x)))
  day <- substring(x, (nchar(x) - 3), nchar(x) - 2)
  
  if(any(as.numeric(day)>31)){
    day<-as.character(sapply(day,function(x){
      if(as.numeric(x)>31){
        x<-substring(x,2,2)
      }
      x
    }))
  }
  
  
  mon<-substring(x,1,nchar(x)-4)
  
  if(any(nchar(mon)==0)){
    mon[which(nchar(mon)==0)]<-substring(x[which(nchar(mon)==0)],1,1)
  }
  
  mon<-as.character(sapply(mon,function(x){
    if(as.numeric(x)<10){
      x<-paste0("0",x)
    }else{
      x
    }
  }))
  
  date<-paste0(year,"-",mon,"-",day)
  return(as.POSIXct(strptime(date,format="%Y-%m-%d")))
}


#'@name coordinatize
#'@title Convert georeferenced data.frames into projected SpatialPointsDataFrames
#'@description Convert georeferenced data.frames into projected SpatialPointsDataFrames
#'@param latname character column name of the "y" coordinate
#'@param lonname character column name of the "x" coordinate
#'@param dt data.frame with two coordinate columns
#'@export
#'@examples \dontrun{
#'dt<-coordinatize(streamget(201002), latname = "lat_dd", lonname = "lon_dd")
#'}
coordinatize <- function(dt, latname = "latdec", lonname = "londec"){
  projstr<-"+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  latlonproj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  dt<-dt[!is.na(dt[,lonname]),]
  dt<-dt[!is.na(dt[,latname]),]
  
  sp::coordinates(dt)<-c(lonname,latname)
  sp::proj4string(dt)<-sp::CRS(latlonproj)
  dt<-sp::spTransform(dt,sp::CRS(projstr))
  return(dt)
}

#'@name logramp
#'@title Create a log scaled color ramp for the grassmap function
#'@description Create a log scaled color ramp for the grassmap function
#'@export
#'@import viridis
#'@import scales 
#'@param n integer number of breaks
#'@param maxrange integer maximum range
#'n <- 9
#'maxrange <- 20
#'logramp(n, maxrange)
logramp <- function(n, maxrange){
  breaks <- round(exp((seq(from=log(1),to=log(maxrange + 1),length=n)))-1, 2)
  displaybreaks <- round(exp((seq(from=log(1),to=log(breaks[length(breaks)-1]),length=6)))-1, 2)
  
  #'scales::show_col(viridis::viridis_pal()(9))
  
  list(breaks = breaks, displaybreaks = displaybreaks, rgb = t(col2rgb(viridis::viridis(n))))
}

#'@name mdy2mmyyyy
#'@title convert m/d/yy to mm/dd/yyyy
#'@description Pads dates in preparation for POSIX coercion
#'@param x character date to be formatted
#'@export
#'@examples
#' x <- "5/5/15"
#' mdy2mmyyyy(x)
mdy2mmyyyy <- function(x){
  
  #strsplit based on "/"
  month <- strsplit(x, "/")[[1]][1]
  if(nchar(month) < 2){
    month <- paste("0", month, sep="")
  }
  day <- strsplit(x, "/")[[1]][2]
  if(nchar(day) < 2){
    day <- paste("0", day ,sep="")
  }
  year <- strsplit(x, "/")[[1]][3]
  if(nchar(year) < 3 & as.numeric(year) < 80){
    year <- paste("20", year, sep="")
  }else{
    if(nchar(year) < 3){
      year <- paste("19", year, sep="")
    }
  }
  
  paste(month, "/", day, "/", year, sep="")
}
