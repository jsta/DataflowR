#'@name chlcoef
#'@title calculation of extracted versus fluoresced chlorophyll coefficients
#'@param yearmon numeric date of survey
#'@param remove.flags logical trim dataset based on flags?
#'@details this function should be interactive
#'@export
#'@importFrom MASS stepAIC
#'@importFrom car vif

chlcoef<-function(yearmon,remove.flags=FALSE){

#suppressMessages(library(MASS))
#library(car)

dt<-grabget(yearmon)

if(remove.flags==TRUE){
  dt<-dt[dt$flags!="s",]
}

#remove all or partial incomplete cases?
#dt2<-dt[complete.cases(dt),]#all incomplete
#dt2<-dt[!is.na(dt$chla),]#partial incomplete

#choose variables
varlist<-c("chla","cdom","chlaiv","phycoe","c6chl","phycoc","c6cdom")
varlist<-varlist[sapply(varlist,function(x) any(!is.na(dt[,x])))]

cormat<-cor(dt[,varlist],use="complete")[-1,1]
varlist<-names(cormat[abs(cormat)>0.4])
lmeq<-as.formula(paste("chla ~ ", paste(varlist,collapse="+")))
fit<-lm(lmeq,data=dt)

if(summary(fit)$adj.r.squared<0.7){#poly fit
  lmeq<-as.formula(paste("chla ~ ", paste("poly(",varlist,",2,raw=T)",collapse="+")))
  fit<-lm(lmeq,data=dt[complete.cases(dt[,varlist]),c("chla",varlist)])
  poly=TRUE
}

if(summary(fit)$r.squared<0.6){
  stop("bad fit. r-squared not high enough.")
}

saic<-MASS::stepAIC(fit)#pick reduced eq according to maximized AIC (remove terms with a the smallest (largest negative score))
rmlist<-as.character(saic$anova$Step)
rmlist<-gsub("-","",rmlist)
rmlist<-gsub(" ","",rmlist)
rmlist<-rmlist[nchar(rmlist)>1]
rmlist<-gsub("poly\\(","",rmlist)
rmlist<-gsub(",2,raw=T\\)","",rmlist)

varlist<-varlist[is.na(match(varlist,rmlist))]
lmeq<-as.formula(paste("chla ~ ", paste(varlist,collapse="+")))
fit<-lm(lmeq,data=dt)

#examine VIF; If GVIF > 5:10, then remove colinear terms, Helsel and Hirsh 2002
if(length(fit$coefficients)>2){
viftest<-car::vif(fit)
if(length(which(viftest>10))>0){
  rmlist<-names(viftest)[which(viftest>5)]
  varlist<-varlist[is.na(match(varlist,rmlist))]
  lmeq<-as.formula(paste("chla ~ ", paste(varlist,collapse="+")))
  if(poly==TRUE){
    lmeq<-as.formula(paste("chla ~ ", paste("poly(",varlist,",2,raw=T)",collapse="+")))
    }
  fit<-lm(lmeq,data=dt[complete.cases(dt[,varlist]),c("chla",varlist)])
}
}

plot(fitted(fit),dt$chla[complete.cases(dt[,varlist])])
abline(0,1)

# #optimize for low values (<1)?####
# #hist(dt2$CHLa)
# dt2<-dt2[dt2$CHLa<2.5,]
# #full eq
# fit<-lm(CHLa~poly(CDOM,2,raw=T)+poly(CHLAiv,2,raw=T)+poly(Phycoerythrin,1,raw=T)+poly(C6Chla,2,raw=T)+poly(C6cdom,2,raw=T),data=dt3)
# #DF only eq
# fit<-lm(CHLa~poly(CDOM,2,raw=T)+poly(CHLA,2,raw=T),data=dt3)
# stepAIC(fit)
# fit2<-lm(CHLa~poly(CDOM,2,raw=T),data=dt3)
# plot(predict(fit2,dt2),dt2$CHLa,xlim=range(dt2$CHLa))#predict based on full dataset 
# plot(predict(fit2,dt3),dt3$CHLa,xlim=range(dt3$CHLa))#predict based on reduced dataset 
# 
# fit2<-lm(CHLa~poly(C6Chla,2,raw=T),data=dt3)
# plot(fitted(fit2),dt3$CHLa)#inspect 1:1
# abline(a=0,b=1)
# #fit to full dataset
# plot(predict(fit2,dt2),dt2$CHLa)#inspect 1:1
# abline(a=0,b=1)
# #####

#retrieve chla coefficients####

coeflist<-read.csv(file.path(fdir,"DF_GrabSamples","extractChlcoef2.csv"),header=T,na.strings="NA")[,-1]
names(coeflist)<-tolower(names(coeflist))
vartemplate<-c("yearmon","survey","date","cdom","chlaiv","phycoe","c6chl","c6cdom","phycoc","cdom2","chlaiv2","phycoe2","c6chl2","c6cdom2","phycoc2","intercept","rsquared","pvalue","model","notes")
outtemp<-data.frame(matrix(NA,nrow=1,ncol=length(vartemplate)))
names(outtemp)<-vartemplate

outcoef<-fit$coefficients
names(outcoef)[1]<-"intercept"
outtemp[match(names(outcoef),names(outtemp))]<-outcoef
outtemp[,"yearmon"]<-yearmon
outtemp[,"date"]<-Mode(dt[,"date"])
outtemp[,"rsquared"]<-summary(fit)$r.squared

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

outtemp[,"pvalue"]<-lmp(fit)

model<-outtemp[4:16]
model[1,]<-round(as.numeric(model[1,]),5)
model<-data.frame(matrix(c(model,names(model)),nrow=2,byrow=T))
model<-model[,!is.na(model[1,])]
model[1,]<-sapply(model[1,],as.character)
intercept<-model[1,ncol(model)]
model<-model[,-ncol(model)]
if(length(model)<3){
model<-data.frame(matrix(unlist(model),nrow=2))#new
}

outtemp[,"model"]<-paste("Chla = ",gsub(" ","+",gsub(",","",toString(apply(model,2,function(x) paste(x[1],x[2],sep="*"))))),"+",intercept,sep="")

if(any(outtemp[,"yearmon"]==coeflist[,"yearmon"],na.rm=T)){
  stop("replace previously recorded fit?")
}else{
  coeflist<-rbind(coeflist,outtemp)
  write.csv(coeflist,file.path(fdir,"DF_GrabSamples","extractChlcoef2.csv"))
}
}
      



#output extracted chl raster
#list.files(paste("DF_Surfaces/",yearmon,"/",sep=""))
# test<-raster(paste("DF_Surfaces/",yearmon,"/","c6chla",".tif",sep=""))
# test2<-raster(paste("DF_Surfaces/",yearmon,"/","c6cdom",".tif",sep=""))
# test<-reclassify(test,c(-Inf,0,NA))
# test2<-reclassify(test2,c(-Inf,0,NA))
# plot(test2)
# 
# testfinal<-test2*(-0.004969)+test*(0.041821)-0.37868
# testfinal<-reclassify(testfinal,c(-Inf,0,0))
# writeRaster(testfinal,filename=paste("DF_Surfaces/",yearmon,"/",yearmon,"chlext",".tif",sep=""),overwrite=T,format="GTiff")
# 
# 
# #calculate extracted chlorophyll####
# # chldt<-merge(fulldataset.over,coeflist,by="SURVEY")
# # chldt$one<-chldt$CDOM_RFU.y*chldt$CDOM_RFU.x
# # chldt$two<-chldt$FLUOR_CHLA_RFU.y*chldt$FLUOR_CHLA_RFU.x
# # chldt$three<-chldt$PHYCOE_RFU.y*chldt$PHYCOE_RFU.x
# # chldt$four<-chldt$C6CHLA_RFU.y*chldt$C6CHLA_RFU.x
# # chldt$five<-chldt$C6CDOM_RFU.y*chldt$C6CDOM_RFU.x
# # chldt$six<-chldt$CDOM_RFU2*chldt$CDOM_RFU.x*chldt$CDOM_RFU.x
# # chldt$seven<-chldt$FLUOR_CHLA_RFU2*chldt$FLUOR_CHLA_RFU.x*chldt$FLUOR_CHLA_RFU.x
# # chldt$eight<-chldt$PHYCOE_RFU2*chldt$PHYCOE_RFU.x*chldt$PHYCOE_RFU.x
# # chldt$nine<-chldt$C6CHLA_RFU2*chldt$C6CHLA_RFU.x*chldt$C6CHLA_RFU.x
# # chldt$ten<-chldt$C6CDOM_RFU2*chldt$C6CDOM_RFU.x*chldt$C6CDOM_RFU.x
# }