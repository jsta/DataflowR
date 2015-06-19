.onLoad<-function(libname=find.package("DataflowR"),pkgname="DataflowR"){
  options("fdir"=matrix(utils::read.delim(system.file("localpath",package="DataflowR"),header=FALSE,stringsAsFactors=FALSE))[[1]][1])
  
  #fdir_asgn()
  #assign("pkg_globals",new.env(),envir=parent.env(environment()))
  
  
#   ldir<<-matrix(utils::read.delim(system.file("localpath",package="DataflowR"),header=FALSE,stringsAsFactors=FALSE))[[1]][1]
#   dir.create(ldir)
}




