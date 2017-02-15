### R code from vignette source 'DataflowR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: DataflowR.Rnw:65-67 (eval = FALSE)
###################################################
## install.packages(path.to.zip, type = "win.binary", repos = NULL,
##                  dependencies = TRUE)


###################################################
### code chunk number 2: DataflowR.Rnw:76-78 (eval = FALSE)
###################################################
## install.packages("devtools")
## devtools::install_github("jsta/DataflowR")


###################################################
### code chunk number 3: DataflowR.Rnw:103-104 (eval = FALSE)
###################################################
## system.file("localpath", package = "DataflowR")


###################################################
### code chunk number 4: DataflowR.Rnw:150-152 (eval = FALSE)
###################################################
## dt <- streamclean(yearmon = 201606, gps = "eu", eummin = 12, c6mmin = 12,
##       tofile = FALSE)


###################################################
### code chunk number 5: DataflowR.Rnw:157-158 (eval = FALSE)
###################################################
## dt <- streamparse(yearmon = 201007, tofile = FALSE)


###################################################
### code chunk number 6: DataflowR.Rnw:168-169 (eval = FALSE)
###################################################
## streamqa(yearmon = 201606, parset = names(streamget(201606))[c(4:12, 16:22)])


###################################################
### code chunk number 7: DataflowR.Rnw:178-179 (eval = FALSE)
###################################################
## dt <- streamget(yearmon = 201606, qa = TRUE)


###################################################
### code chunk number 8: DataflowR.Rnw:190-192 (eval = FALSE)
###################################################
## streaminterp(streamget(yearmon = 201606, qa = TRUE),
##   paramlist = c("salinity.pss"), 201606)


###################################################
### code chunk number 9: DataflowR.Rnw:201-202 (eval = FALSE)
###################################################
## surfplot(rnge = c(201502), params = c("sal"))


###################################################
### code chunk number 10: DataflowR.Rnw:218-219 (eval = FALSE)
###################################################
## grassmap(rnge = 201505, params = c("sal"))


###################################################
### code chunk number 11: DataflowR.Rnw:226-228 (eval = FALSE)
###################################################
## grassmap(rnge = c(201205, 201305), params = c("sal"),
##          basin = "Manatee Bay", numcol = 3, numrow = 3)


###################################################
### code chunk number 12: DataflowR.Rnw:236-237 (eval = FALSE)
###################################################
## grabclean(yearmon = 201410, tofile = FALSE)


###################################################
### code chunk number 13: DataflowR.Rnw:246-247 (eval = FALSE)
###################################################
## grabs <- grabget(rnge = c(201402, 201410))


###################################################
### code chunk number 14: DataflowR.Rnw:254-256 (eval = FALSE)
###################################################
## avmap(yearmon = 201502, params = "sal", tofile = TRUE, percentcov = 0.6,
##       tolerance = 1)


###################################################
### code chunk number 15: DataflowR.Rnw:274-275 (eval = FALSE)
###################################################
## chlcoef(yearmon = 201502, remove.flags = TRUE)


###################################################
### code chunk number 16: DataflowR.Rnw:282-283 (eval = FALSE)
###################################################
## chlmap(yearmon = 201502)


