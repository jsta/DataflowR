# Intro

The DataflowR package is used interally at SFWMD to handle streaming output and discrete grab samples collected as part of the Dataflow monitoring program.

# Installation

## Prereqs

* dbhydroR
* R
* RStudio
* RTools (Windows only)

## Commands

`install.packages(c("rasterVis", "raster", "maptools", "rgeos", "rgrass7",`
`"scales", "viridis", "MASS", "car", "gdata", "ipdw", "latticeExtra", "rgdal", "sp", "zoo"))`

`install.packages(c("devtools", "getPass"))`

`devtools::install_git("http://gitlab.com/jsta/DataflowR.git",`
`  credentials = git2r::cred_user_pass("<username>",`             `  getPass::getPass()))`

# Function Documentation and Examples

see either:

* vignettes/DataflowR.pdf

* inst/doc/DataflowR.pdf

* `??DataflowR`
