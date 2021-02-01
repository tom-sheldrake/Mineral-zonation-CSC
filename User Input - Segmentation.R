########## 0. PRELIMINARY CHECKS AND INFORMATION ##########

# i: Ensure the this file is in the CSC folder
# ii: Ensure that a raster file with the correct name (i.e. the sample name) is contained in "~/Phase classification/{Sample name}/"
# iii: This file will be automatically in the correct folder if youu have used the Phase classification algorithm

### make sure the following packages are installed ###
# install.packages("raster",dependencies = TRUE)
# install.packages("conicfit",dependencies = TRUE)
# install.packages("viridis",dependencies = TRUE)
# install.packages("imager",dependencies = TRUE)
# install.packages("apcluster",dependencies = TRUE)

### Instructions to run this code ###
#Once the user input in section 1 has been finalised, highlight sections 1 & 2, then press run.
#Wait for the relevant crystals to be extracted, update section 3, then highlight sections 3 & 4, then press run.

### Summary of output ###
#The output of the segmentation is saved as an Rdata file, which can be reopened and is used for the correlation of crystal zones.
#The Rdata file will contain all the parameters that have been defined by the user.
#IMPORTANT: Each time the code is re-run, you will overwrite the previous results in "~/Crystal segmentation/{Sample name}/"


########### 1. USER INPUT REQUIRED  - DEFINES THIN SECTION AND PARAMETERS FOR SEGMENTATION ########### 

### Choose the sample ###
Sample_name <- "SK408"

### Pixel size equivalent (in um) ###
pix.mm <- 20

### Minimum crystal size to segment (in um) ###
min.area <- 180^2

### Choose which number phase to analyse, based on the phase map produced by the phase classification ###
### It can be multiple numbers, if individual phases have been classified into multiple groups ###
Phase <- c(1)

### Choose which elements to perform the segmentation on. This could also be BSE, for example ###
Element.choices <- c("Ca","Na")

### Set fraction limit for crystal pixels to use in AP ###
xstal.pc <- 0.9

### Set maximum percentile of diagonal in AP matrix (default is 0.05 -- 5th percentile) ###
percentile.max <- 0.05

### Set number of iterations the SLICAP algorithm is run for convergence ###
N.Iters <- 10

### Set search radius w.r.t. pixel radius. Leave as NA to use automatic setting (*2) ###
search.radius.manual <- NA

### Set spacing of interval. Leave as NA to use automatic setting ###
S.manual <- NA

### Set chemical weighting parameter. Leave as NA to use automatic setting ###
M.manual <- NA

### Set home dir (leave this as default if using RStudio) ###
home.dir <- paste(dirname(rstudioapi::getActiveDocumentContext()$path))

########### 2. LOADS IN THIN SECTION AND EXTRACTS CRYSTALS TO SEGMENT ########### 

# loads packages
require(raster)
require(conicfit)
require(viridis)
library(imager)
library(apcluster)

# Define user-derived functions
'%!in%' <- function(x,y)!('%in%'(x,y))
range01 <- function(x){(x-min(x,na.rm = TRUE))/(max(x,na.rm = TRUE)-min(x,na.rm = TRUE))}
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
score.funct <- function(N,arr,r1,c1,r2,c2) {((arr[r1,c1,N] - arr[r2,c2,N])^2)}
med.dk <- function(x,p) {median(x[p],na.rm = TRUE)}

setwd(home.dir)
setwd("Crystal segmentation/")
source("0-Crystal_extraction.R")

plot(Phase.segment,col=grey.colors(2),main=Sample_name,legend=F)
invisible(lapply(1:length(xstal.to.segment),function(x) {lines(Exterior[[xstal.to.segment[x]]],lty=1,lwd=1.5,col="red")}))
text(labpt[1:length(xstal.to.segment),],labels=1:length(xstal.to.segment),col="red",cex=1, font = 2)

########### 3. SETS WHICH CRYSTALS TO SEGMENT, AND EACH CRYSTAL ARE PLOTTED ########### 

### Defines crystals to segment (if all crystals, leave as NA, otherwise choose crystal numbers ) ###
xstal.choice <- NA    #xstal.choice <- c(2)

### Plot individual crystals (TRUE or FALSE) ###
plot.xstal <- FALSE 

########### 4. PERFORM SLICAP ALGORITHM AND PLOTS SEGMENTED RESULTS ########### 

#Make rasters to store the clusters and labels in 
assign("textures", raster(ras@extent, ncols = ras@extent[2], nrows = ras@extent[4]))
textures[] <- 0 
assign("label.ras", raster(ras@extent, ncols = ras@extent[2], nrows = ras@extent[4]))

xstal.to.segment <- if(is.na(xstal.choice)==TRUE) {xstal.to.segment} else {xstal.choice}

for(N in 1:length(xstal.to.segment)){

  N <- xstal.to.segment[N]
  cat(paste("\nCrystal",N,"\n"))
  
  setwd(home.dir)
  setwd("Crystal segmentation/")
  
  source("1-SLICPreprocess.R")
  source("2-SLIC.R")
  source("3-SLICPostprocess.R")
  source("4-AP.R")

  if(plot.xstal==TRUE) {source("4b-Crystal_plots.R")}
  
  setwd(home.dir)
  setwd("Crystal segmentation/")
  source("5-Output.R")
  
  rm(Ap,Ck,Cl,Ck.grid,Ck.initial,coords,coords.temp,d,Ellipse_geometric,Gradient,id.xstal,inv.xstal,l,sps,superpixel.id)
  rm(list=ls(pattern = "mat_"))
  rm(list=ls(pattern = "centroids.temp"))
  rm(bg,exterior,M,mat,mat.master,N,n.clus,N.mat,ncols,nrows,q,q1,q2,q3,S,search.radius,xmax,xmin,xpix.ts,ymax,ymin,ypix.ts)
}

textures[which(textures[]==0)] <- NA
names(textures) <- "Zones"
textures$Labels <- label.ras$layer
rm(label.ras)

labels.id <- unique(textures$Labels[])
textures$Labels[] <- match(textures$Labels[],labels.id)

Segmented <- stack(Phase.segment,textures,ras)
rm(Phase.segment,textures,ras)

setwd(home.dir)
setwd("Crystal segmentation/")
setwd(Sample_name)

writeRaster(x=Segmented, filename = paste0(Sample_name,"_segmented.grd"),overwrite=TRUE,type="")
plot(Segmented$Zones,col=magma(max(Segmented$Zones[],na.rm=TRUE)))
rm(Segmented,home.dir)

save(list = ls(),file = paste0(Sample_name,"_Zones.Rdata"))






