####### PREPARES DATA FOR SLICAP ALGORITHM FOR INDIVIDUAL CRYSTAL ####### 

### 1a. Extract individual crystal as smaller raster ###

#Specify dimensions of sample
xpix.ts <- ras@nrows
ypix.ts <- ras@ncols

#Select the crystal to do the image segmentation on
poly.temp <- crystals@Polygons[N]
#Extract the coordinates of the border of the selected pixel
coords.temp <- poly.temp[[1]]@coords

#Fits an ellipse to exterior of the crystal
Ellipse_algebraic <- EllipseDirectFit(coords.temp)
Ellipse_geometric <- AtoG(Ellipse_algebraic)$ParG

#Extract the x and y limits of a square plotting region around the crystal
xmin <- min(coords.temp[,1]) ; xmax <- max(coords.temp[,1])
ymin <- min(coords.temp[,2]) ; ymax <- max(coords.temp[,2])

#crop to extent of individual crystal 
e <- extent(x = c(xmin,xmax),y=c(ymin,ymax))
inv.xstal <- crop(ras,e)

#Specify dimensions of individual crystal
ncols <- inv.xstal@ncols
nrows <- inv.xstal@nrows

#plots chemistry of the individual crystal (if user demands)
if(plot.xstal==TRUE) {
setwd(Sample_name)
wd <- ncols*7
ht <- nrows*7
if(wd<ht) {ht <- (ht/wd)*1000 ; wd <- 1000} else {wd <- (wd/ht)*1000 ; ht <- 1000}
plot.rows <- round(length(Element.choices)^0.5,0)
plot.cols <- ceiling(length(Element.choices)^0.5)
png(filename = paste0("Crystal",N,"_elements.png"),width = plot.cols*wd,height = plot.rows*ht,units = "px")
par(mfrow=c(plot.rows,plot.cols))
lapply(Element.choices, function(x) {dev.hold(level = 1L); plot(eval(parse(text = paste0("inv.xstal$",x))),main=x,col=plasma(100)) ; lines(coords.temp[,"x"], coords.temp[,"y"], type = "l",lwd=5,col="white")})
dev.off()
rm(plot.rows,plot.cols)
setwd(home.dir)
setwd("Crystal segmentation/")
}

### 1b. Identify crystal and non-crystal mask ###

#Create a dataframe that contains all of the permetations of coordinates that may be inside or outside the crystal
x.mat <- seq(from = xmin+0.5, to = xmax-0.5, by = 1)
y.mat <- seq(from = ymin+0.5, to = ymax-0.5, by = 1)
gr <- as.data.frame(expand.grid(x = x.mat, y = y.mat))

#Use the point.in.polygon function to check if the coordinates are inside or not
gr$pol <- point.in.polygon(point.x = gr$x, point.y = gr$y, pol.x = coords.temp[, "x"], pol.y = coords.temp[, "y"])
exterior <- which(gr$pol == 0)

#Convert exterior points to a matrix (i.e. rotated)
mat.gr <- gr
mat.gr$x <- ymax - gr$y + 0.5
mat.gr$y <- gr$x - xmin + 0.5

#Flip raster to set exterior and then flip back
ras.temp <- flip(inv.xstal,"y")

#Defines non crystal mask, and sets all values to -99
non.xstal <- which(ras.temp$Phases@data@values!=Phase | is.na(ras.temp$Phases@data@values)==TRUE) 
non.xstal.mask <- unique(c(exterior,non.xstal))
ras.temp@data@values[non.xstal.mask,] <- -99

#Re-assigns new raster to original raster for indiviudal crystal
inv.xstal <- flip(ras.temp,"y")


### 1c. Redefines the chemical content as integers (this is required for the AP algorithm) ###
#contd. Assigns two matrices for the SLIC and AP algorithms ###

#Number of elements used for segmentation
N.mat <- length(Element.choices)

#Extracts matrices for chosen elements for segmentation
lapply(Element.choices, function(x) {assign(paste0("mat_",x),matrix(eval(parse(text = paste0("inv.xstal$",x))),nrow = nrows,ncol = ncols,byrow = TRUE),envir = .GlobalEnv)})

#Identify which values are in the non-xstal mask, and sets values to NA
bg <- which(get(paste0("mat_",Element.choices[1]))==-99,arr.ind = TRUE)
for (i in 1:N.mat){
  mat.temp <- get(ls(pattern = "mat_",pos = sys.frame())[i])
  mat.temp[bg] <- NA
  assign(ls(pattern = "mat_",pos = sys.frame())[i],mat.temp,pos = sys.frame())
  rm(mat.temp)
}

#Normalises values of chosen elements between 0 and 1, and inputs all matrices into a single array
lapply(ls(pattern = "mat_",envir = sys.frame()), function(x) {assign(x,range01(get(x,pos = sys.frame())),envir = sys.frame())})
mat.master <- array(unlist(lapply(ls(pattern = "mat_",pos = sys.frame()),get,pos = sys.frame())),dim = c(nrows,ncols,N.mat))

#Converts NA values in matrices back to -99, which is needed for the SLIC algorithm
for (i in 1:N.mat){
  mat.temp <- get(ls(pattern = "mat_",pos = sys.frame())[i])
  mat.temp[bg] <- -99
  assign(ls(pattern = "mat_",pos = sys.frame())[i],mat.temp,pos = sys.frame())
  rm(mat.temp)
}
#Inputs all numeric matrices into a single array
mat <- array(unlist(lapply(ls(pattern = "mat_",pos = sys.frame()),get,pos = sys.frame())),dim = c(nrows,ncols,N.mat))


### 1d. Caclulates gradient between pixels to help initiate superpixel centroids ###

#Table of all coordinates in individual crystal raster
coords <- expand.grid(seq(1,nrows,1),seq(1,ncols,1))

#Function to calculate gradient
grad.func <- function(x,M) {
  i <- coords[x,1] ; j <- coords[x,2]
  xr <- seq(i-1,i+1,1) ; xr <- xr[which(xr > 0 & xr <= nrow(M))]
  xc <- seq(j-1,j+1,1) ; xc <- xc[which(xc > 0 & xc <= ncol(M))]
  mean.val <- if (N.mat==1) {mean(mat[xr,xc,1],na.rm=TRUE)} else {apply(mat[xr,xc,],MARGIN = 3,mean,na.rm=TRUE)}
  return((mat[i,j,] - mean.val)^2)
}

#Calculates gradient for all coordinates, and sets all background pixels to 99
grad <- lapply(1:nrow(coords), grad.func,M=mat)
grad <- unlist((lapply(grad,mean)))^0.5
grad <- matrix(as.numeric(grad),nrow = nrows,byrow = FALSE)
grad[which(is.na(grad)==TRUE,arr.ind = TRUE)] <- 1
grad[bg] <- 99

#Converts gradient to a raster 
Gradient <- raster(grad) ; extent(Gradient) <- extent(inv.xstal)

### 1e. Defines position of centorids, and weighting between spatial and chemical distances ###

#Calculate optimum interval between initial centroids based on the 2^log10(area of each image)
#Hence it is linear in log-log space: plot(c(100,1000,10000,100000),c(2,4,8,16),log="xy")
S <- round(2^(log(nrows*ncols,10)-1),0)

#Overwrite S in manual spacing is chosen
S <- ifelse(is.na(S.manual)==TRUE,max(S,3),S.manual)

#Calculate the weighting factor for the [x,y] distance in the scoring function D.
#Ratio of S to mean pixels or x & y dimension
M <- (S/mean(c(nrows,ncols)))

#Overwrite M in manual weighting is chosen
M <- ifelse(is.na(M.manual)==TRUE,M,M.manual)

#Defines the x & y values for initial superpixel centorids
nr <- seq(0,(floor((ncols-2)/S))*S,by = S); nr <- nr + floor((ncols-max(nr))/2)
nc <- seq(0,(floor((nrows-2)/S))*S,by = S); nc <- nc + floor((nrows-max(nc))/2)
Ck <- expand.grid(nc,nr) ; Ck <- matrix(unlist(Ck),ncol=2,byrow = FALSE)
Ck.grid <- Ck

#Function to find lowest gradient position around an initial centroid
Ck.lowest.grad <- function(x,M){
  i <- Ck[x,1] ; j <- Ck[x,2]
  xr <- seq(i-1,i+1,1) ; xr <- xr[which(xr > 0 & xr <= nrow(M))]
  xc <- seq(j-1,j+1,1) ; xc <- xc[which(xc > 0 & xc <= ncol(M))]
  if (grad[i,j] > min(grad[xr,xc])) {
  return(as.numeric(expand.grid(xr,xc)[which.min(grad[xr,xc]),]))
  } else
  return(Ck[x,])
}

#Move centroids to lowest gradient position
Ck <- matrix(unlist(lapply(1:nrow(Ck), Ck.lowest.grad,M=mat)),ncol=2,byrow=TRUE)
Ck <- unique.matrix(Ck)
Ck.initial <- Ck


rm(e,Ellipse_algebraic,gr,grad,mat.gr,poly.temp,ras.temp,
   i,nc,non.xstal,non.xstal.mask,nr,x.mat,y.mat)

