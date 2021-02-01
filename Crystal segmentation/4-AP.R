####### AFFINITY PROPAGATION ALGORITHM #######

### 4a. Calculation of dissimilarity matrix ###

#Convert matrix of superpixels (l) to raster (sps) and define extent as same as individual crystal
sps <- raster(l)
extent(sps) <- extent(inv.xstal)

#Define list of superpixel id's
superpixel.id <- sort(unique(as.vector(l)))

#This calculates: (i) the mean chemical value of a superpixel; and (ii) the percentage of pixels belonging to mask 2
xstal.percent <- NULL
centroid <- matrix(,ncol=N.mat,nrow = length(superpixel.id))
for(i in 1:length(superpixel.id)){
  id.temp <- which(l==superpixel.id[i],arr.ind = TRUE)
  val.temp <- mat.master[cbind(id.temp,rep(1,nrow(id.temp)))]
  cen<- if (N.mat==1) {median(val.temp,na.rm=TRUE)} else {apply(mat.master, 3, med.dk,id.temp)}
  if (N.mat==1) {centroid[i] <- cen} else {centroid[i,] <- cen}
  xstal.percent[i] <- length(which(is.na(val.temp)==TRUE))/length(val.temp)
} ; rm(i,cen,id.temp)

#Rounds central values between 0 and 100, to the nearest integer
Cen_vals <- round(100*centroid,0)

#Identifies crystals dominated by non-xstal pixels and removes them from the AP
id.xstal <- which(xstal.percent < xstal.pc)
Cen_vals <- Cen_vals[id.xstal,] 

#This calculates the similarity matrix using the centroids from all superpixels
X <- if (N.mat==1) {length(Cen_vals)} else {nrow(Cen_vals)}
diss <- matrix(,nrow=X,ncol=X)
for (i in 1:X){
  for (j in 1:X){
    cen1 <- if (N.mat==1) {Cen_vals[i]} else {Cen_vals[i,]}
    cen2 <- if (N.mat==1) {Cen_vals[j]} else {Cen_vals[j,]}
    cen.dif <- -1*(sum((cen1 - cen2)^2))
    diss[i,j] <- cen.dif
  }
} ; rm(i,j,X)

### 4b. Perform AP algorithm ###

# Defines percentile for diagonal of dissimilarity matrix
q1 <- (max(Ellipse_geometric[c(3:4)]))/(min(Ellipse_geometric[c(3:4)])) #aspect ratio
q2 <- log(areas[1] / areas[N]) #area realtive to all other crystals
q3 <- log(percentile.max)+2 
q <- exp(q3-q1-q2)

# Performs the AP algorithm, and estimates number of clusters
Ap <- apcluster(diss,q=q) 
n.clus <- length(Ap@exemplars)

#Defines a vector Cl, which we will fill up with results of AP algorithm
Cl <- as.vector(l)

#Identfies background pixels - i.e. with majority non-xstal pixels
bg <-  which(superpixel.id %!in% id.xstal)
Cl[which(l %in% bg)] <- NA

#Calculates which pixel belongs to which cluster
Cl.master <- Cl
for (i in 1:n.clus){
  Cl[which(Cl.master %in% id.xstal[Ap@clusters[[i]]])] <- i
}
rm(Cl.master)

#Redefines Cl as a matrix
Cl <- matrix(Cl,nrows,ncols,byrow = FALSE)

### 4b. Post-processing of AP algorithm results ###

#This step reassigns individual xstal pixels that are on the boundary of clusters, but currently belong to background superpixels
bg.add <- function(x) {
  i <- id.temp[x,1] ; j <- id.temp[x,2]
  xr <- seq(i-1,i+1,1) ; xr <- xr[which(xr > 0 & xr <= nrow(mat.master))]
  xc <- seq(j-1,j+1,1) ; xc <- xc[which(xc > 0 & xc <= ncol(mat.master))]
  return(round(median(Cl[xr,xc],na.rm=TRUE),0))
}
id.temp <- which(is.na(Cl)==TRUE & is.na(mat.master[,,1])==FALSE,arr.ind = TRUE)
if(nrow(id.temp)>0) {Cl[id.temp] <- unlist(lapply(1:nrow(id.temp),bg.add))}
rm(id.temp)

#This makes all original background points NA
Cl[which(is.na(mat.master[,,1])==TRUE,arr.ind = TRUE)] <- NA

#This step converts Cl from matrix to raster to set all exterior pixels to NA
Cl <- raster(Cl)
extent(Cl) <- extent(inv.xstal)
Cl.temp <- flip(Cl,"y")
Cl.temp[exterior] <- NA
Cl <- flip(Cl.temp,"y")
rm(Cl.temp)

rm(Cen_vals,centroid,diss,cen.dif,cen1,cen2,i,val.temp,xstal.percent)
