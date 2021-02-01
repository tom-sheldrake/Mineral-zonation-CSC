####### POST PROCESSING STEP TO REMOVE SINGLE PIXELS AND DUPLICATE CENTROIDS #######

### 1a. Removes single pixels ###

#Defines function to find all pixels that do not belong to the same superpixel as at least one pixel to the [left,right] or [up,down].
single.func <- function(x,Mat) {
  i <- coords[x,1] ; j <- coords[x,2]
  xr <- seq(i-1,i+1,1) ; xr <- xr[which(xr > 0 & xr <= nrow(Mat))]
  xc <- seq(j-1,j+1,1) ; xc <- xc[which(xc > 0 & xc <= ncol(Mat))]
  diff <- Mat[xr,xc] - Mat[i,j]
  test <- length(which(diff==0)) > 1
  if(test==FALSE) {
    lc <- unique(as.vector(Mat[xr,xc])) ; lc <- lc[-which(lc==Mat[x])]
    return(lc[which.min((((Ck[lc,1]-i)^2) + ((Ck[lc,2]-j)^2))^0.5)])} else {return(l[i,j])}
  } 

#Reassigns single pixels to the closest centroid of one of the superpixels to the [left,right] or [up,down].
l <- unlist(lapply(1:nrow(coords), single.func,Mat=l))

#Redefines superpixels, in case a superpixel is removed
id.single <- sort(unique(l))
if((length(id.single) < nrow(Ck))==TRUE) {
  singles <- which(is.na(match(1:nrow(Ck),id.single))==TRUE)
  l <- match(1:nrow(Ck),id.single)[l]
  Ck <- Ck[-singles,]
  rm(singles)
}

#Converts l back to a matrix of superpixels
l <- matrix(as.numeric(l),nrow = nrows,byrow = FALSE)


### 1b. Removes duplicate centroids ###

new.l <- seq(1,nrow(Ck),1)
duplicates <- which(duplicated.matrix(Ck)==TRUE)
if (length(duplicates)>0){
  for (t in 1:length(duplicates)){
    d.replace <- which(Ck[,1]==Ck[duplicates[t],1] & Ck[,2]==Ck[duplicates[t],2])[1]
    l[which(l==duplicates[t],arr.ind = TRUE)] <- d.replace
  }
  Ck <- Ck[-duplicates,]
  new.l <- match(new.l,new.l[-duplicates])
  l <- matrix(new.l[l],nrow = nrows,ncol = ncols)
}

rm(duplicates,id.single,new.l)
