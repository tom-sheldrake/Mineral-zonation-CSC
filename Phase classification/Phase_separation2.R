
### 4a. Reduce to element pairs that have been chosen ###
Results.final <- Results[,final.choice]

#Calculates all distinct combinations of clusters
dis <- distinct(Results.final)
dis$N <- seq(1,nrow(dis),1)

#Joins dis$N onto Results dataframe, and convert to matrix
id <- join(Results.final,dis,by=intersect(names(Results.final), names(dis)))
col.id <- matrix(id$N,ncol = ncols,byrow = FALSE)

rm(id)

### 4b. Performs central pixel assumption to remove mixed pixels
grad <- col.id
for (i in 1:nrows){
  for (j in 1:ncols){
    xr <- seq(i-1,i+1,1) ; xr <- xr[which(xr > 0 & xr <= nrow(col.id))]
    xc <- seq(j-1,j+1,1) ; xc <- xc[which(xc > 0 & xc <= ncol(col.id))]
    mat.temp <- col.id[xr,xc]
    grad[i,j] <- var(as.vector(mat.temp))==0
  }
} ; rm(i,j,xr,xc,mat.temp)

components <- unique(col.id[which(grad==TRUE,arr.ind = TRUE)])
id <- match(unique(as.vector(col.id)),components)
Phases <- id[col.id]
# Phases[is.na(Phases)] <- 0

rm(components,id,grad,dis,col.id,Results.final)

### 3d. Plot final classification results ###

# Sets working directory in Phase classification results
setwd("Phase classification/")
setwd(Sample_name)

par(mfrow=c(1,1))
png(filename = paste(Sample_name,"_phases.png"),width = 1000,height = 1000, units = "px")
plot(raster(matrix(Phases,nrow = nrows,byrow = FALSE)), col=rainbow(length(unique(Phases))), main=paste(Sample_name))
dev.off()



