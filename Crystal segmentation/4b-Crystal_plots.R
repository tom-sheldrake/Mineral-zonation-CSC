####### PRODUCES FIGURE OF MASK, GRADIENT, SUPERPIXELS, & ZONING FOR EACH CRYSTAL (IF USER DEMANDS) ####### 

### 4b-a. Defining plotting dimensions, and opens .png plotting ###
setwd(Sample_name)

#Plot borders
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Creates polygon line for exterior of all superpixels
for(i in 1:length(superpixel.id)){
  assign(paste0("poly.sp",i),rasterToPolygons(sps,fun=function(x){x==superpixel.id[i]},n=4,dissolve = TRUE),envir = globalenv())
} ; rm(i)
assign("superpixel.crystals",bind(lapply(ls(pattern = "poly.sp"),get,envir = globalenv())),envir = globalenv())
rm(list=ls(pattern = "poly.sp"),envir=globalenv())

#Creates .png file
png(filename = paste0("Crystal",N,".png"),width = wd,height = ht,units = "px")
par(mfrow=c(2,3))

### 4b-b. Plot of crystal mask ###
plot(inv.xstal$Phases,zlim=c(Phase,Phase),legend=F,col=c("lightgrey","darkgrey"))
plot(Which(inv.xstal$Phases==-99),add=T)
lines(coords.temp[,"x"], coords.temp[,"y"], type = "l",lwd=3)

### 4b-c. Plot of gradient ###
plot(log(Gradient),main="ln(Gradient)",col=magma(100))
points(Ck.initial[,2]+xmin-0.5,nrows-Ck.initial[,1]+0.5+ymin,pch=21,bg=rgb(0.7,0.7,0.7,1),col="red",cex=1.5)
points(Ck.grid[,2]+xmin-0.5,nrows-Ck.grid[,1]+0.5+ymin,pch=21,bg=rgb(0.7,0.7,0.7,1),col="black",cex=1.5)
lines(coords.temp[,"x"], coords.temp[,"y"], type = "l",lwd=3)

### 4b-d. Plot of superpixels with centroids from each iteration of SLIC algorithm ###
plot(sps,main="Superpixels",col=viridis(length(superpixel.id)))
for (i in 1:(length(ls(pattern = "centroids.temp"))-1)){
  points.temp <- get(paste0("centroids.temp",i))
  points(points.temp,col=rgb(1,1,1,0.5),pch=21,bg=rgb(1,1,1,0.5),cex=1.5)
  rm(points.temp)
}
final.centroids <- get(paste0("centroids.temp",length(ls(pattern = "centroids.temp"))))
points(final.centroids,col="red",pch=21,bg=rgb(1,1,1,0.5),cex=1.5)

### 4b-e. Plot of superpixels with final centroids numbered ###
plot(sps,main="Superpixels",col=viridis(length(superpixel.id))) ; lines(superpixel.crystals,col="darkgrey")
points(final.centroids,col="red",pch=21,bg=rgb(1,1,1,1),cex=3)
text(final.centroids[,1],final.centroids[,2],seq(1,length(superpixel.id),1),cex=1)
lines(coords.temp[,"x"], coords.temp[,"y"], type = "l",lwd=3)

### 4b-f. Plot results of AP algorithm ###
plot(Cl,col=terrain.colors(n.clus+1)[1:n.clus],main="Zoning") ; lines(superpixel.crystals,col="darkgrey")
##This identifies the exemplar superpixels
for (i in 1:n.clus){
  c <- rasterToPolygons(sps,fun=function(x){x==superpixel.id[id.xstal[Ap@exemplars[i]]]},n=4,dissolve = TRUE)
  lines(c,col="red",lwd=1.5)
}
lines(coords.temp[,"x"], coords.temp[,"y"], type = "l",lwd=3)

dev.off()


rm(i,c,superpixel.crystals,wd,ht)




