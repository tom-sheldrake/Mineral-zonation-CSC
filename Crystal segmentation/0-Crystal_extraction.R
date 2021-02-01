####### EXTRACTS ALL CRYSTALS IN SAMPLE FOR SEGMENTATION ####### 

### Oa. Set working directory and load data ###
setwd(paste0("../Phase classification/",Sample_name))

#Assign the sample
assign("ras", stack(paste0(Sample_name,".grd")))


### Ob. Identify phase and extract crystals ###

#Find values in the raster that match the chosen phase P
Phase.segment <- Which(ras$Phases %in% Phase)
names(Phase.segment) <- "Crystals"
if (exists("crystals")==FALSE) {
  #raster to polygons the plag
  p1 <- rasterToPolygons(Phase.segment,n = 4,dissolve = TRUE,fun = function(x) {x==1},)
  #Unlist the polygons
  crystals <- unlist(p1@polygons[[1]])
  rm(p1)
}

#Identify if polygon is a hole
holes <- unlist(lapply(crystals@Polygons, function(x) {x@hole==TRUE}))

# Remove polygon if it is a hole
crystals@Polygons[which(holes==TRUE)] <- NULL

#Calculate area of each crystal polygon, and re-order crystals in decreasing area
areas <- unlist(lapply(crystals@Polygons, function(x) {eval((text=x@area))}))
areas.order <- order(areas,decreasing = TRUE)
crystals <- Polygons(crystals@Polygons[areas.order],ID=NA)
areas <- unlist(lapply(crystals@Polygons, function(x) {eval((text=x@area))}))
rm(areas.order)

#Get exterior of each crystal
Exterior <- lapply(crystals@Polygons, function(x) {eval((text=x@coords))})

#Identify crystals for segmentation based on minimum size
min.px.area <- (ifelse(((min.area^0.5)/pix.mm)<9, 81,(min.area^0.5)/pix.mm ))^2
xstal.to.segment <- which(areas > min.px.area)

#Labels of crystals to be segmentated
labpt <- matrix(,nrow=length(crystals@Polygons),ncol=2)
for (i in 1:length(crystals@Polygons)){
  labpt[i,] <- crystals@Polygons[[i]]@labpt  
}


### Oc. Create plot of all crystals, and identify those to be segmented ###
setwd("../../Crystal segmentation/")
if(dir.exists(Sample_name)==FALSE) {dir.create(Sample_name)}
setwd(Sample_name)

par(mfrow=c(1,1))
png(filename = paste0(Sample_name,".png"),width = 1000,height = 1000,units = "px")
plot(Phase.segment,col=grey.colors(2),main=Sample_name,legend=F)
lapply(1:length(xstal.to.segment),function(x) {lines(Exterior[[xstal.to.segment[x]]],lty=1,lwd=1.5,col="red")})
text(labpt[1:length(xstal.to.segment),],labels=1:length(xstal.to.segment),col="red",cex=1, font = 2)
dev.off()

rm(i)
