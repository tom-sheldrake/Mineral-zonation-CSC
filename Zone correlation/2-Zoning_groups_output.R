####### CREATES PLOTS AND SAVES DATA FOR SELECTED NUMBER OF ZONING GROUPS ####### 

setwd(results.dir)

### 1a. Boxplot split into zoning groups

#Splits hierarchical cluster into groups, and assigns zone to respective group
branches <- as.data.frame(cutree(ht$tree_row, k = br))
branches$ID <- names(correl_var)
colnames(branches) <- c("branch", "ID")
rownames(branches) <- NULL

#Orders zones according to distance matrix
branches.order <- order(branches$branch)

#Sets plotting limits of shading boxes and zoning group numbers
box.ids <- c(-0.5,which(diff(branches$branch[branches.order])>0),nrow(branches)+1)
box.mids <- box.ids[-length(box.ids)] + diff(box.ids)/2
box.labels <- order(unique(branches$branch[ht$tree_row$order]))

#Function for boxplots
plot.boxes <- function(ymin,ymax,ylab1) {
  
  par(mai = c(1,1,1,1))
  plot(1, type = "n", ylim = c(ymin,ymax), xlim = c(0, nrow(branches)+0.5),
       xlab = "", xaxt = "n",ylab=ylab1,
       font = 2, cex.lab = 2, cex.axis = 2, yaxs = "i", xaxs = "i")
  
  axis(side = 1, lwd = 3, labels = FALSE, tick = F)
  axis(side = 2, lwd = 3, labels = FALSE, lwd.ticks = 3)
  
  lapply(1:max(branches$branch), function(x) {
    
    rect(xleft = box.ids[x]+0.5, xright = box.ids[x+1]+0.5,
                                  ytop = ymax, ybottom = ymin, col = rainbow(max(branches$branch))[x])})
  
  text(x= box.mids, y= ymax-((ymax-ymin)*0.9),labels=box.labels)
  
  #Plot the boxplots
  boxplot(lapply(strsplit(branches$ID[branches.order], split = " "), function (x) {eval(parse(text = paste0("correl_var$",x)))}),
          range =1.5, outline = F, add = T, xaxt = "n", yaxt = "n", whisklwd = 1.5,col="black",border=c("black"),medcol="white")
  
}

#Plot box plots
png(filename = paste0("boxes",br,".png"),height = 500,width = nrow(branches)*10)
plot.boxes(ymin = ymin1,ymax = ymax1,ylab1=correl_func_name)
dev.off()

### 1b. Heatmap

#Function to allow heatmap to be saved as .png
save_pheatmap_png <- function(x, filename ,width = 1000,height = 1000,units = "px") {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#Plot and save heatmap as .png
save_pheatmap_png(ht, paste0("Heatmap_",br,".png"))

### 2c. Readme file with settings for correlation

fileConn<-file("README.txt")
writeLines(c(
  paste0("Samples = ",list(Samples)),
  paste0("Elements = ",list(Correlating_elements)),
  paste0("Correl eq = ",list(strRep(correl_eq,"[[n]][[i]]","")))
),fileConn)
close(fileConn)

### 2d. Creates and saves raster files of results

for(i in 1:length(Samples)){
  lab.temp <- eval(parse(text = paste0(Samples[i],"_textures$Labels")))
  lab.temp[] <- branches[match(lab.temp[], branches$ID), "branch"]
  names(lab.temp) <- "Groups"
  assign("temp",stack(get(paste0(Samples[i],"_textures")),lab.temp,get(paste0(Samples[i],"_correl"))))
  id.rm <- unique(unlist(lapply(c("_", "Labels"),function(x) {grep(names(temp),pattern=x)})))
  temp <- temp[[-id.rm]]
  assign(paste0(Samples[i],"_zoning"),get("temp"))
  rm(lab.temp,id.rm,temp)
}
invisible(lapply(Samples, function(x) {writeRaster(filename = paste0(x,"_zoning",br),x=get(paste0(x,"_zoning")),overwrite=TRUE)} ))

### 2e. Plots zoning groups for each sample

plot.rows <- round(length(Samples)^0.5,0)
plot.cols <- ceiling(length(Samples)^0.5)
png(file=paste0(br,"zones.png"),width=1000*plot.cols,height=1000*plot.rows)
  par(mfrow=c(plot.rows,plot.cols))
  lapply(seq_along(Samples), function(x) {
  ras.temp <- eval(parse(text=paste0(Samples[x],"_zoning","$Groups")))
  id.temp <- unique(ras.temp$Groups[])[which(is.na(unique(ras.temp$Groups[])) == FALSE)]
  col.id <- rainbow(n = br)[min(id.temp):max(id.temp)]
  plot(ras.temp$Groups,col=col.id,main=Samples[x],legend=F,axes=F,box=T)})
  dev.off()
  
rm(plot.rows,plot.cols,fileConn,i)
save(list=ls()[-((grep("Raster",unlist(lapply(ls(),function(x) {class(get(x))}))))-1)],file = paste0(br,"groups.Rdata"))
