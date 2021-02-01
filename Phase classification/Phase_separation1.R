
### 3a. Set working directory and load data ###
setwd(paste0("Data/",Sample_name))

#Extract elements in analyses and load data
Element_raw_data <- list.files(pattern = ".txt") 
# ps.id <- grep("Phases",Element_raw_data) # this is a check to see if the 
# if(length(ps.id)>0) {Element_raw_data <- Element_raw_data[-ps.id]}
invisible(lapply(Element_raw_data,function(x) {assign(unlist(strsplit(x,split = ".txt")),read.table(x,header=FALSE),envir = .GlobalEnv)}))
Elements <- unlist(strsplit(Element_raw_data,split = ".txt"))
Elements <- c(Elements[which(Elements=="Al")], Elements[which(Elements!="Al")])
invisible(lapply(Elements, function(x) {assign(x,value = matrix(unlist(get(x)),nrow = nrow(get(x)),byrow = FALSE),envir = .GlobalEnv)}))

rm(Element_raw_data)

# Defines dimensions of thin section ###
nrows <- nrow(get(Elements[1]))
ncols <- ncol(get(Elements[1]))

### 3b. Performs the finite mixture model for the first pair of elements ###

X1 <- (1*get(Elements1[1]))+(2*get(Elements2[1]))
X2 <- (2*get(Elements1[1]))+(1*get(Elements2[1]))
dfall <- cbind(as.vector(X1),as.vector(X2))

phases <- flexmix(dfall~1, k = nK, model = FLXmclust(diagonal = FALSE))

Results <-  eval(parse(text = "phases@cluster"))
if(class(Results)=="data.frame") {Results <- Results[,1]}

Con<- ifelse(phases@converged==TRUE,"Convergence","No convergence")
print(paste("1 of",length(Elements1),"-",Con))
Con.all <- Con

rm(X1,X2,phases,Con)

### 3c. Performs the finite mixture model for all other pairs of elements  ###

for (i in 2:length(Elements1)){
  
  mat.temp1 <- get(Elements1[i])
  mat.temp2 <- get(Elements2[i])
  X1 <- (1*mat.temp1)+(2*mat.temp2)
  X2 <- (2*mat.temp1)+(1*mat.temp2)
  dfall <- cbind(as.vector(X1),as.vector(X2))
  phases <- flexmix(dfall~1, k = nK, model = FLXmclust(diagonal = FALSE))
  Results.temp <-  eval(parse(text = "phases@cluster"))
  Results <- cbind(Results,Results.temp)
  
  Con <- ifelse(phases@converged==TRUE,"Convergence","No convergence")
  print(paste(i,"of",length(Elements1),"-",Con))
  Con.all <- c(Con.all,Con)
  
  rm(mat.temp1,mat.temp2,X1,X2,dfall,phases,Results.temp,Con)
}
rm(i)

Results <- data.frame(Results)
colnames(Results) <- paste0(Elements1,"_",Elements2,"_phases")

#Maximum number of clusters identified for each pairwise combination of elements
max.Nk <- apply(Results,2,max)

### 3d. Plot fmm results ###
#Sets working directory in Phase classification results
setwd("../../Phase classification/")
if(dir.exists(Sample_name)==FALSE) {dir.create(Sample_name)}
setwd(Sample_name)

plot.rows <- round(length(max.Nk)^0.5,0)
plot.cols <- ceiling(length(max.Nk)^0.5)
png(filename = paste0(Sample_name,"_fmm.png"),height = plot.rows*1000,width = plot.cols*1000)
par(mfrow=c(plot.rows,plot.cols))
lapply(seq(1,length(max.Nk),1), function(x) {plot(raster(matrix(Results[,x],ncol=nrows,byrow=FALSE)),col=rainbow(max.Nk[x]),main=paste(x,colnames(Results)[x]))})
dev.off()



