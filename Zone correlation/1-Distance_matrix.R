####### CALCULATES CORRELATING PARAMETER AND DISTANCE MATRIX BETWEEN ALL ZONES ####### 

#Rearranges Samples in alphabetical order
Samples <- sort(Samples,decreasing = FALSE)

### 1a. Extract all zones in all crystals in all samples ###
all_zones <- lapply(ls(pattern = "_textures",envir = sys.frame()), function(x) {unique(eval(parse(text=paste0(x,"$Labels"))))})
names(all_zones) <- Samples
N <- unlist(lapply(all_zones,length))
zone_ID <- lapply(seq_along(Samples), function(n,i) {lapply(eval(parse(text=paste0("all_zones$",n[[i]]))), function(x,n2) {which(eval(parse(text=paste0(n2,"_textures$Labels@data@values")))==x)},n2=n[[i]])},n=Samples)

### 1b. Extract correlating elements ###

#Function to extract elements used for calibration and perform calibration (if required by user)
if(Calibration==TRUE) {
      correl_vars <- function(E) {lapply(seq_along(Samples), function(n,i) {lapply(zone_ID[[i]], function(x,n2,n3) {n3[1] + n3[2]*(eval(parse(text=paste0(n2,"_textures","$",E)))[x])},n2=n[[i]],n3=eval(parse(text=paste0(n[[i]],"_calibration","$",E))))},n=Samples)}
} else {
      correl_vars <- function(E) {lapply(seq_along(Samples), function(n,i) {lapply(zone_ID[[i]], function(x,n2) {eval(parse(text=paste0(n2,"_textures","$",E)))[x]},n2=n[[i]])},n=Samples)}
}

#Extract correlating elements
for (i in 1:length(Correlating_elements)){
  E <- Correlating_elements[i]
  assign(E,correl_vars(E))
}

### 1c. Calculate correlating variable ###

#Define function to calculate correlating variable
correl_eq <- strRep(correl_eq,".","[[n]][[i]]")
correl_func <- function(i,n) {eval(parse(text=correl_eq))}

#Calculate correlating variable for all samples in raster format
for (j in 1:length(Samples)){
  var.temp2 <- unlist(zone_ID[[j]])
  mat.temp2 <- eval(parse(text = paste0(Samples[j],"_textures$Zones")))
  mat.temp2[var.temp2] <- unlist(lapply(seq(1,length(zone_ID[[j]]),1),correl_func,j))
  names(mat.temp2) <- correl_func_name
  assign(paste0(Samples[j],"_correl"),mat.temp2)
  rm(mat.temp2,var.temp2)
}

#Unlist names of all zones to input into correlating variable into single list
all_zones <- unlist(all_zones)

#Create single list of all zones (correlating variable)
correl_var <- pblapply(seq(1,length(Samples),1),function(x) {lapply(seq(1,N[x],1), correl_func,n=x)})
correl_var <- unlist(correl_var,recursive = FALSE)
names(correl_var) <- all_zones

### 1d. Comparison of all zones to create distance matrix ###

#Set limits for plotting and k-s test
ymax1 <- max(pretty(unlist(correl_var),h=5))
ymin1 <- min(pretty(unlist(correl_var),h=5))
x.seq <- seq(ymax1,ymin1,by = -0.005)
options(warn=-1)

#create matrix of all combinations
combos <- expand.grid(1:sum(N),1:sum(N))

#function to calculate mean distance in cumulative probability distribution
chem.dist <- function(x) {
  zone.combo <- combos[x,]
  a <- correl_var[[as.numeric(zone.combo[1])]]
  b <- correl_var[[as.numeric(zone.combo[2])]]

  a <- a[which(a!=0)]
  b <- b[which(b!=0)] 
  
  a.cut <- cut(a,x.seq)
  a.cut <- table(a.cut) / length(a)
  a.cut <- 1- cumsum(a.cut)
  
  b.cut <- cut(b,x.seq)
  b.cut <- table(b.cut) / length(b)
  b.cut <- 1- cumsum(b.cut)
  
  stat.temp <- mean(((a.cut-b.cut)^2)^0.5)
  
  return(stat.temp)
}

#Reads in correlation matrix, or calculates new matrix
if(is.na(cor.input)==FALSE) {cor <- read.table(cor.input,header = TRUE)} else {

#calculates distance and converts into a matrix
cor  <- unlist(lapply(1:nrow(combos),chem.dist))
cor  <- matrix(cor,ncol = sum(N),byrow = TRUE)
options(warn=0)

#assigns row and column names to names of zones
rownames(cor) <- all_zones ; colnames(cor) <- all_zones
}
### 1e. Create ouput and calculate c-index ###

#Creates folder for specific date and time
results.dir <- date()
dir.create(results.dir)

#Saves distance matrix
setwd(results.dir)
write.table(cor,"cor.txt",row.names=FALSE)

#Function to calculate and plot C-Index
C.IND <- function(a,b){
diss_matrix<- dist(cor, method = "euclidean", diag=FALSE)
hc <- hclust(d = diss_matrix,method = "average")
NT <- (nrow(cor)*(nrow(cor)-1))/2
C_IND <- NULL
Smin <- NULL
Smax <- NULL
SW <- NULL
for(i in a:b){
  ind.tmp <- i-a+1
  memb <- cutree(hc, k = i)
  N_memb <- as.numeric(table(memb))
  NW <- (N_memb*(N_memb-1))/2
  NW <- sum(NW)
  S <- sort(as.numeric(diss_matrix),decreasing = FALSE)[1:NW] ; Smin[ind.tmp] <- sum(S)
  S <- sort(as.numeric(diss_matrix),decreasing = TRUE)[1:NW] ; Smax[ind.tmp] <- sum(S)
  within.clus <- NULL
  for(j in 1:i){
    id.temp <- which(memb==j)
    coords<- expand.grid(id.temp,id.temp)
    within.clus[j] <- sum(diss_matrix[as.matrix(coords)])
  }
  SW[ind.tmp] <- sum(within.clus)
  C_IND[ind.tmp] <- (SW[ind.tmp] - Smin[ind.tmp]) / (Smax[ind.tmp] - Smin[ind.tmp])
}

png(filename = paste0("C-Index.png"),width = 1000,height = 500,units = "px")
par(mfrow=c(1,1))
plot(seq(a,b,1),C_IND,type="b",lwd=2,cex=2,main="C-Index",ylab="C-Index",xlab="Zoning groups")
dev.off()

plot(seq(a,b,1),C_IND,type="b",lwd=2,cex=2,main="C-Index",ylab="C-Index",xlab="Zoning groups")

}

rm(combos,x.seq,i,j,E)
rm(list = Correlating_elements)


