####### SLIC ALGORITHM ####### 

### 2a. SLIC algorithm ###

#Defines the initial superpixel matrix (which superpixel each pixel belongs to)
l <- matrix(-1,nrows,ncols)
#Defines the initial minimum distance between each pixel and a centroid
d <- matrix(Inf,nrows,ncols)

pb = txtProgressBar(min = 0, max = N.Iters, initial = 0) 

for (j in 1:N.Iters){
  for (i in 1:nrow(Ck)){
    center <- Ck[i,] ; r.c <- center[1] ; c.c <- center[2]
    search.radius <- if(is.na(search.radius.manual)==TRUE) {round(S*2)} else (search.radius.manual)
    r <- seq(r.c-search.radius,r.c+search.radius,1) ; r <- r[which(r>0)] ; r <- r[which(r<=nrows)]
    c <- seq(c.c-search.radius,c.c+search.radius,1) ; c <- c[which(c>0)] ; c <- c[which(c<=ncols)]
    sp <- expand.grid(r,c)
    r.p <- as.numeric(sp[,1])
    c.p <- as.numeric(sp[,2])
    d1 <- as.numeric(unlist(lapply(seq(1,nrow(sp),1), function(x) { sum(unlist(lapply(1:N.mat, function(N) {score.funct(arr = mat,r1 = r.c,c1 = c.c,r2 = r.p[x],c2 = c.p[x],N=N) })))})))
    d1 <- d1^0.5
    d2 <- as.numeric(unlist(lapply(seq(1,nrow(sp),1), function(x) { ((r.p[x]-r.c)^2 + (c.p[x]-c.c)^2)^0.5 })))
    D <-  (d1^2 + ((d2^2)/S)*(M^2))^0.5 
    test <- D < d[as.matrix(sp)] 
    test.true <- which(test==TRUE)
    l[as.matrix(sp[test.true,])] <- i
    d[as.matrix(sp[test.true,])] <- D[test.true]
    rm(D)
  }
  l <- matrix(match(1:nrow(Ck),sort(unique(as.vector(l))))[l],ncol=ncol(l),byrow = FALSE)  #This ensures empty indices are removed
  l.temp <- as.vector(l)
  matches <- splitAt(order(l.temp),which(diff(l.temp[order(l.temp)])>0)+1)
  new.cen <- lapply(seq(1,nrow(Ck),1),function(x) {round(apply(coords[matches[[x]],],2,mean),0)})
  new.cen1 <- lapply(seq(1,nrow(Ck),1),function(x) {which.min(apply((sweep(coords[matches[[x]],],2,STATS = new.cen[[x]],FUN = "-"))^2,1,sum))})
  Ck <- coords[as.numeric(names(unlist(new.cen1))),]
  med <- lapply(seq(1,nrow(Ck),1),function(x) {median(unlist(mat[as.matrix(coords[matches[[x]],])]),na.rm = TRUE)})
  med <- lapply(seq(1,nrow(Ck),1),function(x) {unlist(lapply(1:N.mat,function(y,x1) {median(mat[as.matrix(cbind2(coords[matches[[x1]],],rep(y,length(matches[[x1]]))))],na.rm=TRUE)},x1=x))})
  for (k in 1:nrow(Ck)){mat[Ck[k,1],Ck[k,2],] <- med[[k]]} ; rm(k)
  # Ck <- Ck[which(is.nan(Ck[,1])==FALSE),]
  Ck <- as.matrix(Ck)
  # print(paste0(j," of ",N.Iters))
  
  centroids.temp <- data.frame(Ck[,2]+xmin-0.5,ymax-Ck[,1]+0.5)
  assign(paste0("centroids.temp",j),centroids.temp,pos = .GlobalEnv)
  rm(centroids.temp)
  
  setTxtProgressBar(pb,j)
} ; rm(i,j)

rm(matches,med,new.cen,new.cen1,sp,pb,
   c,c.c,c.p,center,d1,d2,l.temp,r,r.c,r.p,test,test.true)

