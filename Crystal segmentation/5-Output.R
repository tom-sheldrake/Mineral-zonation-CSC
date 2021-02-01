######## PREPARES DATA TO BE SAVED INTO OUTPUT RASTER IN RDATA FILE ######

### 5a. Creation of output raster ###

#Converts NA to zero so does not overwrite previous crystals
Cl[is.na(Cl)==TRUE] <- 0 

###this converts cluster from raster to matrix
# Cl.mat <- matrix(Cl$layer,nrow = nrows,ncol = ncols,byrow = FALSE)
e <- extent(inv.xstal)
textures[e] <- textures[e] + Cl$layer

#Input the labels into a layer
for (j in 1:n.clus){
  id.temp <- which(Cl$layer[]==j)
  label.ras[e][id.temp] <- 99
  label.ras[e][id.temp] <- paste0(Sample_name,"Zone",N,letters[j])
}

rm(id.temp,e,j)

