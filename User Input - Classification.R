########## 0. PRELIMINARY CHECKS AND INFORMATION ##########

# i: Ensure the this file is in the CSC folder
# ii: Ensure raw data is in correct folder (e.g., CSC/Data/SK394C)
# iii: Ensure the raw data is in a text file, in matrix form, with no row or column headers

### make sure the following packages are installed ###
# install.packages("flexmix",dependencies = TRUE)
# install.packages("dplyr",dependencies = TRUE)
# install.packages("plyr",dependencies = TRUE)
# install.packages("raster",dependencies = TRUE)

### Instructions to run this code ###
#Once the user input in section 1 has been finalised, highlight sections 1 & 2, then press run.
#The output will be a file in "~/Phase classification/{Sample_name}/" containing a map of clusters for each element combination (e.g. "SK394C_fmm.png).
#Whether the algorithm convergences for each combination of elements will be printed in the R console.
#Use the output in the R console and the.png figure to make a final choice of which element combinations to use, and update section 3.
#Highlight section 3 and 4, then press run.

### Summary of output ###
#The output of the classifcation is saved as an R raster-stack in two files (.grd & .gri).
#Additionally, the final phase map will be saved as a .png file and a summary .txt file will contain the results and choices for the finite mixture model
#IMPORTANT: Each time the code is re-run, you will overwrite the previous results in "~/Crystal segmentation/{Sample name}/"


########## 1. USER INPUT REQUIRED - DEFINED THIN SECTION AND ELEMENTS TO USE ########## 

### Sample name as it appears in the Data folder ###
Sample_name <- "SK394C"

### Pairs of elements to perform the classification with ###
Elements1 <- c("Al","Al","Al","Al","Al","Si","Si","Si","Si") #Tetrahedral elements
Elements2 <- c("Ca","Fe","K","Mg", "Si","Ca","Fe","K", "Mg") #Mobile cations

### Maximum number of components identified by finite mixture model ###
nK <- 7

### Set home dir (leave this as default if using RStudio) ###
home.dir <- paste(dirname(rstudioapi::getActiveDocumentContext()$path))

########## 2. RUN THE FINITE MIXTURE MODEL ########## 

#Loads packages  
require(flexmix) 
require(dplyr)  
require(plyr)  
require(raster) 

#Sets working directory
setwd(home.dir)
source("Phase classification/Phase_separation1.R")

########## 3. IDENTIFY THE FINAL COMBINATIONS OF ELEMENTS TO MAKE THE PHASE MAP ########## 

### User input - Use combinations of elements to use ### 
final.choice <- c(1,2,3,4,6,7,8,9)


########## 4. PRODUCES PHASE MAP AND SAVES RESULTS ########## 

#Sets working directory & runs model
setwd(home.dir)
source("Phase classification/Phase_separation2.R")

invisible(lapply(seq(1,length(Elements),1), function(x) {assign(Elements[x],value = as.vector(get(Elements[x])),envir = .GlobalEnv)}))

data <- mapply(function(x) {get(Elements[x])},seq(1,length(Elements),1))
colnames(data) <- Elements

Output <- cbind(Phases,data,Results)

row.coord <- rep(0.5:(nrows-0.5),each=(ncols))
col.coord <- rep((ncols-0.5):0.5,nrows)
spg <- data.frame(row.coord,col.coord,Output)
coordinates(spg) <- ~ row.coord + col.coord
gridded(spg) <- TRUE
assign(x = Sample_name,value = stack(spg),envir = .GlobalEnv)
writeRaster(x=get(Sample_name), filename = paste0(Sample_name,".grd"),overwrite=TRUE,type="")

Included <- seq(1,length(Elements1),1) %in% final.choice

Summary <- rbind(Elements1,Elements2,Con.all,Included,max.Nk)
colnames(Summary) <- seq(1,length(Elements1),1)
write.table(Summary,"Summary.txt")





