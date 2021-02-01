########## 0. PRELIMINARY CHECKS AND INFORMATION ##########

# i: Ensure the this file is in the CSC folder
# ii: Ensure that the segmentation for all thin sections has been run an is contained in "~/Phase classification/{Sample name}/" as an Rdata file
# iii: Ensure that if Calbration == TRUE, a .txt file is located in "~/Data/{Sample name}" with the same name structure (e.g., "SK394C_calibration.txt")
# iv: This structure of this file should be the same as in the example, with the relevant elements in the columns

### make sure the following packages are installed ###
# install.packages("raster",dependencies = TRUE)
# install.packages("pheatmap",dependencies = TRUE)
# install.packages("pracma",dependencies = TRUE)
# install.packages("pbapply",dependencies = TRUE)


### Instructions to run this code ###
#Once the user input in section 1 has been finalised, highlight sections 1 & 2, then press run.
#Wait for the code to fiNish running, then plot the C-Index.
#Use the C-Index result to choose the total number of zoning groups then highlight and run section 4.
#Section 4 can be re-run to change the number of zoning groups, which can also be informed by the distance matrix that is plotted.
#Once a final number of zoning groups has been chosen, run and highlight section 5.


### Summary of output ###
#The output of the correlation is saved as an folder with the date and time in "~/Zone correlation/{Sample name}/"
#Each time this script is run from the beginning, a new folder with be created.
#If only sections 4 & 5 are re-run, but for a different number of zoning groups, all results and figures will be output in the same output folder.
#Each output folder contains, for each number of total zoning groups chosen:
#1. Plot of all zoning groups in all thin sections.
#2. Bar plot of each zone in each zoning group, using the same colour scale as the thn section map.
#3. Distance matrix split into zoning groups (sequnetially from 1 at top left)
#4. Plot of the C-Index
#5. R raster-stacks for each thin section (".gri" & ".grd")


########### 1. USER INPUT REQUIRED - DEFINES THIN SECTION AND PARAMETER TO CORRELATE ########### 

### Choose samples to correlate - it can also be only one sample ###
Samples <- c("SK394C","SK394A")

### Choose elements or measures of intensity (e.g.,BSE) that are used to calculate the correlating variable ###
Correlating_elements <- c("Ca","Na")

### Define the correlating variable using the above chosen elements or intensities ###
### If a single element or intensity it does not need to be an equation ###
### Ensure that a full stop follows after every element or intensity ###
correl_eq <- "Ca. / (Ca. + Na.)"

### Define the correlating variable name ###
correl_func_name <- "Anorthite"

### State whether the data will be calibrated. If yes then type TRUE, ortherwise FALSE. ###
### If the segmentation is already performed on calibrated data this variable should be set to FALSE. ###
Calibration <- TRUE

### Set home dir (leave this as default if using RStudio) ###
home.dir <- paste(dirname(rstudioapi::getActiveDocumentContext()$path))

### Set file path for previous correlation matrix (within "Zone correlation" folder), or set to NA to calculate new matrix ###
cor.input <- NA
# cor.input <- "Mon Jan 18 14:27:15 2021/cor.txt"

########## 2. LOAD IN THE THIN SECTIONS AND CALCULATE DISTANCE MATRIX FOR ALL ZONES ########## 

library(raster)
library(pheatmap)
library(pracma)
library(pbapply)
options(scipen = 999)

'%!in%' <- function(x,y)!('%in%'(x,y))

setwd(home.dir)

setwd("Crystal segmentation/")
for (i in 1:length(Samples)){
  load(paste0(Samples[i],"/",Samples[i],"_Zones.Rdata"))
  assign(paste0(Samples[i],"_crystals"),crystals)
  assign("ras.temp",stack(paste0(Samples[i],"/",Samples[i],"_segmented.grd")))
  
  #Converts numeric labels to text labels
  label.names <- labels.id[ras.temp$Labels[]]
  ras.temp$Labels[which(is.na(label.names)==FALSE)] <- label.names[which(is.na(label.names)==FALSE)]
  ras.temp$Labels[which(is.na(label.names)==TRUE)] <- NA
  
  assign(paste0(Samples[i],"_textures"),ras.temp)
  rm(ras.temp)
}

rm(list = ls()[-grep(ls(),pattern = paste0(paste0(Samples,"_|",collapse = ""),"correl_|","Samples|","Correlating_|","home.dir|","Calibration|","cor.input"))])

if(Calibration==TRUE) {
setwd("../Data/")
for (i in 1:length(Samples)){
assign(paste0(Samples[i],"_calibration"), read.table(paste0(Samples[i],"/",Samples[i],"_calibration.txt")))
}
}

setwd(home.dir)
setwd("Zone correlation/")
source("1-Distance_matrix.R")

########## 3. USER INPUT REQUIRED - CALCULATE AND PLOT C-INDEX ########## 

### Set a & b to the minimum (> 2) and maximum possible number of zoning groups, and plot C-Index ###
C.IND(a=2,b=25)


########## 4. USER INPUT REQUIRED - CHOOSE NUMBER OF ZONING GROUPS  ########## 

### Choose final number of zoning groups ###
br <- 6

#Plot a heatmap using pretty heatmap
ht <- pheatmap(1-cor, cutree_rows = br, cutree_cols = br, display_numbers = F,show_rownames  = F,show_colnames =  F,clustering_method = "average")


########## 5. SAVES OUTPUT FIGURES AND RESULTS  ########## 

setwd(home.dir)
setwd("Zone correlation/")
source("2-Zoning_groups_output.R")

