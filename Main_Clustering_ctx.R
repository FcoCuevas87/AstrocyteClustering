rm(list=ls())
library("SpatialPack")
library("geoR")
library("imager")
library("tiff")
library("factoextra")
library("dendextend")

setwd("C:/Users/franc/Desktop/R_Codes/Yorka/Archivos/Codes/Github/")
source("Plot_function.R")
source("best_cluster.R")
source("Hierarchical_plot_ctx.R")

#setwd("./Images/Photos/")
files.vec <- c("GLT1_1.tif")
#setwd("./Images/Masks/")
mask_vec <- c("GLT1_1")
name_img <- strsplit(files.vec, ".tif")[[1]]

color.val =  c("red", "blue", "green3", "orange", "purple", "cyan", "pink")
setwd("./Results/Cortex")
if(!dir.exists(name_img)) dir.create(name_img)
setwd("../../")
Hierarchical_plot(Imagen1 = files.vec, mask_folder = mask_vec, 
                            criteria = "1l", color.val = color.val)
