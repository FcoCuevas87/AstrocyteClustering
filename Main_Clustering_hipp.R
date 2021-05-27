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
source("Hierarchical_plot_hipp.R")

#setwd("./Images/Photos/")
files.vec <- c("GFAP_3.tif")
#setwd("./Images/Masks/")
mask_vec <- c("GFAP_3")
name_img <- strsplit(files.vec, ".tif")[[1]]

color.val =  c("red", "blue", "green3", "orange", "purple", "cyan", "pink")
setwd("./Results/Hippocampus")
if(!dir.exists(name_img)) dir.create(name_img)
setwd("../../")
Hierarchical_plot(Imagen1 = files.vec, mask_folder = mask_vec, 
                  criteria = "1l", color.val = color.val)
