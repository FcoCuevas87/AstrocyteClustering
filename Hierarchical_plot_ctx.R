library('cluster')
Hierarchical_plot <- function(Imagen1, mask_folder,  criteria, color.val =  c("red", "blue", "green3", "orange", "purple", "cyan", "pink")){
  name_img <- strsplit(Imagen1, ".tif")[[1]]
  
  Imagen.raw1 <- as.cimg( t(readTIFF(paste("./Images/Photos/",Imagen1,sep = "")) ))
  Imagen.raw1 <- renorm(Imagen.raw1, min = 0, max = 1)
  
  Imagen.tmp <- Imagen.raw1*0
  
  ecdf.list <- indexes <- list()
  
  k <- 0
  names.vec <- NULL
  mask.list <- list()
  local.mask.name.vec <- NULL
  ini.folder <- getwd()
  for(j in 1:4){
    folder.name <- paste( "./Images/Masks/",
                          mask_folder, "/L", j, "/", sep="")
    setwd(folder.name)
    Mask.names <- dir()
    indexes[[j]] <- 1:length(Mask.names)
    print(length(Mask.names))
    X <- NULL
    
    for(i in 1:length(Mask.names)){
      k <- k + 1

      local.mask.name.vec[k] <- local.mask.name <- strsplit(Mask.names[i], ".tif")[[1]]
      
      Mask.tmp <- readTIFF(Mask.names[i])[,,1]
      Mask.tmp[Mask.tmp >= 0.5] <- 1
      Mask.tmp[Mask.tmp < 0.5] <- 0
      
      Mask.raw1 <- as.cimg( t(Mask.tmp))
      
      Circulos1 <- Imagen.raw1*0
      Circulos1[,,1,1] <- Imagen.raw1[,,1,1]* (1 - Mask.raw1[,,1,1])
      
      Mask.tmp <- matrix(1-Mask.raw1[,,1,1], nrow = nrow(Mask.tmp), ncol = ncol(Mask.tmp) )
      x.mean <- mean( which(colSums(Mask.tmp)  != 0), trim = 0.5 )
      y.mean <- mean( which(rowSums(Mask.tmp)  != 0), trim = 0.5 )
      
      mask.list[[k]] <- c(y.mean,x.mean, 25)
      
      ecdf.list[[k]] <- ecdf( Circulos1[Circulos1[,,1,1] != 0] )
      names.vec[k] <- paste("R",j,"c",local.mask.name,sep = "")
    }
    setwd(ini.folder)
  }
 
  setwd(paste("./Results/Cortex/",name_img,sep=""))
  
  ks.dist <- matrix(0, ncol = length(ecdf.list), nrow = length(ecdf.list))
  colnames(ks.dist) <- names.vec
  rownames(ks.dist) <- names.vec
  
  for(i in 1:length(ecdf.list)){
    for(j in 1:length(ecdf.list)){
      ecdfx <- ecdf.list[[i]](seq(0,1,l=1000))
      ecdfy <- ecdf.list[[j]](seq(0,1,l=1000))
      ks.dist[i,j] <- max( abs(ecdfx - ecdfy) )
    }
  }
  
  indexes[[2]] <- indexes[[2]] + max( indexes[[1]] )
  indexes[[3]] <- indexes[[3]] + max( indexes[[2]] )
  indexes[[4]] <- indexes[[4]] + max( indexes[[3]] )
  
  ks.dist1 <- ks.dist[indexes[[1]], indexes[[1]] ]
  ks.dist2 <- ks.dist[indexes[[2]], indexes[[2]] ]
  ks.dist3 <- ks.dist[indexes[[3]], indexes[[3]] ]
  ks.dist4 <- ks.dist[indexes[[4]], indexes[[4]] ]
  
  
  cluster.ind <- c("complete", "single", "average", "ward.D2")
  
  best1 <- which.max(best_cluster(ks.dist1))
  best2 <- which.max(best_cluster(ks.dist2))
  best3 <- which.max(best_cluster(ks.dist3))
  best4 <- which.max(best_cluster(ks.dist4))
  bestT <- which.max(best_cluster(ks.dist))
  
  dendoT <- hclust(as.dist(ks.dist), method = cluster.ind[bestT])
  dendo1 <- hclust(as.dist(ks.dist1), method = cluster.ind[best1])
  dendo2 <- hclust(as.dist(ks.dist2), method = cluster.ind[best2])
  dendo3 <- hclust(as.dist(ks.dist3), method = cluster.ind[best3])
  dendo4 <- hclust(as.dist(ks.dist4), method = cluster.ind[best4])

  dendoT.old <- hclust(as.dist(ks.dist), method = "complete")
  dendo1.old <- hclust(as.dist(ks.dist1), method = "complete")
  dendo2.old <- hclust(as.dist(ks.dist2), method = "complete")
  dendo3.old <- hclust(as.dist(ks.dist3), method = "complete")
  dendo4.old <- hclust(as.dist(ks.dist4), method = "complete")
  
  Kobs.new <- rep(0,5)
  Kobs.old <- rep(0,5)
  # Optimal
  plot.obj <- plot_dendro(dendo1, ks.dist1, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.new[1] <- plot.obj$K
  ggsave( paste(name_img,"_L1.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo2, ks.dist2, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.new[2] <- plot.obj$K
  ggsave( paste(name_img,"_L2.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo3, ks.dist3, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.new[3] <- plot.obj$K
  ggsave( paste(name_img,"_L3.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo4, ks.dist4, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.new[4] <- plot.obj$K
  ggsave( paste(name_img,"_L4.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendoT, ks.dist , k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.new[5] <- plot.obj$K
  ggsave( paste(name_img,"_all.pdf", sep = "") )
  
  # Complete
  plot.obj <- plot_dendro(dendo1.old, ks.dist1, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[1] <- plot.obj$K
  ggsave( paste(name_img,"_L1_old.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo2.old, ks.dist2, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[2] <- plot.obj$K
  ggsave( paste(name_img,"_L2_old.pdf", sep = "") )
  
  
  plot.obj <- plot_dendro(dendo3.old, ks.dist3, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[3] <- plot.obj$K
  ggsave( paste(name_img,"_L3_old.pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo4.old, ks.dist4, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[4] <- plot.obj$K
  ggsave( paste(name_img,"_L4_old.pdf", sep = "") )
  
  
  plot.obj <- plot_dendro(dendoT.old, ks.dist, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[5] <- plot.obj$K
  ggsave( paste(name_img,"_all_old.pdf", sep = "") )
  
  pdf("tanglegram11.pdf")
  dend_11 <- dendlist(as.dendrogram(dendo1.old), as.dendrogram(dendo1))
  dend_11 %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()
  
  pdf("tanglegram22.pdf")
  dend_22 <- dendlist(as.dendrogram(dendo2.old), as.dendrogram(dendo2))
  dend_22 %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()

  pdf("tanglegram33.pdf")  
  dend_33 <- dendlist(as.dendrogram(dendo3.old), as.dendrogram(dendo3))
  dend_33 %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()
  
  pdf("tanglegram44.pdf")  
  dend_44 <- dendlist(as.dendrogram(dendo4.old), as.dendrogram(dendo4))
  dend_44 %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()
  
  pdf("tanglegramTT.pdf")  
  dend_TT <- dendlist(as.dendrogram(dendoT.old), as.dendrogram(dendoT))
  dend_TT %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()
  
  ent.ind11 <- dend_11 %>% untangle(method = "step1side") %>% entanglement
  ent.ind22 <- dend_22 %>% untangle(method = "step1side") %>% entanglement
  ent.ind33 <- dend_33 %>% untangle(method = "step1side") %>% entanglement
  ent.ind44 <- dend_44 %>% untangle(method = "step1side") %>% entanglement
  ent.indTT <- dend_TT %>% untangle(method = "step1side") %>% entanglement
  
  sink("Best_cluster.txt")
  print("Agglomerative coefficient")
  print(best_cluster(ks.dist1))
  print( paste( "Entaglement index 1: ", ent.ind11) )
  
  print("Agglomerative coefficient")
  print(best_cluster(ks.dist2))
  print( paste( "Entaglement index 2: ", ent.ind22) )
  
  print("Agglomerative coefficient")
  print(best_cluster(ks.dist3))
  print( paste( "Entaglement index 3: ", ent.ind33) )
  
  print("Agglomerative coefficient")
  print(best_cluster(ks.dist4))
  print( paste( "Entaglement index 4: ", ent.ind44) )
  
  print("Agglomerative coefficient")
  print(best_cluster(ks.dist))
  print( paste( "Entaglement index all: ", ent.indTT) )
  
  #2d representation
  library(MASS)
  
  dend.new <- as.dendrogram(dendoT)
  dend.old <- as.dendrogram(dendoT.old)
  cluster.new <- dendextend::cutree(dend.new, k = Kobs.new[5])
  cluster.old <- dendextend::cutree(dend.old, k = Kobs.old[5])
  tree_order.new <- stats::order.dendrogram(dend.new)
  tree_order.old <- stats::order.dendrogram(dend.old)
  
  clustab.new <- table(cluster.new)[unique(cluster.new[tree_order.new])]
  clustab.old <- table(cluster.old)[unique(cluster.old[tree_order.old])]
  col.seq.new <- rep(1:Kobs.new[5],clustab.new)
  col.seq.old <- rep(1:Kobs.old[5],clustab.old)
  
  mds <- isoMDS(ks.dist)
  #remove {type = "n"} to see dots
  group.new <- cutree(dendoT,Kobs.new[5])
  group.old <- cutree(dendoT.old,Kobs.old[5])
  
  pdf(paste("2d_rep_new_",name_img,".pdf",sep = ""))
  plot(mds$points[tree_order.new,], pch = 20, cex = 3, col = adjustcolor( color.val[col.seq.new], alpha = 0.3), xlab = "X", ylab = "Y") 
  text(mds$points[tree_order.new,], labels = names(cluster.new)[tree_order.new], col = color.val[col.seq.new])
  dev.off()
  
  pdf(paste("2d_rep_old_",name_img,".pdf",sep = ""))
  plot(mds$points[tree_order.old,], pch = 20, cex = 3, col = adjustcolor( color.val[col.seq.old], alpha = 0.3), xlab = "X", ylab = "Y") 
  text(mds$points[tree_order.old,], labels = names(cluster.old)[tree_order.old], col = color.val[col.seq.old])
  dev.off()
  
  tiff(paste(name_img,"_colors.tiff",sep = ""), width = 1024, height = 1024)
  plot(Imagen.raw1, axes = FALSE)
  
  k <- 0
  for(i in tree_order.new){
    k <- k + 1
    y.o   <- mask.list[[i]][1]
    x.o   <- mask.list[[i]][2]
    rad.o <- mask.list[[i]][3]
    circles(y.o, x.o, rad.o  , fg = color.val[col.seq.new[k]])
    circles(y.o, x.o, rad.o+1, fg = color.val[col.seq.new[k]])
    circles(y.o, x.o, rad.o+2, fg = color.val[col.seq.new[k]])
    #text(y.o, x.o,i)
  }
  dev.off()
  
  print("Old number of clusters")
  print(Kobs.old)
  print("New number of clusters")
  print(Kobs.new)
  sink()
  setwd(ini.folder)
}
