library('cluster')
Hierarchical_plot_hipocampo <- function(Imagen1, mask_folder,  criteria, color.val =  c("red", "blue", "green3", "orange", "purple", "cyan", "pink")){
  name_img <- strsplit(Imagen1, ".tif")[[1]]
  
  Imagen.raw1 <- as.cimg( t(readTIFF(paste("./Images/Photos/",Imagen1,sep = "")) ))
  Imagen.raw1 <- renorm(Imagen.raw1, min = 0, max = 1)
  
  ecdf.list <- indexes <- list()
  
  k <- 0

  inner.list <- c("Oriens", "Pyramidale", "Radiatum")
  label_num <- as.numeric(gsub("[^0-9.]", "",  mask_folder))
  names.vec <- NULL
  mask.list <- list()
  local.mask.name.vec <- NULL
  
  ini.folder <- getwd()
  for(j in 1:3){
    folder.name <- paste( "./Images/Masks/",
                          mask_folder, "/", inner.list[j], "/", sep="")
    setwd(folder.name)
    Mask.names <- dir()
    indexes[[j]] <- 1:length(Mask.names)
    print(length(Mask.names))
    X <- NULL
    for(i0 in 1:length(Mask.names)){
      k <- k + 1
      local.mask.name.vec[k] <- local.mask.name <- strsplit(Mask.names[i0], ".tif")[[1]]
      
      Mask.tmp <- readTIFF(Mask.names[i0])[,,1]
      Mask.tmp[Mask.tmp >= 0.5] <- 1
      Mask.tmp[Mask.tmp < 0.5] <- 0
      
      Mask.raw1 <- as.cimg( t(Mask.tmp))
      
      Circulos1 <- Imagen.raw1*0
      Circulos1[,,1,1] <- Imagen.raw1[,,1,1]* (1 - Mask.raw1[,,1,1])
      
      Mask.tmp <- matrix(1-Mask.raw1[,,1,1], 1024,1024 )
      x.mean <- mean( which(colSums(Mask.tmp)  != 0), trim = 0.25 )
      y.mean <- mean( which(rowSums(Mask.tmp)  != 0), trim = 0.25 )
      
      mask.list[[k]] <- c(y.mean,x.mean, 25)
      
      ecdf.list[[k]] <- ecdf( Circulos1[Circulos1[,,1,1] != 0] )
      names.vec[k] <- paste(substr(inner.list[j],1,3),"_c",local.mask.name,sep = "")
    }
    setwd(ini.folder)
  }
  
  setwd(paste("./Results/Hippocampus/",name_img,sep=""))
  
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
  
  ks.dist1 <- ks.dist[indexes[[1]], indexes[[1]] ]
  ks.dist2 <- ks.dist[indexes[[2]], indexes[[2]] ]
  ks.dist3 <- ks.dist[indexes[[3]], indexes[[3]] ]
  
  cluster.ind <- c("complete", "single", "average", "ward.D2")
  
  best1 <- which.max(best_cluster(ks.dist1))
  best2 <- which.max(best_cluster(ks.dist2))
  best3 <- which.max(best_cluster(ks.dist3))
  bestT <- which.max(best_cluster(ks.dist))
  
  dendoT <- hclust(as.dist(ks.dist ), method = "ward.D2")
  dendo1 <- hclust(as.dist(ks.dist1), method = "ward.D2")
  dendo2 <- hclust(as.dist(ks.dist2), method = "ward.D2")
  dendo3 <- hclust(as.dist(ks.dist3), method = "ward.D2")
  
  dendoT.old <- hclust(as.dist(ks.dist ), method = "complete")
  dendo1.old <- hclust(as.dist(ks.dist1), method = "complete")
  dendo2.old <- hclust(as.dist(ks.dist2), method = "complete")
  dendo3.old <- hclust(as.dist(ks.dist3), method = "complete")
  
  Kobs.new <- Kobs.old <- NULL
  
  #Optimal
  plot.d <- plot_dendro(dendo1, ks.dist1, k.max = 10, criteria, color.val)
  Kobs.new[1] <- plot.d$K
  plot.d$plot
  ggsave( paste(name_img,"_",inner.list[1],".pdf", sep = "") )
  
  plot.d <- plot_dendro(dendo2, ks.dist2, k.max = 10, criteria, color.val)
  Kobs.new[2] <- plot.d$K
  plot.d$plot
  ggsave( paste(name_img,"_",inner.list[2],".pdf", sep = "") )
  
  plot.d <- plot_dendro(dendo3, ks.dist3, k.max = 10, criteria, color.val)
  Kobs.new[3] <- plot.d$K
  plot.d$plot
  ggsave( paste(name_img,"_",inner.list[3],".pdf", sep = "") )
  
  plot.d <- plot_dendro(dendoT, ks.dist, k.max = 10, criteria, color.val)
  Kobs.new[4] <- plot.d$K
  plot.d$plot
  ggsave( paste(name_img,"_all.pdf", sep = "") )
  
  #Complete
  plot.obj <- plot_dendro(dendo1.old, ks.dist1, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[1] <- plot.obj$K
  ggsave( paste(name_img,"_old_",inner.list[1],".pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo2.old, ks.dist2, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[2] <- plot.obj$K
  ggsave( paste(name_img,"_old_",inner.list[2],".pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendo3.old, ks.dist3, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[3] <- plot.obj$K
  ggsave( paste(name_img,"_old_",inner.list[3],".pdf", sep = "") )
  
  plot.obj <- plot_dendro(dendoT.old, ks.dist, k.max = 10, criteria, color.val)
  plot.obj$plot
  Kobs.old[4] <- plot.obj$K 
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
  
  pdf("tanglegramTT.pdf")  
  dend_TT <- dendlist(as.dendrogram(dendoT.old), as.dendrogram(dendoT))
  dend_TT %>% untangle(method = "step1side") %>% 
    tanglegram(common_subtrees_color_branches = TRUE)
  dev.off()
  
  ent.ind11 <- dend_11 %>% untangle(method = "step1side") %>% entanglement
  ent.ind22 <- dend_22 %>% untangle(method = "step1side") %>% entanglement
  ent.ind33 <- dend_33 %>% untangle(method = "step1side") %>% entanglement
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
  print(best_cluster(ks.dist))
  print( paste( "Entaglement index all: ", ent.indTT) )
  sink()
  #2d representation
  library(MASS)
  
  dend.new <- as.dendrogram(dendoT)
  dend.old <- as.dendrogram(dendoT.old)
  cluster.new <- dendextend::cutree(dend.new, k = Kobs.new[4])
  cluster.old <- dendextend::cutree(dend.old, k = Kobs.old[4])
  tree_order.new <- stats::order.dendrogram(dend.new)
  tree_order.old <- stats::order.dendrogram(dend.old)
  
  clustab.new <- table(cluster.new)[unique(cluster.new[tree_order.new])]
  clustab.old <- table(cluster.old)[unique(cluster.old[tree_order.old])]
  col.seq.new <- rep(1:Kobs.new[4],clustab.new)
  col.seq.old <- rep(1:Kobs.old[4],clustab.old)
  
  mds <- isoMDS(ks.dist)

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
  }
  
  dev.off()
  setwd(ini.folder)  
}

