localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

# dendo.obj <- dendoT
# dist.mat <- ks.dist

plot_dendro <- function(dendo.obj, dist.mat, k.max = 15, criteria = "1l", color.val){
  
  k.max.teo <- floor(0.5*nrow(dist.mat))
  k.max <- pmin(k.max.teo, k.max)
  
  avg.silhouette <- rep(0,length(2:k.max))
  for(k in 2:k.max){
    si4 <- silhouette(cutree(dendo.obj, k = k), dist.mat)
    si4 <- summary(si4)
    avg.silhouette[k-1] <- si4$avg.width
  }
  
  if(criteria == "max"){ 
    max.sil <- which.max(avg.silhouette)
  }
  
  if(criteria == "1l"){ 
    max.sil<- min(localMaxima(avg.silhouette))
  }
  
  Kobs <- (2:k.max)[max.sil]
  
  plot.obj <- fviz_dend( as.dendrogram(dendo.obj), k =  Kobs, # Cut in four groups
            cex = 0.5,                 # label size
            k_colors = color.val[1:Kobs],
            color_labels_by_k = TRUE,  # color labels by groups
            ggtheme = theme_gray(),     # Change theme
            main = paste(Kobs, "clusters"),
            ylab = "Distance",
            rect_fill = TRUE
  )
  
  return(list(K = Kobs, plot = plot.obj, avg.sil = avg.silhouette))
}
