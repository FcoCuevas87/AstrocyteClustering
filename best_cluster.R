library("cluster")
best_cluster <- function(dist.mat){
  agg.coef.com <- agnes(as.dist(dist.mat), method = "complete")
  agg.coef.sin <- agnes(as.dist(dist.mat), method = "single")
  agg.coef.avg <- agnes(as.dist(dist.mat), method = "average")
  agg.coef.war <- agnes(as.dist(dist.mat), method = "ward")
  
  
  agg.coef <- c("complete" = agg.coef.com$ac,
                "single" = agg.coef.sin$ac,
                "average" = agg.coef.avg$ac,
                "ward" = agg.coef.war$ac)
  return(agg.coef)
}
