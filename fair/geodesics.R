library(rgl) 
library(igraph)

geodesicsFun <- function(name="", vertexNum)
{
  object <- readOBJ(name)
  vertex.xyz <- t(object[[1]][1:3,])
  edges <- unique(rbind(t(object[[2]][1:2,]),t(object[[2]][2:3,]),t(object[[2]][c(3,1),])))
  edgesLen <- sqrt(apply((vertex.xyz[edges[,1],] - vertex.xyz[edges[,2],])^2,1,sum))
  distance <- distances(graph.edgelist(edges), weights=edgesLen)
  normalizeDistance <- distance/max(distance)
  
  geodesicsColor <- ceiling(normalizeDistance[vertexNum,] * 15)+1
  
  output <- list(vertex.xyz, normalizeDistance, geodesicsColor)
  return (output)
}

vertexNum <- 30

#bunny geodesics
bunnyRet = geodesicsFun("bunnytest.obj", vertexNum)

open3d()
plot3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3], xlab = "bunny geodesics")
spheres3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3],radius=0.03,col=bunnyRet[[3]])

#sphere geodesics
sphereRet = geodesicsFun("spheretest.obj", vertexNum)

open3d()
plot3d(sphereRet[[1]][, 1], sphereRet[[1]][, 2], sphereRet[[1]][, 3], xlab = "sphere geodesics")
spheres3d(sphereRet[[1]][, 1], sphereRet[[1]][, 2], sphereRet[[1]][, 3],radius=0.03,col=sphereRet[[3]])

#difference
differentGeodesics <- (bunnyRet[[2]] - sphereRet[[2]])
differentGeodesics <- differentGeodesics + max(differentGeodesics)
differentGeodesicsColor <- ceiling(differentGeodesics[vertexNum,] * 15)+1

open3d()
plot3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3], xlab = "bunny geodesics sphere geodesics difference")
spheres3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3],radius=0.03,col=differentGeodesicsColor)

#bunny vertex, sphere geodesics
open3d()
plot3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3], xlab = "bunny vertex, sphere geodesics")
spheres3d(bunnyRet[[1]][, 1], bunnyRet[[1]][, 2], bunnyRet[[1]][, 3],radius=0.03,col=sphereRet[[3]])

#sphere vertex, bunny geodesics
open3d()
plot3d(sphereRet[[1]][, 1], sphereRet[[1]][, 2], sphereRet[[1]][, 3], xlab = "sphere vertex, bunny geodesics")
spheres3d(sphereRet[[1]][, 1], sphereRet[[1]][, 2], sphereRet[[1]][, 3],radius=0.03,col=bunnyRet[[3]])

