library(rgl) # package for 3d object handling
# reading the bunny.obj in your fair zip
bunny <- readOBJ("bunny.obj")

library(igraph) # package for graph theory

# 3D coordinates of vertices in the shape of n x 3 matrix
V.xyz <- t(bunny[[1]][1:3,])

# Enumerate edges of triangle in n x 2 matrix shape
Edges <- rbind(t(bunny[[2]][1:2,]),t(bunny[[2]][2:3,]),t(bunny[[2]][c(3,1),]))
# Remove duplicates of edges
Edges <- unique(Edges)
# length of edges
Edge.len <- sqrt(apply((V.xyz[Edges[,1],] - V.xyz[Edges[,2],])^2,1,sum))
# make a graph object
g <- graph.edgelist(Edges)
# distance on the graph
d <- distances(g,weights=Edge.len)

### Post-spherization
# Spherization maps all the vertices on the bunny on the S2 sphere
# Along the geodesics on the S2 sphere, we can "re-measure" the distance.
# The geodesics on the S2 sphere is drawn back on the bunny and
# along the geodesics, we can measure the "distance between vertices on the bunny".
# We can make the distance matrix of all vertices pairs along this back-to-the-bunny geodesics.
# Can you compare the two distance matrices; one is the distance matrix that the R codes above generates and the other is the back-to-the-bunny distance matrix.
# Showing distance from a point.

d. <- d/max(d)
# Coloring the vertices with distance from the n-th vertex
n <- 1
col <- rgb(d.[n,],1-d.[n,],1)
col2 <- ceiling(d.[n,] * 15)+1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col2)

# When you get different distance matrix
# you can replace the object d to the new distance matrix
# and draw the similar contours.
# Then, you can visually compare two distance matrices.


n <- 1
col <- rgb(d.[n,],1-d.[n,],1)
col2 <- ceiling(d.[n,] * 15)+1
col3 <- col2 %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3)

n2 <- 4000
col_ <- rgb(d.[n2,],1-d.[n2,],1)
col2_ <- ceiling(d.[n2,] * 15)+1
col3_ <- col2_ %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3_)

plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=as.numeric(col3 == col3_)+1)
