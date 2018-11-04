my.cell.sim.2 <- function (n = 5, k = 5, n.mesh = 32, n.step = 50, scale.shift = 0.1, 
    scale.rotation1 = 0.1, scale.rotation2 = 0.1, scale.shear = 0.1, 
    file.name = "hoge") 
{
    A. <- matrix(runif(n^2), n, n)
    A.[1, 1] <- k
    B <- matrix(rnorm(n^2), n, n)
    xxx <- my.spherical.harm.mesh(A = A., B = B, n = n.mesh)
    #plot3d(xxx$v)
    #segments3d(xxx$v[c(t(xxx$edge)), ])
    M <- diag(rep(1, 4))
    M[1:3, 4] <- runif(3) * scale.shift
    theta1 <- runif(1) * scale.rotation1
    theta2 <- runif(1) * scale.rotation2
    Mm <- diag(rep(1, 3))
    Mm[1:2, 1:2] <- my.2d.rot(theta1) %*% Mm[1:2, 1:2]
    Mm[2:3, 2:3] <- my.2d.rot(theta2) %*% Mm[2:3, 2:3]
    M[1:3, 1:3] <- Mm
    M[4, 1:3] <- rnorm(3) * scale.shear
    for (i in 1:n.step) {
        A. <- A. + rnorm(n^2, 0, 0.05)
        xxx <- my.spherical.harm.mesh(A = A., n = n.mesh)
        xxxx <- cbind(xxx$v, rep(1, length(xxx$v[, 1])))
        rot.xxxx <- t(M %*% t(xxxx))
        rot.xxxx. <- rot.xxxx[, 1:3]/rot.xxxx[, 4]
        #plot3d(rot.xxxx.)
        #segments3d(rot.xxxx.[c(t(xxx$edge)), ])
        file.out <- paste(file.name, i, ".obj", sep = "")
        if(i == n.step){
					my.write.obj(rot.xxxx., xxx$f, file.out)
				}
        
    }
    file.out.param <- paste(file.name, "_param.txt", sep = "")
    fileConn <- file(file.out.param, "w")
    cat(paste("n=", n, "\\n", sep = ""), file = fileConn)
    cat(paste("k=", k, "\\n", sep = ""), file = fileConn)
    cat(paste("n.mesh=", n.mesh, "\\n", sep = ""), file = fileConn)
    cat(paste("n.step=", n.step, "\\n", sep = ""), file = fileConn)
    cat(paste("scale.shift=", scale.shift, "\\n", sep = ""), file = fileConn)
    cat(paste("scale.rotation1=", scale.rotation1, "\\n", sep = ""), 
        file = fileConn)
    cat(paste("scale.rotation2=", scale.rotation2, "\\n", sep = ""), 
        file = fileConn)
    cat(paste("scale.shear=", scale.shear, "\\n", sep = ""), file = fileConn)
    close(fileConn)
}