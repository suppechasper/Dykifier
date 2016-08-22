process.data <- function(infile = "data.mat", outfile="intersections.rda"){

  library(R.matlab)
  library(FNN)
  X <- readMat(infile)

  xlb <- -Inf #440000
  xub <- Inf #490000
  ylb <- -Inf #4990000
  yub <- Inf #5040000
  fname <- outfile

  dims <- dim(X$D)


  segments <- list()
  means <- matrix( NA, ncol=2, nrow=dims[2] )
  dirs <- matrix( NA, ncol=2, nrow=dims[2] )
  sdevs <- matrix( NA, ncol=2, nrow=dims[2] )
  slength <- rep(0, dims[2] )
  rtype <- c()
  #impl <- matrix(NA, ncol=3, nrow=dims[2])
  index <- 1
  for(i in 1:dims[2] ){
    if( min(X$D[,i,1]$X, na.rm=TRUE ) > xlb & 
       max(X$D[,i,1]$X, na.rm=TRUE ) < xub & 
       min(X$D[,i,1]$Y, na.rm=TRUE ) > ylb & 
       max(X$D[,i,1]$Y, na.rm=TRUE ) < yub ){ 
      segments[[index]] <- cbind( t(X$D[,i,1]$X), t(X$D[,i,1]$Y) )
      segments[[index]] <- segments[[index]][ complete.cases(segments[[index]]), ]
      pca <- prcomp(segments[[index]])
      means[i,] <- pca$center
      dirs[i, ] <- pca$rotation[,1]

      #segm = t(segments[[index]]) - means[i, ]
      #l <- dirs[i, ] %*% segm
      #sdevs[i, ] <- max(l) - min(l)
      sdevs[i, ] <- pca$sdev
      rtype[index] = X$D[,i,1]$host.rock

      for( j in 2:nrow( segments[[index]] ) ){
        slength[i] = slength[i] + sqrt( sum( (segments[[index]][j, ] 
                                            - segments[[index]][j-1, ] )^2 ) )

      }
      
      index <- index + 1
    }
    #impl[i, ] = c( pca$rotation[,2], -pca$rotation[,2] %*% means[i, ] )
  }

  means <- means[complete.cases(means), ]
  dirs <- dirs[complete.cases(dirs), ]
  sdevs <- sdevs[complete.cases(sdevs), ]
  slength <- slength[ !is.na(slength) ]

  ns <- length(segments)



  c1 <- matrix(NA, ncol=2, nrow=ns)
  c2 <- matrix(NA, ncol=2, nrow=ns )
  for(i in 1:ns ){
    c1[i, ] <- segments[[i]][1, ]
    c2[i, ] <- segments[[i]][ nrow(segments[[i]]), ]
  }

  k <- 100
  coords <- rbind(c1, c2) 
  KNN <- get.knn( coords, k)

  angles1 <- matrix(NA, ncol=k, nrow=nrow(KNN$nn.dist) )
  angles2 <- matrix(NA, ncol=k, nrow=nrow(KNN$nn.dist) )
  dists   <- matrix(NA, ncol=k, nrow=nrow(KNN$nn.dist) )

  for(i in 1:nrow(dists)  ){
    angle <- rep(NA, k) 
    dist <- rep(NA, k) 
    for(j in 1:k ){
      si <- i %% ns
      if(si == 0){
        si = ns
      }
      sj <- KNN$nn.index[i, j] %% ns
      if(sj == 0){
        sj = ns
      }
      if(si == sj){
        angles1[i, j] = Inf
        angles2[i, j] = Inf
        dists[i, j] = Inf
      }
      else{
        v1 <- means[si,] - coords[i,] 
        v2 <- coords[i,] - coords[ KNN$nn.index[i, j], ]
        v3 <- coords[ KNN$nn.index[i, j], ] - means[sj, ]
        v1 <- v1 / sqrt( sum( v1^2) )
        v2 <- v2 / sqrt( sum( v2^2) )
        v3 <- v3 / sqrt( sum( v3^2) )
        angles1[i, j] <- v1 %*% v2
        angles2[i, j] <- v2 %*% v3

        dists[i, j] <- KNN$nn.dist[i,j]
      }
    }
  }

  save( segments, means, dirs, sdevs, KNN, dists, angles1, angles2, coords,
       slength, file=fname )

}
