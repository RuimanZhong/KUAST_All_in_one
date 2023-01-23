latt_generation <- function (xlim, ylim, by, mu, nu, scl, sig2, seed = NULL) {
  if (!is.null(seed)) {
    print('no seed')
  }
  
  x <- seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)
  y <- seq(from = ylim[1] + (by / 2), to = ylim[2] - (by / 2), by = by)
  coord <- expand.grid(x = x, y = y)
  
  # Define the covariance structure
  cov <- RMmatern(nu = nu, var = sig2, scale = scl) + RMtrend(mean = mu)
  
  # Generate data
  process <- RFsimulate(model = cov, x = coord)
  process <- data.frame(x = process@coords[, 1], y = process@coords[, 2], z = process$variable1)
  
  # Create a gridded spatial object from "process"
  coordinates(process) <- ~ x + y
  gridded(process) <- TRUE
  
  raster(process)
}

punifsample <- function(n,...){
  dp <- rpoint(n)
  #plot(dp)
  # Extract data at points
  s <- as.data.frame(cbind(x = dp$x,y = dp$y, value = raster::extract(r, cbind(dp$x, dp$y))))
  depoint = s %>% st_as_sf(coords = c("x", "y"), dim = "XY")  %>% st_cast("POINT")
  
  return(depoint)
}
pimsample <- function(n,sig.err, ...){
  im_r <- as.im(r)
  # sample from im_r based on value of im_r
  dp <- rpoint(n,im_r)
  #plot(dp)
  # Extract data at points
  s <- as.data.frame(cbind(x = dp$x,y = dp$y, value = raster::extract(r, cbind(dp$x, dp$y))))
  depoint = s %>% st_as_sf(coords = c("x", "y"), dim = "XY")  %>% st_cast("POINT")
  
  return(depoint)
}

pcoxsample <- function(n = NULL,rs, ...){
 
  lambda <- rs
  values(lambda) <- exp(beta0 + beta1* values(rs))
   
    if (!is.null(n)){
      loct <- rpoint(n = n, f = as.im(lambda), win = owin(xrange = xlim, yrange = ylim))
    
  } else {
    
    intensity_func <- function (x, y, intensity, ...) {
      raster::extract(intensity, cbind(x, y))
    }
    
    loct <- rpoispp(lambda = intensity_func, intensity = lambda)
    }
 
  loct
}

pgcsample <- function(n = NULL, rs,...){
  x <- seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)
  y <- seq(from = ylim[1] + (by / 2), to = ylim[2] - (by / 2), by = by)

  sim1 <- simgc(locs = cbind(x,y), sim.n = n, marginal = gm.gc(shape = 1, rate = 1),
                corr = matern.gc(range = scl0, kappa = nu, nugget = 0))
}

datagenerator <- function(loct,r, ...){
  Y = NULL
  n_points <- loct$n
  
  for (i in (1:n_points)) {
    Y[i] <- rnorm(n = 1, 
                  mean = raster::extract(r, rbind(c(loct$x[i], loct$y[i]))),
                  sd = sqrt(sig.err))
  }
  return(as.data.frame(cbind(value = Y, x= loct$x, y = loct$y)))
}
# Generate areas
areasample <- function(n, r){
  dearea <- aggregate(r, 50/n)
  #plot(dearea)
  # convert raster 2 to sp to be able to use extract()
  dearea <- rasterToPolygons(dearea)
  # Extract data at areas
  # sum values raster 1 in raster 2 with weights proportional to the overlapping area
  re <- raster::extract(r, dearea, weights = TRUE, normalizeWeights = TRUE)
  # this returns values and weights.
  # I cannot use extract() with fun = mean when weights = TRUE so I do it with values*weights
  dearea$value <- sapply(re, function(m){sum(m[, 1]*m[, 2], na.rm = TRUE)})
  
  dearea = st_as_sf(dearea)[,c(2,3)]
  
  return(dearea)
}

